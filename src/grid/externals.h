typedef struct _External External;

struct _External {
  char * name;    // the name of the variable
  void * pointer; // a pointer to the data
  int type;       // the type of the variable
  int nd;         // the number of pointer dereferences or attribute offset or enum constant
  char reduct;    // the reduction operation
  char global;    // is it a global variable?
  char constant;  // is it a constant?
  void * data;    // the dimensions (int *) for arrays or the code (char *) for functions
  scalar s;       // used for reductions on GPUs
  External * externals, * next;
  int used;
};

#if _GPU
static bool is_external_constant (const External * g) {
  return g->constant && g->type == sym_INT && !g->data;
}

static bool is_normal_variable (const External * g)
{
  if (g->name[0] == '.') return false;
  if (g->type == sym_function_declaration || g->type == sym_function_definition) return false;
  if (g->type == sym_INT && g->global &&
      (!strcmp (g->name, "N") ||
       !strcmp (g->name, "nl") ||
       !strcmp (g->name, "bc_period_x") ||
       !strcmp (g->name, "bc_period_y")))
    return false;
  return true;
}

static bool is_external_variable (const External * g)
{
  if (is_external_constant (g) || !is_normal_variable (g))
    return false;
  if (g->type == sym_INT ||
      g->type == sym_LONG ||
      g->type == sym_FLOAT ||
      g->type == sym_DOUBLE ||
      g->type == sym__COORD ||
      g->type == sym_COORD ||
      g->type == sym_BOOL ||
      g->type == sym_VEC4) {
    assert (!g->nd);
    assert (!g->data || ((int *)g->data)[1] == 0);
    return true;
  }
  return false;
}
#endif // _GPU
