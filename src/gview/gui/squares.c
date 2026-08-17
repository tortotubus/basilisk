/**
# UI for squares() */

typedef struct {
  UIDrawItem ui;

  char color[64];
  double min, max, spread;
  bool linear;
} UISquares;

static void squares_properties (UIDrawItem *ui)
{
  UISquares *squares = (UISquares *) ui;

  IgInputText ("color", squares->color, sizeof(squares->color), 0);
  IgInputDouble ("min", &squares->min);
  IgInputDouble ("max", &squares->max);
  IgInputDouble ("spread", &squares->spread);
  IgCheckbox ("linear", &squares->linear);
}

static char *squares_command (UIDrawItem *ui)
{
  UISquares *squares = (UISquares *) ui;
  char *command;

  if (asprintf (&command,
                "squares (color = \"%s\", min = %g, max = %g, spread = %g, linear = %d);",
                squares->color, squares->min, squares->max, squares->spread, squares->linear))
    return NULL;

  return command;
}

static const char *squares_summary (UIDrawItem *ui)
{
  UISquares *squares = (UISquares *) ui;

  return squares->color;
}

static
UIDrawItem *ui_squares_new (unsigned int id)
{
  UISquares *squares = calloc (1, sizeof(UISquares));

  squares->ui.id = id;
  squares->ui.name = "squares";
  squares->ui.visible = true;
  squares->ui.properties = squares_properties;
  squares->ui.summary = squares_summary;
  squares->ui.command = squares_command;

  squares->linear = false;

  return &squares->ui;
}
