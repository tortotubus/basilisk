#include <stddef.h>

typedef struct {
  float *position;
  float *normal;
  float *color;
  unsigned int *index;

  size_t nvertex;
  size_t nindex;

  int primitive;
} DisplayFrame;

int  display_start (void);
void display_publish (DisplayFrame *frame);
void display_stop (void);
