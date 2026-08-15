#include <stddef.h>
#include <stdbool.h>

typedef struct {
  float *position;
  float *normal;
  float *color;
  unsigned int *index;
  size_t position_len;
  size_t normal_len;
  size_t color_len;
  size_t index_len;
  int primitive;
} DisplayFrame;

int  display_start (void);
void display_publish (DisplayFrame *frame);
void display_stop (void);

typedef struct {
  bool lighting;
  bool rotation;  
} GviewParameters;

extern GviewParameters Gview;
