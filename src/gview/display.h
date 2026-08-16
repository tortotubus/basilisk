#include <stddef.h>
#include <stdbool.h>

typedef struct {
  float *position;
  float *normal;
  float *color;
  unsigned int *index, *lindex;
  size_t position_len;
  size_t normal_len;
  size_t color_len;
  size_t index_len, lindex_len;
  int primitive;
} DisplayFrame;

int  display_start (void);
void display_publish (DisplayFrame *frame);
void display_stop (void);
void display_wait (void);

typedef struct {
  bool lighting;
  bool rotation;
  bool wait;
} GviewParameters;

extern GviewParameters Gview;
