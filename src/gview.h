#include "view.h"
#include "vertexbuffer.h"
#include "gview/display.h"
#pragma autolink -L$BASILISK/gview -lgview -lglfw

static DisplayFrame * basilisk_display_frame (void)
{
  DisplayFrame *frame;
  Array *position = VertexBuffer.position;
  Array *normal = VertexBuffer.normal;
  Array *color = VertexBuffer.color;
  Array *index = VertexBuffer.index;

  frame = malloc (sizeof(DisplayFrame));
  if (!frame)
    return NULL;

  frame->position = position->p;
  frame->normal = normal->p;
  frame->color = color->p;
  frame->index = index->p;

  frame->position_len = position->len;
  frame->normal_len = normal->len;
  frame->color_len = color->len;
  frame->index_len = index->len;
  frame->primitive = (VertexBuffer.type != 0);

  position->p = NULL;
  normal->p = NULL;
  color->p = NULL;
  index->p = NULL;

  array_free (position);
  array_free (normal);
  array_free (color);
  array_free (index);

  VertexBuffer.position = NULL;
  VertexBuffer.normal = NULL;
  VertexBuffer.color = NULL;
  VertexBuffer.index = NULL;

  return frame;
}

static void basilisk_display_publish (void)
{
  DisplayFrame *frame = basilisk_display_frame();
  if (frame)
    display_publish (frame);
}

/**
The Basilisk macro wraps the drawing commands.

The contents of the macro body are the commands which
generate the vertex buffer. */

macro gview()
{
  vertex_buffer_setup();
  {...}
  basilisk_display_publish();
}

event display_init (i = 0)
{
  display_start();
}

event display_cleanup (t = end)
{
  display_stop();
}
