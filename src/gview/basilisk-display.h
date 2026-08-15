#ifndef BASILISK_DISPLAY_H
#define BASILISK_DISPLAY_H

#include "vertexbuffer.h"
#include "display.h"


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

  frame->nvertex = VertexBuffer.nvertex;
  frame->nindex =
    index->len/sizeof(unsigned int);

  frame->primitive =
    VertexBuffer.type == 0 ?
    DISPLAY_LINES : DISPLAY_TRIANGLES;

  /*
   * Transfer ownership of the allocated buffers.
   *
   * array_free() must not free the buffers, so detach
   * the pointers first.
   */
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
  DisplayFrame *frame = basilisk_display_frame ();

  if (frame)
    display_publish (frame);
}


/*
 * The Basilisk macro wraps the drawing commands.
 *
 * The contents of the macro body are the commands which
 * generate the vertex buffer.
 */

macro display_thread ()
{
  vertex_buffer_setup();

  {
    ...
  }

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

#endif
