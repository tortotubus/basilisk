#include "view.h"
#include "vertexbuffer.h"
#include "gview/display.h"
#include "gview/gui/ui.h"
#pragma autolink -L$BASILISK/gview -lgview -lglfw -limgui -lstb -lstdc++

static DisplayFrame * basilisk_display_frame (void)
{
  DisplayFrame *frame;
  Array *position = VertexBuffer.position;
  Array *normal = VertexBuffer.normal;
  Array *color = VertexBuffer.color;
  Array *index = VertexBuffer.index;
  Array *lindex = VertexBuffer.lindex;

  frame = malloc (sizeof(DisplayFrame));
  if (!frame)
    return NULL;

  frame->position = position->p;
  frame->normal = normal->p;
  frame->color = color->p;
  frame->index = index->p;
  frame->lindex = lindex->p;

  frame->position_len = position->len;
  frame->normal_len = normal->len;
  frame->color_len = color->len;
  frame->index_len = index->len;
  frame->lindex_len = lindex->len;
  frame->primitive = (VertexBuffer.type != 0);

  position->p = NULL;
  normal->p = NULL;
  color->p = NULL;
  index->p = NULL;
  lindex->p = NULL;

  array_free (position);
  array_free (normal);
  array_free (color);
  array_free (index);
  array_free (lindex);

  VertexBuffer.position = NULL;
  VertexBuffer.normal = NULL;
  VertexBuffer.color = NULL;
  VertexBuffer.index = NULL;
  VertexBuffer.lindex = NULL;

  return frame;
}

static void basilisk_display_publish (void)
{
  DisplayFrame *frame = basilisk_display_frame();
  if (frame)
    display_publish (frame);
}

event display_init (i = 0)
{
  Gview.rotation = Gview.lighting = dimension > 2;
  display_start();
}

static void gview_redraw()
{
  const char * commands = ui_commands();
  if (commands) {
    vertex_buffer_setup();
    char * copy = strdup (commands);
    char * line = copy, * end;
    int l = 1;
    while ((end = strchr (line, '\n'))) {
      *end = '\0';
      fprintf (stderr, "%d: %s\n", l++, line);
      process_line (line);
      line = end + 1;
    }
    free (copy);
    basilisk_display_publish();
  }
}

event display_ui_update (i++)
{
  gview_redraw(); // fixme: should limit the refresh rate at 50Hz
}

event display_cleanup (t = end)
{
  if (Gview.wait) {
    DisplayStatus redraw;
    while ((redraw = display_redraw()) != DISPLAY_STOP) {
      if (redraw == DISPLAY_REDRAW)
        gview_redraw();
      usleep (20000); // 50 Hz max
    }
  }
  display_stop();
}
