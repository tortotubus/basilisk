#include "gview/display.h"
#include <pthread.h>

static pthread_t gview_thread;

static void * gview_main (void *arg)
{
  // display_set_debug (true);
  if (display_main ("local", 0, NULL))
    exit (1);
  return NULL;
}

static void gview_terminate()
{
  if (pid() == 0) {
    display_quit();
    pthread_join (gview_thread, NULL);
  }
}

@init_solver
void gview_init()
{
  if (pid() == 0 &&
      pthread_create (&gview_thread, NULL, gview_main, NULL) != 0) {
    perror ("pthread_create");
      exit (1);
  }
  free_solver_func_add (gview_terminate);
}

#define DISPLAY_JS ""
#pragma autolink $BASILISK/gview/transfer_local.o -L$BASILISK/gview -ldisplay \
  -Lgl -lglutils -lglfw -limgui -lstb -lstdc++ -lm
#define DISPLAY -1
#include "display.h"
