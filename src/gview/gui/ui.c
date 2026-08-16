#include "bindings/dcimgui.h"
#include "bindings/dcimgui_impl_glfw.h"
#include "bindings/dcimgui_impl_opengl3.h"

/**
# Configuration */

static const float font_size = 18.f;
static const float items_width = 260.f;
static const float ui_top = 20.f;

/**
# Initialisation and termination */

static void ui_init (GLFWwindow *window)
{
  IgCreateContext (NULL);
  ImGuiIO *io = IgGetIO();
  ImFontAtlas_ClearFonts (io->Fonts);
  ImFontAtlas_AddFontFromFileTTF (io->Fonts,
                                  "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
                                  font_size, NULL, NULL);
  cImGui_ImplGlfw_InitForOpenGL (window, true);
  cImGui_ImplOpenGL3_Init();
}

static void ui_destroy()
{
  cImGui_ImplOpenGL3_Shutdown();
  cImGui_ImplGlfw_Shutdown();
  IgDestroyContext (NULL);
}

/**
# Drawing items */

#include "bindings/dcimgui.h"
#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  unsigned int id;
  char *name;
  char *parameters;
  bool visible;
} DisplayDrawItem;

static DisplayDrawItem *draw_items;
static int ndraw_items;
static unsigned int next_draw_item_id = 1;
static unsigned int selected_draw_id;

static void ui_add_draw_item (const char *name)
{
  DisplayDrawItem item = {
    .id = next_draw_item_id++,
    .name = strdup (name),
    .parameters = strdup (""),
    .visible = true
  };

  draw_items = realloc (draw_items,
                         (ndraw_items + 1)*sizeof(DisplayDrawItem));
  draw_items[ndraw_items++] = item;
}

static int ui_find_draw_item (unsigned int id)
{
  for (int i = 0; i < ndraw_items; i++)
    if (draw_items[i].id == id)
      return i;

  return -1;
}

static void ui_delete_draw_item (int index)
{
  free (draw_items[index].name);
  free (draw_items[index].parameters);

  memmove (&draw_items[index], &draw_items[index + 1],
           (ndraw_items - index - 1)*sizeof(DisplayDrawItem));
  ndraw_items--;

  if (ndraw_items == 0) {
    free (draw_items);
    draw_items = NULL;
  }

  selected_draw_id = 0;
}

static void ui_draw_menu (void)
{
  if (IgBeginMenu ("Draw")) {
    if (IgMenuItem ("squares"))
      ui_add_draw_item ("squares");

    if (IgMenuItem ("draw_vof"))
      ui_add_draw_item ("draw_vof");

    if (IgMenuItem ("isoline"))
      ui_add_draw_item ("isoline");

    IgEndMenu ();
  }
}

static void ui_draw_list (int width, int height)
{
  IgSetNextWindowPos ((ImVec2){width - items_width, ui_top},
                      ImGuiCond_Always);
  IgSetNextWindowSize ((ImVec2){items_width, height - ui_top},
                       ImGuiCond_Always);

  if (IgBegin ("Draw", NULL,
               ImGuiWindowFlags_NoMove |
               ImGuiWindowFlags_NoResize |
               ImGuiWindowFlags_NoCollapse)) {

    for (int i = 0; i < ndraw_items; i++) {
      bool selected = draw_items[i].id == selected_draw_id;

      IgPushIDInt (draw_items[i].id);

      IgCheckbox ("##visible", &draw_items[i].visible);
      IgSameLine ();

      if (IgSelectableEx (draw_items[i].name, selected, 0,
                          (ImVec2){0, 0}))
        selected_draw_id = draw_items[i].id;

      if (IgBeginPopupContextItem ()) {
        if (IgMenuItem ("Properties"))
          selected_draw_id = draw_items[i].id;

        if (IgMenuItem ("Delete")) {
          ui_delete_draw_item (i);
          IgEndPopup ();
          IgPopID ();
          break;
        }

        IgEndPopup ();
      }
  
      IgPopID ();
    }

    if (selected_draw_id && IgIsKeyPressed (ImGuiKey_Delete)) {
      int index = ui_find_draw_item (selected_draw_id);
      if (index >= 0)
        ui_delete_draw_item (index);
    }
  }

  IgEnd ();
}

/**
# Draw the user interface.

The UI belongs exclusively to the display thread. */

#include "bindings/dcimgui.h"
#include <GLFW/glfw3.h>

/**
Draw the user interface. */

static void ui (GLFWwindow *window)
{
  cImGui_ImplOpenGL3_NewFrame();
  cImGui_ImplGlfw_NewFrame();
  IgNewFrame();

  int width, height;
  
  glfwGetWindowSize (window, &width, &height);

  if (IgBeginMainMenuBar ()) {
    if (IgBeginMenu ("File")) {
      if (IgMenuItem ("Quit"))
        glfwSetWindowShouldClose (window, GLFW_TRUE);
      IgEndMenu ();
    }

    if (IgBeginMenu ("View")) {
      IgMenuItem ("Reset view");
      IgEndMenu ();
    }

    ui_draw_menu();

    IgEndMainMenuBar ();
  }

  ui_draw_list (width, height);

  IgRender ();
  cImGui_ImplOpenGL3_RenderDrawData (IgGetDrawData ());
}
