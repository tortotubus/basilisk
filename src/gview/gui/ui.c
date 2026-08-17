#include "bindings/dcimgui.h"
#include "bindings/dcimgui_impl_glfw.h"
#include "bindings/dcimgui_impl_opengl3.h"
#include "ui.h"

/**
# Configuration */

static const int window_width = 1024, window_height = 768;
static const float font_size = 18.f;
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
# Draw item */

typedef struct UIDrawItem UIDrawItem;

struct UIDrawItem {
  unsigned int id;
  const char *name;
  bool visible;

  bool (*properties) (UIDrawItem *);
  const char *(*summary) (UIDrawItem *);
  char *(*command) (UIDrawItem *);
};

#include "../../draw_ui.h"

typedef struct {
  UIDrawItem **items;
  int nitems;
} UIDrawList;

static UIDrawList draw_list;
static unsigned int next_draw_item_id = 1;
static unsigned int selected_draw_id;

static bool properties_open;

UIDrawItem *ui_squares_new (unsigned int id);

static void ui_add_draw_item (UIDrawItem *item)
{
  draw_list.items = realloc (draw_list.items,
                             (draw_list.nitems + 1)*
                             sizeof(UIDrawItem *));

  draw_list.items[draw_list.nitems++] = item;
}

static int ui_find_draw_item (unsigned int id)
{
  for (int i = 0; i < draw_list.nitems; i++)
    if (draw_list.items[i]->id == id)
      return i;

  return -1;
}

static void ui_delete_draw_item (int index)
{
  unsigned int id = draw_list.items[index]->id;

  if (selected_draw_id == id) {
    selected_draw_id = 0;
    properties_open = false;
  }

  free (draw_list.items[index]);

  memmove (&draw_list.items[index],
           &draw_list.items[index + 1],
           (draw_list.nitems - index - 1)*sizeof(UIDrawItem *));

  draw_list.nitems--;

  if (draw_list.nitems == 0) {
    free (draw_list.items);
    draw_list.items = NULL;
  }
}

static bool ui_draw_menu (void)
{
  bool changed = false;
  if (IgBeginMenu ("Draw")) {
    if (IgMenuItem ("squares")) {
      ui_add_draw_item (ui_squares_new (next_draw_item_id++));
      changed = true;
    }
    IgEndMenu ();
  }
  return changed;
}

static bool ui_draw_list (int width, int height)
{
  bool changed = false;

  static float items_width = 250;
  static float items_width_start;

  int lwidth = items_width;

  /* Draw window */
  IgSetNextWindowPos ((ImVec2){0, ui_top}, ImGuiCond_Always);
  IgSetNextWindowSize ((ImVec2){lwidth, height - ui_top},
                       ImGuiCond_Always);
  if (IgBegin ("Draw", NULL,
               ImGuiWindowFlags_NoMove |
               ImGuiWindowFlags_NoTitleBar |
               ImGuiWindowFlags_NoResize |
               ImGuiWindowFlags_NoCollapse)) {

    if (selected_draw_id && IgIsKeyPressed (ImGuiKey_Delete)) {
      int index = ui_find_draw_item (selected_draw_id);

      if (index >= 0) {
        ui_delete_draw_item (index);
        changed = true;
      }
    }

    if (IgBeginTable ("draw_items", 3,
                      ImGuiTableFlags_BordersInnerV |
                      ImGuiTableFlags_Resizable |
                      ImGuiTableFlags_SizingFixedFit)) {

      IgTableSetupColumnEx ("",
                            ImGuiTableColumnFlags_WidthFixed,
                            30.f, 0);

      IgTableSetupColumnEx ("Draw",
                            ImGuiTableColumnFlags_WidthFixed,
                            100.f, 0);

      IgTableSetupColumnEx ("Properties",
                            ImGuiTableColumnFlags_WidthStretch,
                            0.f, 0);

      for (int i = 0; i < draw_list.nitems; i++) {
        UIDrawItem *item = draw_list.items[i];
        bool selected = item->id == selected_draw_id;

        IgPushIDInt (item->id);

        IgTableNextRow ();

        if (selected)
          IgTableSetBgColor (ImGuiTableBgTarget_RowBg1,
                             IM_COL32 (80, 80, 120, 255), -1);
        
        /* Visibility */
        IgTableNextColumn ();
        if (IgCheckbox ("##visible", &item->visible))
          changed = true;
        
        /* Name */
        IgTableNextColumn ();
        IgPushStyleColor (ImGuiCol_HeaderHovered, IM_COL32 (0, 0, 0, 0));
        IgPushStyleColor (ImGuiCol_HeaderActive, IM_COL32 (0, 0, 0, 0));
        if (IgSelectableEx (item->name, false, 0,
                            (ImVec2){0, 0})) {
          selected_draw_id = item->id;
          properties_open = true;
        }
        IgPopStyleColorEx (2);
 
        /* Summary */
        IgTableNextColumn ();
        IgPushStyleColor (ImGuiCol_HeaderHovered, IM_COL32 (0, 0, 0, 0));
        IgPushStyleColor (ImGuiCol_HeaderActive, IM_COL32 (0, 0, 0, 0));
        if (IgSelectableEx (item->summary (item), false, 0,
                            (ImVec2){0, 0})) {
          selected_draw_id = item->id;
          properties_open = true;
        }
        IgPopStyleColorEx (2);
 
        IgPopID ();
      }

      IgEndTable ();
    }
  }

  IgEnd ();
  
  /**
  Draggable separator. */

  const int separator_width = 10;
  
  IgPushStyleVarImVec2 (ImGuiStyleVar_WindowMinSize, (ImVec2){0.f, 0.f});
  IgSetNextWindowPos ((ImVec2){lwidth - separator_width/2.f, ui_top},
                      ImGuiCond_Always);
  IgSetNextWindowSize ((ImVec2){separator_width, height - ui_top},
                       ImGuiCond_Always);
  ImVec2 old_padding = IgGetStyle()->WindowPadding;
  IgGetStyle()->WindowPadding = (ImVec2){0.f, 0.f};
  if (IgBegin ("##draw_separator", NULL,
               ImGuiWindowFlags_NoDecoration |
               ImGuiWindowFlags_NoMove |
               ImGuiWindowFlags_NoResize |
               ImGuiWindowFlags_NoTitleBar |
               ImGuiWindowFlags_NoScrollbar |
               ImGuiWindowFlags_NoCollapse | 
               // ImGuiWindowFlags_NoSavedSettings |
               ImGuiWindowFlags_NoBringToFrontOnFocus)) {

    IgButtonEx ("##separator", IgGetContentRegionAvail());

    if (IgIsItemActivated ())
      items_width_start = items_width;

    if (IgIsItemActive ()) {
      ImVec2 delta = IgGetMouseDragDelta (0, -1.f);

      items_width = items_width_start + delta.x;

      if (items_width < 150)
        items_width = 150;
      if (items_width > 600)
        items_width = 600;
    }

  }

  IgEnd();
  
  IgPopStyleVar();
  IgGetStyle()->WindowPadding = old_padding;

  return changed;
}

static bool ui_draw_properties (void)
{
  int index = ui_find_draw_item (selected_draw_id);

  if (index < 0 || !properties_open)
    return false;

  UIDrawItem *item = draw_list.items[index];
  bool changed = false;
  
  IgSetNextWindowSize ((ImVec2){400.f, 200.f}, ImGuiCond_FirstUseEver);
  if (IgBegin (item->name, &properties_open, 0))
    changed = item->properties (item);

  IgEnd ();

  return changed;
}

/**
# Draw the user interface.

The UI belongs exclusively to the display thread. */

#include "bindings/dcimgui.h"
#include <GLFW/glfw3.h>

static bool ui (GLFWwindow *window)
{
  bool changed = false;
  
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
    
    if (IgBeginMenu ("Edit")) {
      int index = ui_find_draw_item (selected_draw_id);
      bool have_selection = index >= 0;

      if (IgMenuItemEx ("Delete", NULL, false, have_selection)) {
        ui_delete_draw_item (index);
        changed = true;
      }

      IgEndMenu ();
    }
    
    if (IgBeginMenu ("View")) {
      IgMenuItem ("Reset view");
      IgEndMenu ();
    }

    if (ui_draw_menu())
      changed = true;

    IgEndMainMenuBar ();
  }

  if (ui_draw_list (width, height))
    changed = true;
  
  if (ui_draw_properties ())
    changed = true;
  
  IgRender ();
  cImGui_ImplOpenGL3_RenderDrawData (IgGetDrawData ());

  return changed;
}

/**
# Returns the drawing commands */

const char *ui_commands (void)
{
  static char *commands;
  static size_t size;

  size_t length = 0;

  if (!commands) {
    size = 1024;
    commands = malloc (size);
  }

  commands[0] = '\0';
  for (int i = 0; i < draw_list.nitems; i++) {
    UIDrawItem *item = draw_list.items[i];
    if (!item->visible)
      continue;
    char *command = item->command (item);
    if (!command)
      continue;
    size_t n = strlen (command);
    if (length + n + 2 > size) {
      while (length + n + 2 > size)
        size *= 2;
      commands = realloc (commands, size);
    }
    memcpy (&commands[length], command, n);
    length += n;
    commands[length++] = '\n';
    commands[length] = '\0';
    free (command);
  }

  return commands;
}
