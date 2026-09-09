/**
# Gview

To install

~~~bash
sudo apt install libimgui-dev
~~~
*/

#include "transfer.h"
#include "../gl/trackball.h"

#include "../grid/gpu/glad.h"
#include <GLFW/glfw3.h>

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#ifdef __EMSCRIPTEN__
# include <emscripten.h>
# include <emscripten/html5.h>

static EM_BOOL resize_callback (int eventType,
                                const EmscriptenUiEvent *event,
                                void *userData)
{
  GLFWwindow *window = userData;
  double width, height;

  emscripten_get_element_css_size ("#canvas-container", &width, &height);
  glfwSetWindowSize (window, (int) width, (int) height);

  return EM_TRUE;
}
#endif // __EMSCRIPTEN__

#include "parse_json.h"
#include "display.h"

#if UI
# include "gui/ui.h"
#endif

/**
# Parameters */

const float mouse_zoom_sensitivity = 10.f; // smaller is more sensitive
const unsigned int reconnect_ms = 100; // try to reconnect every 0.1 sec
const GLfloat light_direction[3] = {0.f, 0.f, 50.f};

/**
We limit the size or number of messages in the queue to between 50
and 100MB and 15 and 30 messages. Note that this means that messages
larger than 100MB can still be received, but only one at a time. */

const size_t high_bytes = 100*1024*1024, low_bytes = 50*1024*1024;
const size_t high_messages = 30, low_messages = 15;

/**
# Display */

typedef struct {
  unsigned int id;
  char * command;
  int dimension, type;
  bool visible;
  GLuint vbo[4], len[4];
  float color[3];
} DisplayItem;

typedef struct {
  GLFWwindow *window;
  DisplayItem ** items;
  struct transfer *transfer;

  int play; // -2: no run/pause controls, -1: pause, 0: run, 1: step
  
  int dimension;

  bool debug;

  int connected; // -1: waiting  0: disconnected  1: connected
} DisplayState;

static DisplayState Display = { .debug = false };

void display_set_debug (bool debug) {
  Display.debug = debug;
}

void display_set_connected (int connected)
{
  Display.connected = connected;
#if UI
  ui_set_connected (connected);
#endif
}

void display_quit()
{
  glfwSetWindowShouldClose (Display.window, GLFW_TRUE);
}

#define debug(...) do { if (Display.debug) {                            \
      fprintf (stderr, "%s:%d:%s: ", __FILE__, __LINE__, __func__);     \
      fprintf (stderr, __VA_ARGS__), fflush(stderr);                    \
    }} while(0)

/**
OpenGL state belongs exclusively to the display thread. */

static GLuint vao;

static GLint lighting_loc, light_direction_loc;
static GLint modelview_loc, projection_loc;

static GLuint program;

/**
# Camera */

DisplayCamera Camera;

static void matrix_identity (float *m)
{
  memset (m, 0, 16*sizeof(float));

  m[0] = 1.f;
  m[5] = 1.f;
  m[10] = 1.f;
  m[15] = 1.f;
}

static void matrix_perspective (float *m,
                                float fov, float aspect,
                                float near, float far)
{
  float f = 1.f/tanf (fov*(float) M_PI/360.f);

  memset (m, 0, 16*sizeof(float));

  m[0] = f/aspect;
  m[5] = f;
  m[10] = (far + near)/(near - far);
  m[11] = -1.f;
  m[14] = 2.f*far*near/(near - far);
}

static void display_quaternion_matrix (float m[16], const float q[4])
{
  float x = q[0], y = q[1], z = q[2], w = q[3];

  matrix_identity (m);

  m[0] = 1.f - 2.f*(y*y + z*z);
  m[1] = 2.f*(x*y - z*w);
  m[2] = 2.f*(x*z + y*w);

  m[4] = 2.f*(x*y + z*w);
  m[5] = 1.f - 2.f*(x*x + z*z);
  m[6] = 2.f*(y*z - x*w);

  m[8] = 2.f*(x*z - y*w);
  m[9] = 2.f*(y*z + x*w);
  m[10] = 1.f - 2.f*(x*x + y*y);
}

static void display_zoom (double factor)
{
  Camera.tz *= (float) factor;

  if (Camera.tz > -Camera.near)
    Camera.tz = -Camera.near;
  else if (Camera.tz < -100.f)
    Camera.tz = -100.f;
}

static void display_scroll (GLFWwindow *window, double xoffset, double yoffset)
{
#if UI
  if (ui_want_capture_mouse())
    return;
#endif
  
  display_zoom (powf (1.1f, (float) -yoffset));
}

static void display_mouse_button (GLFWwindow *window, int button, int action, int mods)
{
#if UI
  if (action == GLFW_PRESS && ui_want_capture_mouse())
    return;
#endif

  if (Display.dimension > 2 &&
      button == GLFW_MOUSE_BUTTON_LEFT) {
    Camera.rotating = action == GLFW_PRESS;
    if (Camera.rotating)
      glfwGetCursorPos (window, &Camera.mouse_x, &Camera.mouse_y);
  }
  else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
    Camera.panning = action == GLFW_PRESS;
    if (Camera.panning)
      glfwGetCursorPos (window, &Camera.mouse_x, &Camera.mouse_y);
  }
  else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
    Camera.zooming = action == GLFW_PRESS;
    if (Camera.zooming)
      glfwGetCursorPos (window, &Camera.mouse_x, &Camera.mouse_y);
  }
}

static void display_cursor_pos (GLFWwindow *window, double x, double y)
{
#if UI
  if (ui_want_capture_mouse() &&
      !Camera.rotating && !Camera.panning && !Camera.zooming)
    return;
#endif

  if (Camera.rotating) {
    int width, height;
    glfwGetWindowSize (window, &width, &height);
    float dquat[4];    
    gl_trackball (dquat,
                  (2.0*Camera.mouse_x -            width)/width,
                  (          height - 2.0*Camera.mouse_y)/height,
                  (           2.0*x -              width)/width,
                  (          height -              2.0*y)/height);
    gl_add_quats (dquat, Camera.rotation, Camera.rotation);

    Camera.mouse_x = x;
    Camera.mouse_y = y;
  }
  else if (Camera.panning) {
    int width, height;
    glfwGetWindowSize (window, &width, &height);

    if (width && height) {
      double dx = x - Camera.mouse_x;
      double dy = y - Camera.mouse_y;
      
      Camera.tx += -Camera.tz*(float) dx/width;
      Camera.ty -= -Camera.tz*(float) dy/height;
      
      Camera.mouse_x = x;
      Camera.mouse_y = y;
    }
  }
  else if (Camera.zooming) {
    double dy = y - Camera.mouse_y;

    if (dy)
      display_zoom (powf (1.1f, (float) dy/mouse_zoom_sensitivity));

    Camera.mouse_y = y;
  }
}

static void display_modelview (float m[16])
{
  float r[16];

  display_quaternion_matrix (r, Camera.rotation);

  /**
  modelview = translation * rotation * scale
  
  OpenGL matrices are column-major, so scaling the columns
  of the rotation matrix gives rotation * scale. */

  m[0]  = Camera.sx*r[0];
  m[1]  = Camera.sx*r[1];
  m[2]  = Camera.sx*r[2];
  m[3]  = 0.f;

  m[4]  = Camera.sy*r[4];
  m[5]  = Camera.sy*r[5];
  m[6]  = Camera.sy*r[6];
  m[7]  = 0.f;

  m[8]  = Camera.sz*r[8];
  m[9]  = Camera.sz*r[9];
  m[10] = Camera.sz*r[10];
  m[11] = 0.f;

  m[12] = Camera.tx;
  m[13] = Camera.ty;
  m[14] = Camera.tz;
  m[15] = 1.f;
}

/**
# Shader compilation */

static GLuint compile_shader (GLenum type, const char *source)
{
  GLuint shader;
  GLint status;
  char log[4096];

  shader = glCreateShader (type);
  glShaderSource (shader, 1, &source, NULL);
  glCompileShader (shader);
  glGetShaderiv (shader, GL_COMPILE_STATUS, &status);
  if (!status) {
    glGetShaderInfoLog (shader, sizeof(log), NULL, log);
    fprintf (stderr,
             "display: shader compilation failed:\n%s\n",
             log);
    glDeleteShader (shader);
    return 0;
  }
  
  return shader;
}

static GLuint create_program (void)
{
  // Use OpenGL ES 3.0 / WebGL 2.0 shaders for Emscripten
#ifdef __EMSCRIPTEN__
  static const char *vertex_source =
    "#version 300 es\n"
    "precision highp float;\n"
    "\n"
    "in vec3 position;\n"
    "in vec3 normal;\n"
    "in vec3 color;\n"
#else // Use OpenGL >= 4.3 for the native version
  static const char *vertex_source =
    "#version 430\n"
    "\n"
    "layout (location = 0) in vec3 position;\n"
    "layout (location = 1) in vec3 normal;\n"
    "layout (location = 2) in vec3 color;\n"
#endif
    "\n"
    "uniform mat4 modelview;\n"
    "uniform mat4 projection;\n"
    "uniform bool lighting;\n"
    "\n"
    "out vec3 vertex_normal;\n"
    "out vec3 vertex_color;\n"
    "\n"
    "void main ()\n"
    "{\n"
    "  gl_Position = projection * modelview * vec4 (position, 1.);\n"
    "  if (lighting)\n"
    "    vertex_normal = mat3 (modelview)*normal;\n"
    "  vertex_color = color;\n"
    "}\n";
#ifdef __EMSCRIPTEN__
  static const char *fragment_source =
    "#version 300 es\n"
    "precision highp float;\n"
    "\n"
    "in vec3 vertex_normal;\n"
    "in vec3 vertex_color;\n"
#else // Native
  static const char *fragment_source =
    "#version 430\n"
    "\n"
    "in vec3 vertex_normal;\n"
    "in vec3 vertex_color;\n"
#endif
    "\n"
    "uniform bool lighting;\n"
    "uniform vec3 light_direction;\n"
    "\n"
    "out vec4 frag_color;\n"
    "\n"
    "void main ()\n"
    "{\n"
    "  float intensity = 1.;\n"
    "  if (lighting) {\n"
    "    vec3 n = normalize (vertex_normal);\n"
    "    if (!gl_FrontFacing)\n"
    "      n = -n;\n"
    "    float diffuse = max (dot (n, normalize (light_direction)), 0.);\n"
    "    intensity = 0.2 + 0.8*diffuse;\n"
    "  }\n"
    "  frag_color = vec4 (intensity*vertex_color, 1.);\n"
    "}\n";

  GLuint vertex;
  GLuint fragment;
  GLuint result;
  GLint status;
  char log[4096];

  vertex = compile_shader (GL_VERTEX_SHADER, vertex_source);
  fragment = compile_shader (GL_FRAGMENT_SHADER, fragment_source);

  if (!vertex || !fragment) {
    if (vertex)
      glDeleteShader (vertex);
    if (fragment)
      glDeleteShader (fragment);
    return 0;
  }

  result = glCreateProgram ();

  glAttachShader (result, vertex);
  glAttachShader (result, fragment);
  glLinkProgram (result);

  glGetProgramiv (result, GL_LINK_STATUS, &status);
  if (!status) {
    glGetProgramInfoLog (result, sizeof(log), NULL, log);
    fprintf (stderr, "display: shader linking failed:\n%s\n", log);
    glDeleteProgram (result);
    result = 0;
  }

  glDeleteShader (vertex);
  glDeleteShader (fragment);

  return result;
}

/**
# OpenGL initialization */

void display_front_view()
{
  Camera.rotation[0] = 0;
  Camera.rotation[1] = 0;
  Camera.rotation[2] = 0;
  Camera.rotation[3] = 1;
}

void display_back_view()
{
  Camera.rotation[0] = 0;
  Camera.rotation[1] = 1;
  Camera.rotation[2] = 0;
  Camera.rotation[3] = 0;
}

#define ISQ2 0.707106781186547f

void display_left_view()
{
  Camera.rotation[0] = 0;
  Camera.rotation[1] = -0.707107;
  Camera.rotation[2] = 0;
  Camera.rotation[3] = 0.707107;
}

void display_right_view()
{
  Camera.rotation[0] = 0;
  Camera.rotation[1] = 0.707107;
  Camera.rotation[2] = 0;
  Camera.rotation[3] = 0.707107;
}

void display_top_view()
{
  Camera.rotation[0] = -0.707107;
  Camera.rotation[1] = 0;
  Camera.rotation[2] = 0;
  Camera.rotation[3] = 0.707107;
}

void display_bottom_view()
{
  Camera.rotation[0] = 0.707107;
  Camera.rotation[1] = 0;
  Camera.rotation[2] = 0;
  Camera.rotation[3] = 0.707107;
}

void display_iso_view()
{
  Camera.rotation[0] = -0.353553;
  Camera.rotation[1] = 0.353553;
  Camera.rotation[2] = 0.146447;
  Camera.rotation[3] = 0.853553;
}

void display_reset_view()
{
  Camera.tx = Camera.ty = 0.f;
  Camera.tz = -3.f;

  Camera.sx = Camera.sy = Camera.sz = 1.f;

  Camera.rotation[0] = 0.f;
  Camera.rotation[1] = 0.f;
  Camera.rotation[2] = 0.f;
  Camera.rotation[3] = 1.f;

  Camera.fov = 24.f;
  Camera.near = 0.01f;
  Camera.far = 100.f;

  Camera.rotating = false;
  Camera.panning = false;
  Camera.zooming = false;
}

static void display_gl_init (void)
{
  program = create_program ();

  if (!program)
    exit (1);

  lighting_loc = glGetUniformLocation (program, "lighting");
  light_direction_loc = glGetUniformLocation (program, "light_direction");

  modelview_loc = glGetUniformLocation (program, "modelview");
  projection_loc = glGetUniformLocation (program, "projection");

  // OpenGL ES 3.0 / WebGL 2.0 compatible VAO setup
#ifdef __EMSCRIPTEN__
  // WebGL 2.0: Use traditional vertex attribute setup
  glGenVertexArrays (1, &vao);
  glBindVertexArray (vao);

  // Enable and set up vertex attributes
  glEnableVertexAttribArray (0);
  glEnableVertexAttribArray (1);
  glEnableVertexAttribArray (2);

  glBindVertexArray (0);
#else
  // OpenGL 4.3: Use DSA
  glGenVertexArrays (1, &vao);

  glBindVertexArray (vao);

  glEnableVertexAttribArray (0);
  glVertexAttribFormat (0, 3, GL_FLOAT, GL_FALSE, 0);
  glVertexAttribBinding (0, 0);

  glEnableVertexAttribArray (1);
  glVertexAttribFormat (1, 3, GL_FLOAT, GL_FALSE, 0);
  glVertexAttribBinding (1, 1);

  glEnableVertexAttribArray (2);
  glVertexAttribFormat (2, 3, GL_FLOAT, GL_FALSE, 0);
  glVertexAttribBinding (2, 2);

  glBindVertexArray (0);
#endif

  glEnable (GL_DEPTH_TEST);
#ifndef __EMSCRIPTEN__
  glEnable (GL_MULTISAMPLE);
#endif

  display_reset_view();
}

/**
Upload a newly received frame */
#if 0
static void display_upload (DisplayFrame *frame)
{
  glBindVertexArray (vao);

  glBindBuffer (GL_ARRAY_BUFFER, position_vbo);
  glBufferData (GL_ARRAY_BUFFER, frame->position_len,
                frame->position, GL_DYNAMIC_DRAW);

  glBindBuffer (GL_ARRAY_BUFFER, normal_vbo);
  glBufferData (GL_ARRAY_BUFFER, frame->normal_len,
                frame->normal, GL_DYNAMIC_DRAW);

  constant_normal_enabled = frame->normal_len == 3*sizeof(float);

  if (constant_normal_enabled)
    memcpy (constant_normal_value, frame->normal, 3*sizeof(float));
  else {
    glBindBuffer (GL_ARRAY_BUFFER, normal_vbo);
    glBufferData (GL_ARRAY_BUFFER, frame->normal_len,
                  frame->normal, GL_DYNAMIC_DRAW);
  }
 
  constant_color_enabled = frame->color_len == 3*sizeof(float);

  if (constant_color_enabled)
    memcpy (constant_color_value, frame->color, 3*sizeof(float));
  else {
    glBindBuffer (GL_ARRAY_BUFFER, color_vbo);
    glBufferData (GL_ARRAY_BUFFER, frame->color_len,
                  frame->color, GL_DYNAMIC_DRAW);
  }

  glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, index_vbo);
  glBufferData (GL_ELEMENT_ARRAY_BUFFER, frame->index_len,
                frame->index, GL_DYNAMIC_DRAW);
  index_count = frame->index_len/sizeof(unsigned int);
  
  glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, lindex_vbo);
  glBufferData (GL_ELEMENT_ARRAY_BUFFER, frame->lindex_len,
                frame->lindex, GL_DYNAMIC_DRAW);
  lindex_count = frame->lindex_len/sizeof(unsigned int);

  glBindVertexArray (0);
}
#endif

/**
Render the current VBO */

static void display_render (GLFWwindow *window)
{
  int xpos = 0, ypos = 0, width, height;

#if UI
  ui_viewport (window, &xpos, &ypos, &width, &height);
#else
  glfwGetFramebufferSize (window, &width, &height);
#endif
  
  if (!height)
    return;

  glViewport (xpos, ypos, width, height);

  float aspect = (float) width/height;
  float modelview[16];
  float projection[16];
  
  matrix_perspective (projection,
                      Camera.fov, aspect,
                      Camera.near, Camera.far);
  
  display_modelview (modelview);

  if (Camera.bg[3] == 0) {
    if (Display.dimension == 2)
      Camera.bg[0] = 1, Camera.bg[1] = 1, Camera.bg[2] = 1, Camera.bg[3] = 1;
    else
      Camera.bg[0] = 0.3, Camera.bg[1] = 0.4, Camera.bg[2] = 0.6, Camera.bg[3] = 1;
  }
  glClearColor (Camera.bg[0], Camera.bg[1], Camera.bg[2], Camera.bg[3]);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glUseProgram (program);

  glUniformMatrix4fv (modelview_loc, 1, GL_FALSE, modelview);
  glUniformMatrix4fv (projection_loc, 1, GL_FALSE, projection);
  glUniform3fv (light_direction_loc, 1, light_direction);

  glBindVertexArray (vao);

#define foreach_item(item)                                              \
  if (Display.items)                                                    \
    for (DisplayItem ** _i = Display.items, * item = *_i; item; _i++, item = *_i)

  int dimension = 0;

  /**
  We need to draw in reverse order to respect visibility somehow... */
  
  int nitems = 0;
  foreach_item (item) nitems++;
  for (int i = nitems - 1; i >= 0; i--) {
    DisplayItem * item = Display.items[i];
    
    if (item->vbo[0] && item->visible) {
#ifdef __EMSCRIPTEN__
      // OpenGL ES 3.0 / WebGL 2.0: Traditional vertex attribute setup
      glBindBuffer (GL_ARRAY_BUFFER, item->vbo[0]);
      glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE, 0, 0);
      glEnableVertexAttribArray (0);

      if (item->len[1] == item->len[0]) { // normals
        glBindBuffer (GL_ARRAY_BUFFER, item->vbo[1]);
        glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray (1);
      }
      else { // constant normal
        glDisableVertexAttribArray (1);
        float constant_normal[3] = {0,0,1};
        glVertexAttrib3fv (1, constant_normal);
      }

      if (item->len[2] == item->len[0]) { // colors
        glBindBuffer (GL_ARRAY_BUFFER, item->vbo[2]);
        glVertexAttribPointer (2, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray (2);
      }
      else { // constant color
        glDisableVertexAttribArray (2);
        glVertexAttrib3fv (2, item->color);
      }
#else
      // OpenGL 4.3: Use DSA
      glBindVertexBuffer (0, item->vbo[0], 0, 3*sizeof(float)); // vertices

      if (item->len[1] == item->len[0]) { // normals
        glEnableVertexAttribArray (1);
        glBindVertexBuffer (1, item->vbo[1], 0, 3*sizeof(float));
      }
      else { // constant normal
        glDisableVertexAttribArray (1);
        float constant_normal[3] = {0,0,1};
        glVertexAttrib3fv (1, constant_normal);
      }

      if (item->len[2] == item->len[0]) { // colors
        glEnableVertexAttribArray (2);
        glBindVertexBuffer (2, item->vbo[2], 0, 3*sizeof(float));
      }
      else { // constant color
        glDisableVertexAttribArray (2);
        glVertexAttrib3fv (2, item->color);
      }
#endif
      
      glUniform1i (lighting_loc, item->dimension > 2 && item->len[1] != 0);
 
      glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, item->vbo[3]);
        
      if (item->type == 0)
        glDrawElements (GL_LINES, (GLsizei) (item->len[3]/sizeof(unsigned int)),
                        GL_UNSIGNED_INT, NULL);
      else {
        glEnable (GL_POLYGON_OFFSET_FILL);
        glPolygonOffset (1.f, 1.f);
        glDrawElements (GL_TRIANGLES, (GLsizei) (item->len[3]/sizeof(unsigned int)),
                        GL_UNSIGNED_INT, NULL);
        glDisable (GL_POLYGON_OFFSET_FILL);
      }
      if (item->dimension > dimension)
        dimension = item->dimension;
    }
  }

  glBindVertexArray (0);

  if (dimension && dimension != Display.dimension) {
    debug ("dimension: %d %d\n", dimension, Display.dimension);
    bool reset_view = (Display.dimension > dimension);
    Display.dimension = dimension;
#if UI
    ui_set_dimension (dimension);
#endif
    if (reset_view)
      display_reset_view();
  }
}

/**
# Display items */

static void send_item (const char * command, struct transfer *transfer, char op)
{
  char * s;
  assert (asprintf (&s, "%c%s", op, command) > 0);
  transfer_send (transfer, s, strlen (s), false);
  free (s);
}

static void show_item_list()
{
  debug ("**************************\n");
  foreach_item (i)
    debug ("%d: %s\n", i->id, i->command);
  debug ("--------------------------\n");
}

static DisplayItem * add_item (const char * command, struct transfer *transfer)
{
  DisplayItem * item = calloc (1, sizeof (DisplayItem));
  item->command = strdup (command);
  item->visible = true;

  int len = 0;
  foreach_item (i) len++;
  
  Display.items = realloc (Display.items, (len + 2)*sizeof (DisplayItem *));
  Display.items[len] = item;
  item->id = len;
  Display.items[len + 1] = NULL;
  
  debug ("adding item %d:'%s'\n", item->id, item->command);
  
  send_item (item->command, transfer, '+');
  
  show_item_list();
  return item;
}

static void remove_item (DisplayItem * item, struct transfer *transfer)
{
  debug ("removing item '%s'\n", item->command);
  
  if (transfer)
    send_item (item->command, transfer, '-');

  int len = 0, found = -1;
  foreach_item (i) {
    if (i == item)
      found = len;
    len++;
  }
  assert (found >= 0);
  
  for (int i = found; i < len; i++) {
    Display.items[i] = Display.items[i+1];
    if (Display.items[i])
      Display.items[i]->id = i;
  }

  glDeleteBuffers (4, item->vbo);
  
  free (item->command);
  free (item);
  
  show_item_list();
}

#if UI
void display_add_item (const char * command)
{
  add_item (command, Display.transfer);
}

void display_remove_item (unsigned int id)
{
  if (Display.items && Display.items[id])
    remove_item (Display.items[id], Display.transfer);
}

void display_update_item (unsigned int id, const char * command)
{
  DisplayItem * item = Display.items[id];
  send_item (item->command, Display.transfer, '-');
  free (item->command);
  item->command = strdup (command);
  send_item (item->command, Display.transfer, '+');
}

void display_item_visible (unsigned int id, bool visible)
{
  DisplayItem * item = Display.items[id];
  if (item->visible) {
    assert (!visible);
    debug ("item '%s' is not visible\n", item->command);
    send_item (item->command, Display.transfer, '-');
  }
  else {
    assert (visible);
    debug ("item '%s' is visible\n", item->command);
    send_item (item->command, Display.transfer, '+');
  }
  item->visible = visible;
}
#endif

static void update_item (const char * command,
                         int dimension,
                         int type,
                         unsigned int len[4],
                         char * data)
{
  glBindVertexArray (vao);
  
  foreach_item (item) {
    if (!strcmp (item->command, command)) {
      item->dimension = dimension;
      item->type = type;
      if (!item->vbo[0])
        glGenBuffers (4, item->vbo);
      for (int i = 0; i < 3; i++) {
        glBindBuffer (GL_ARRAY_BUFFER, item->vbo[i]);
        glBufferData (GL_ARRAY_BUFFER, len[i], data, GL_DYNAMIC_DRAW);
        if (i == 2 && len[i] == 12) // constant color
          memcpy (item->color, data, len[i]);
        data += len[i];
      }
      glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, item->vbo[3]);
      glBufferData (GL_ELEMENT_ARRAY_BUFFER, len[3], data, GL_DYNAMIC_DRAW);
      for (int i = 0; i < 4; i++)
        item->len[i] = len[i];
#if UI
      ui_set_error (item->id, NULL);
#endif
    }
  }

  glBindVertexArray (0);
}

/**
A helper function. */

static
void sort_by_weight(char *lines[], int weight[], int nlines)
{
  for (int i = 0; i < nlines - 1; i++) {
    int min_idx = i;
    for (int j = i + 1; j < nlines; j++) {
      if (weight[j] < weight[min_idx])
        min_idx = j;
    }
    if (min_idx != i) {
      int tmp = weight[i];
      weight[i] = weight[min_idx];
      weight[min_idx] = tmp;

      char *tmp2 = lines[i];
      lines[i] = lines[min_idx];
      lines[min_idx] = tmp2;
    }
  }
}

/**
Display control */

static
const char * display_control_value (DisplayControlStrings control, const char * key)
{
  for (int i = 0; i < control.n; i++)
    if (!strcmp (control.keys[i], key))
      return control.values[i];
  return NULL;
}

void display_control_command (const char * command)
{
  transfer_send (Display.transfer, command, strlen (command), false);
}

/**
# Server messages */

#define READ_DATA(var, lsize)          \
  assert (size >= lsize);              \
  memcpy (var, data, lsize),           \
    data += lsize,                     \
    size -= lsize
#define READ_VAR(var) READ_DATA(&var, sizeof(var))

static int poll_messages (struct transfer *transfer)
{
  struct transfer_message *message = transfer_receive (transfer);
  if (!message)
    return 0;
    
  /**
  ## Connection closed */
    
  if (message->type == TRANSFER_MESSAGE_CLOSED) {
    debug ("Connection closed\n");
    display_set_connected (reconnect_ms > 0 ? -1 : 0);
      transfer_message_free (message);
      return -1;
  }

  /**
  ## Binary messages */
    
  if (message->binary) {
    if (!Display.connected) {
      fprintf (stderr, "%s:%d: unexpected binary message\n", __FILE__, __LINE__);
      transfer_message_free (message);
        return -1;
    }
    
    debug ("Received binary message: %zu bytes\n",
           message->size);
    size_t size = message->size;
    char * data = message->data;
    unsigned int commandlen; READ_VAR (commandlen);
    unsigned int errorlen; READ_VAR (errorlen);
    char command[commandlen + 1]; READ_DATA(command, commandlen);
    command[commandlen] = '\0';
    debug ("Received buffer for command '%s'\n", command);
      
    // padding to a multiple of four
    int paddedlen = 4*ceil(commandlen/4.);
    for (int i = 0; i < paddedlen - commandlen; i++) {
      assert (size >= 1);
      assert (*data == '\0');
      data++, size--;
    }

    /**
    ### Status */

    if (!strcmp (command, "status();")) {
      if (errorlen) {
        char status[errorlen + 1]; READ_DATA(status, errorlen);
        status[errorlen] = '\0';
#if UI
        debug ("    status: %s\n", status);
        ui_set_status (status);
#else
        fprintf (stderr, "status: %s           \r", status);
#endif
      }
      else {
#if UI
        debug ("    status: empty\n");
        ui_set_status (NULL);
#else
        fprintf (stderr, "status:                  \r");
#endif
      }
    }
    
    /**
    ### Errors */
      
    else if (errorlen) {
      char error[errorlen + 1]; READ_DATA(error, errorlen);
      error[errorlen] = '\0';
      debug ("    error: %s\n", error);
#if UI
      foreach_item (item)
        if (!strcmp (item->command, command))
          ui_set_error (item->id, error);
#else
      fprintf (stderr, "%s\nerror: %s\n", command, error);
#endif
    }

    /**
    ### VertexBuffer */
      
    else {
      int dimension; READ_VAR (dimension);
      int type; READ_VAR (type);
      unsigned int len[4]; READ_VAR(len);
      debug ("    dimension: %d type: %d pos: %d normal: %d color: %d index: %d\n",
             dimension, type, len[0], len[1], len[2], len[3]);
      assert (size == len[0] + len[1] + len[2] + len[3]);
      update_item (command, dimension, type, len, data);
    }
  }

  /**
  ## Text messages */
    
  else {
    debug ("Received: %.*s\n",
           (int) message->size,
           (char *) message->data);
    
    char nullterminated[message->size + 1], * msg = nullterminated;
    memcpy (nullterminated, message->data, message->size);
    nullterminated[message->size] = '\0';
      
    if (Display.connected <= 0 && msg[0] != '{') {
      fprintf (stderr, "%s:%d: error: expecting a connection message\n",
               __FILE__, __LINE__);
      transfer_message_free (message);
      return -1;
    }
    
    switch (msg[0]) {

    /**
    ### The initialisation message from the running Basilisk code */
        
    case '{': {

      /**
      A sanity check. */

      if (Display.connected > 0) {
        fprintf (stderr, "%s:%d: error: unexpected connection message\n",
                 __FILE__, __LINE__);
        transfer_message_free (message);
        return -1;
      }
      if (message->size < 5 || strncmp (msg, "{\n  \"", 5)) {
        fprintf (stderr, "%s:%d: error: invalid connection message\n", __FILE__, __LINE__);
        transfer_message_free (message);
        return -1;
      }
      display_set_connected (1);

      /**
      The server has just started and does not know about our items,
      so we send them. */

      foreach_item (item)
        if (item->visible)
          send_item (item->command, transfer, '+');

      /**
      We also ask for status updates. */

      send_item ("status();", transfer, '+');
      
      break;
    }
        
    /**
    ### The controls from the running Basilisk code

    This is JSON but with a simple format e.g.
    "Run/Pause": { "type": "int", "value": -1, "min": -1, "max": 1 },
    */

    case '#': {

      if (!strcmp (msg, "#{}")) { // no controls
        Display.play = -2;
#if UI
        ui_set_running (Display.play);
#endif
        break;
      }
      
      while (msg) {
        msg = strchr (msg, '"'); if (!msg) goto error;
        char * end = strchr (msg, '\n'); if (!end) end = msg + strlen(msg);
        *(end - 1) = '\0';

        char *name, *keys[16], *values[16];
        int n = parse_entry (msg, &name, keys, values, 16);
        if (n < 0) goto error;

        DisplayControlStrings control = { name, keys, values, n };

        const char * value;
        if (!strcmp (name, "Run/Pause") && n == 4 &&
            (value = display_control_value (control, "value"))) {            
          Display.play = atoi(value);
#if UI
          ui_set_running (Display.play);
#endif
        }
#if UI
        else
          ui_add_control (control);
#endif 

        if (*end == '\0')
          break;
        msg = end + 1;
      }
      break;
      error:
      fprintf (stderr, "%s:%d: error parsing controls\n%s\n",
               __FILE__, __LINE__, nullterminated);

      break;
    }
        
    /**
    ### Default display commands */

    case '@':
      if (Display.items && Display.items[0]) {
        // ignore default display if some geometries are already present
        transfer_message_free (message);
        return 1;
      }
      msg++; // skip @
      // fall through
        
    /**
    ## Commands */
      
    default: {
      char * lines[80];
      int orders[80];
      int nlines = 0;
      char * line = msg, * end;
      while ((end = strchr (line, ';'))) {
        if (nlines >= 80) {
          fprintf (stderr, "%s:%d: error: too many commands (> 80)\n",
                   __FILE__, __LINE__);
          break;
        }
        
        char c = *++end;
        *end = '\0';

        /**
        Replace ' with " */
        
        for (char *p = line; *p; p++)
          if (*p == '\'')
            *p = '"';

        /**
        'order' parameter. */

        int order = -1;
        char * sorder = strstr (line, "order");
        if (sorder) {
          if (strchr (" \t\n,(", *(sorder - 1))) {
            char * s = sorder + 5;
            while (*s != '\0' && strchr (" \t\n", *s)) s++;
            if (*s == '=') {
              s++;
              while (*s != '\0' && strchr (" \t\n", *s)) s++;
              char * value = s;
              while (*s != '\0' && *s >= '0' && *s <= '9') s++;
              if (strchr (" \t\n)", *s)) {
                char c = *s; *s = '\0';
                order = atoi (value);
                debug ("order: '%s' %d\n", sorder, order);
                for (char * s = sorder; *s != '\0'; s++) *s = ' ';
                *s = c;
              }
            }
          }
        }

        lines[nlines] = line;
        orders[nlines++] = order;
          
        *end = c;
        line = end;
      }

      /**
      Sort in increasing order. */

      sort_by_weight (lines, orders, nlines);
      
      for (int i = 0; i < nlines; i++) {
#if UI
        ui_add_item (lines[i]);
#else
        add_item (lines[i], transfer);
#endif
      }
    }
    }
  }
      
  transfer_message_free (message);
  
  return 1;
}

/**
# Keyboard callbacks */

void display_set_running (int play)
{
  if (Display.connected > 0) {
    Display.play = play;
    char * s;
    assert (asprintf (&s, "#Run/Pause:%d", Display.play) > 0);
    display_control_command (s);
    debug ("key_callback(): %s\n", s);
    free (s);
#if UI
    ui_set_running (Display.play);
#endif
  }
}

static void key_callback (GLFWwindow *window, int key, int scancode,
                          int action, int mods)
{
#if UI
  if (action == GLFW_PRESS && ui_want_capture_keyboard())
    return;
#endif
  
  switch (key) {
  case GLFW_KEY_SPACE:
    if (action == GLFW_PRESS)
      display_set_running (Display.play ? 0 : -1);
    break;

  case GLFW_KEY_S:
    if ((action == GLFW_PRESS || action == GLFW_REPEAT))
      display_set_running (1);
    break;

  case GLFW_KEY_F:
    if (action == GLFW_PRESS)
      display_front_view();
    break;
    
  case GLFW_KEY_B:
    if (Display.dimension > 2 && action == GLFW_PRESS) {
      if (mods & GLFW_MOD_CONTROL)
        display_bottom_view();
      else
        display_back_view();
    }
    break;
    
  case GLFW_KEY_L:
    if (Display.dimension > 2 && action == GLFW_PRESS)
      display_left_view();
    break;
    
  case GLFW_KEY_R:
    if (action == GLFW_PRESS) {
      if (mods & GLFW_MOD_CONTROL)
        display_reset_view();
      else if (Display.dimension > 2)
        display_right_view();
    }
    break;

  case GLFW_KEY_T:
    if (Display.dimension > 2 && action == GLFW_PRESS)
      display_top_view();
    break;
    
  case GLFW_KEY_I:
    if (Display.dimension > 2 && action == GLFW_PRESS)
      display_iso_view();
    break;
    
  case GLFW_KEY_V:
    if (Display.dimension > 2 && action == GLFW_PRESS)
      for (int i = 0; i < 4; i++)
        printf ("Camera.rotation[%d] = %g;\n", i, Camera.rotation[i]);
    break;
    
  }
}

void display_connect (const char * host, int port)
{
  struct transfer *transfer = transfer_open (host, port, reconnect_ms);
  if (!transfer) {
    fprintf (stderr, "gview: could not connect to %s:%d\n", host, port);
    exit (1);
  }
  Display.transfer = transfer;  
  assert (!transfer_set_throttle (transfer, high_bytes, low_bytes, high_messages, low_messages));
  display_set_connected (-1);
}

void display_disconnect()
{
  transfer_close (Display.transfer);
  Display.transfer = NULL;
  display_set_connected (0);
}

/**
# Loading parameters */

typedef int Colormap;

static const Colormap
  jet = 0,
  cool_warm = 1,
  gray = 2,
  randomap = 3,
  blue_white_red = 4;

typedef struct { double x, y, z; } coord;

typedef void * MapFunc, * bview;

#ifndef HUGE
# define HUGE 1e30
#endif

#include "../parse.h"

bool display_load_data (const char *filename, const char *data)
{
  char *copy = strdup (data);
  if (!copy)
    return false;

  char *line = copy;
  while (line) {
    char *next = strchr (line, '\n');
    if (next)
      *next++ = '\0';

    if (*line) {
      char *command = line;
      while (strchr (" \t", *command))
        command++;

      if (strncmp (command, "view", 4)) {
#if UI
        ui_add_item (command);
#else
        add_item (command, Display.transfer);
#endif
      }
      else {
        char *s = mystrtok (remove_blanks (command), "(");
        if (!s) {
          fprintf (stderr, "gview: error reading '%s'\n", filename);
          free (copy);
          return false;
        }

        unsigned int width = 800, height = 800;
        float bg[3] = {-1, -1, -1};
        Params params[] = {
          {"tx", pfloat, &Camera.tx},
          {"ty", pfloat, &Camera.ty},
          {"tz", pfloat, &Camera.tz},
          {"quat", pfloat, Camera.rotation, 4},
          {"width", punsigned, &width},
          {"height", punsigned, &height},
          {"bg", pfloat, bg, 3},
          {NULL}
        };

        if (!parse_params (params)) {
          free (copy);
          return false;
        }

        if (bg[0] >= 0) {
          memcpy (Camera.bg, bg, 3*sizeof(float));
          Camera.bg[3] = 1.;
        }

#if UI
        ui_set_window_size (Display.window, width, height);
#else
        glfwSetWindowSize (Display.window, width, height);
#endif
      }
    }

    line = next;
  }

  free (copy);
  return true;
}

bool display_load (const char *filename)
{
  FILE *fp = fopen (filename, "r");
  if (!fp) {
    perror (filename);
    return false;
  }

  fseek (fp, 0, SEEK_END);
  long size = ftell (fp);
  rewind (fp);

  char *data = malloc (size + 1);
  if (!data) {
    fclose (fp);
    return false;
  }

  if (fread (data, 1, size, fp) != (size_t) size) {
    free (data);
    fclose (fp);
    return false;
  }

  data[size] = '\0';
  fclose (fp);

  bool result = display_load_data (filename, data);
  free (data);
  return result;
}

/**
# GLFW/OpenGL thread */

static void main_loop (void)
{
  if (Display.window) {
    glfwPollEvents();
    if (Display.transfer)
      while (poll_messages (Display.transfer) > 0);
    display_render (Display.window);
#if UI
    ui_display (Display.window);
#endif
    glfwSwapBuffers (Display.window);
  }
}

int display_main (const char * host, int port, const char * gv)
{
  display_connect (host, port);
#if UI
  ui_set_host (host, port);
#endif

  GLFWwindow *window;

  if (!glfwInit()) {
    fprintf (stderr, "gview: glfwInit() failed\n");
    return 1;
  }

  glfwWindowHint (GLFW_SAMPLES, 4);
#ifdef __EMSCRIPTEN__
  glfwWindowHint (GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
  glfwWindowHint (GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint (GLFW_CONTEXT_VERSION_MINOR, 0);
#endif
  window = glfwCreateWindow (800, 800, "gview", NULL, NULL);
  
  if (!window) {
    fprintf (stderr, "gview: glfwCreateWindow() failed\n");
    glfwTerminate();
    return 1;
  }
  
  Display.window = window;

  glfwMakeContextCurrent (window);

  if (!gladLoadGLLoader ((GLADloadproc) glfwGetProcAddress)) {
    fprintf (stderr, "gview: failed to initialize GLAD\n");
    glfwDestroyWindow (window);
    glfwTerminate ();
    return 1;
  }
  
  glfwSetScrollCallback (window, display_scroll);
  glfwSetMouseButtonCallback (window, display_mouse_button);
  glfwSetCursorPosCallback (window, display_cursor_pos);
  glfwSetKeyCallback (window, key_callback);
 
  glfwSwapInterval (0); // not sure whether zero or one

  display_gl_init();

#if UI
  ui_init (window);
  ui_display (window);
  ui_set_window_size (window, 800, 800);
#endif

  display_set_connected (reconnect_ms > 0 ? -1 : 0);

  /**
  Load parameter file (if '-g FILE' is used). */
  
  if (gv) {
    if (!display_load (gv))
      return 1;
    Display.dimension = -1; // do not reset view
  }

  /**
  # Main loop */

#ifdef __EMSCRIPTEN__
  emscripten_set_resize_callback (EMSCRIPTEN_EVENT_TARGET_WINDOW,
                                  window,
                                  EM_FALSE,
                                  resize_callback);
  // Use Emscripten's main loop for WebAssembly
  emscripten_set_main_loop (main_loop, 0, 1);
#else
  // Native GLFW loop
  while (!glfwWindowShouldClose (window)) {
#if FPS
    struct timeval tv;
    gettimeofday (&tv, NULL);
#endif
    main_loop();
#if FPS  
    struct timeval tvend;
    gettimeofday (&tvend, NULL);
    fprintf (stderr, "fps: %.0f  \r", 1./((tvend.tv_sec - tv.tv_sec) + 
                                        (tvend.tv_usec - tv.tv_usec)/1e6));
#endif
  }
#endif

  /**
  # Closing and cleanup */
  
  foreach_item (item)
    remove_item (item, Display.transfer);
  free (Display.items);
  Display.items = NULL;

  if (Display.transfer)
    send_item ("status();", Display.transfer, '-');
  
  display_disconnect();

#if UI
  ui_destroy();
#endif

  glDeleteVertexArrays (1, &vao);
  glDeleteProgram (program);

  glfwDestroyWindow (window);
  glfwTerminate();

  return 0;
}
