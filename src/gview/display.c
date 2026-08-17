/**
# Gview

To install

~~~bash
sudo apt install libimgui-dev
~~~
*/

#include "display.h"

#include "../grid/gpu/glad.h"
#include <GLFW/glfw3.h>

#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

GviewParameters Gview = {0};

typedef struct {
  pthread_t thread;
  pthread_mutex_t mutex;

  DisplayFrame *ready;

  bool stop;
  bool joined;
} DisplayState;

static DisplayState Display;

const float mouse_zoom_sensitivity = 10.f; // smaller is more sensitive

/**
OpenGL state belongs exclusively to the display thread. */

static GLuint vao;
static GLuint position_vbo;
static GLuint normal_vbo;
static GLuint color_vbo;
static GLuint index_vbo, lindex_vbo;

static GLint constant_color;
static GLint use_constant_color;
static bool constant_color_enabled;
static float constant_color_value[3];

static GLint constant_normal;
static GLint use_constant_normal;
static bool constant_normal_enabled;
static float constant_normal_value[3];

static GLint lighting;
static GLint light_direction;

static GLuint program;

static size_t index_count, lindex_count;

/**
# Camera

This is deliberately independent of Basilisk. We will replace/extend
this with the jsview-style interaction later. */

typedef struct {
  float rotation[4];
  float modelview[16];
  float projection[16];
  float zoom;
  float pan_x, pan_y;
  double mouse_x, mouse_y;
  bool rotating;
  bool panning;
  bool zooming;
} DisplayCamera;

static DisplayCamera Camera;

static void matrix_identity (float *m)
{
  memset (m, 0, 16*sizeof(float));

  m[0] = 1.f;
  m[5] = 1.f;
  m[10] = 1.f;
  m[15] = 1.f;
}

static void matrix_ortho (float *m,
                          float left, float right,
                          float bottom, float top,
                          float near, float far)
{
  matrix_identity (m);

  m[0] = 2.f/(right - left);
  m[5] = 2.f/(top - bottom);
  m[10] = -2.f/(far - near);

  m[12] = -(right + left)/(right - left);
  m[13] = -(top + bottom)/(top - bottom);
  m[14] = -(far + near)/(far - near);
}

static void display_arcball (GLFWwindow *window, double x, double y, float p[3])
{
  int width, height;
  float aspect, nx, ny, r;

  glfwGetWindowSize (window, &width, &height);

  aspect = (float) width/height;
  nx = (2.f*(float) x/width - 1.f)*aspect;
  ny = 1.f - 2.f*(float) y/height;

  r = nx*nx + ny*ny;

  if (r <= 1.f)
    p[2] = sqrtf (1.f - r);
  else {
    float scale = 1.f/sqrtf (r);
    nx *= scale;
    ny *= scale;
    p[2] = 0.f;
  }

  p[0] = nx;
  p[1] = ny;

  float length = sqrtf (p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

  p[0] /= length;
  p[1] /= length;
  p[2] /= length;
}

static void display_arcball_rotate (float q[4], const float a[3], const float b[3])
{
  float cross[3] = {
    a[1]*b[2] - a[2]*b[1],
    a[2]*b[0] - a[0]*b[2],
    a[0]*b[1] - a[1]*b[0]
  };
  float dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

  q[0] = 1.f + dot;
  q[1] = cross[0];
  q[2] = cross[1];
  q[3] = cross[2];

  float length = sqrtf (q[0]*q[0] + q[1]*q[1] +
                        q[2]*q[2] + q[3]*q[3]);

  if (length > 0.f) {
    q[0] /= length;
    q[1] /= length;
    q[2] /= length;
    q[3] /= length;
  }
}

static void display_quaternion_multiply (float r[4], const float a[4], const float b[4])
{
  float q[4];

  q[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
  q[1] = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2];
  q[2] = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1];
  q[3] = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0];

  memcpy (r, q, sizeof(q));
}

static void display_quaternion_matrix (float m[16], const float q[4])
{
  float w = q[0], x = q[1], y = q[2], z = q[3];

  matrix_identity (m);

  m[0] = 1.f - 2.f*(y*y + z*z);
  m[1] = 2.f*(x*y + z*w);
  m[2] = 2.f*(x*z - y*w);

  m[4] = 2.f*(x*y - z*w);
  m[5] = 1.f - 2.f*(x*x + z*z);
  m[6] = 2.f*(y*z + x*w);

  m[8] = 2.f*(x*z + y*w);
  m[9] = 2.f*(y*z - x*w);
  m[10] = 1.f - 2.f*(x*x + y*y);
}

static void display_zoom (GLFWwindow *window, double factor)
{
  int width, height;
  double x, y;
  float aspect, old_zoom, new_zoom;
  float nx, ny;

  glfwGetCursorPos (window, &x, &y);
  glfwGetWindowSize (window, &width, &height);

  if (!width || !height)
    return;

  old_zoom = Camera.zoom;
  new_zoom = old_zoom*factor;

  if (new_zoom < 0.001f)
    new_zoom = 0.001f;
  else if (new_zoom > 100.f)
    new_zoom = 100.f;

  if (new_zoom == old_zoom)
    return;

  aspect = (float) width/height;
  nx = 2.f*(float) x/width - 1.f;
  ny = 1.f - 2.f*(float) y/height;

  Camera.pan_x += nx*aspect*(new_zoom - old_zoom);
  Camera.pan_y += ny*(new_zoom - old_zoom);
  Camera.zoom = new_zoom;
}

static void display_scroll (GLFWwindow *window, double xoffset, double yoffset)
{
  (void) xoffset;
  display_zoom (window, powf (1.1f, (float) -yoffset));
}

static void display_mouse_button (GLFWwindow *window, int button, int action, int mods)
{
  (void) mods;

  if (Gview.rotation &&
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
  int width, height;
  double dx, dy;
  
  if (Camera.rotating) {
    float a[3], b[3], q[4], rotation[4];

    display_arcball (window, Camera.mouse_x, Camera.mouse_y, a);
    display_arcball (window, x, y, b);
    display_arcball_rotate (q, a, b);

    display_quaternion_multiply (rotation, q, Camera.rotation);
    memcpy (Camera.rotation, rotation, sizeof(rotation));

    Camera.mouse_x = x;
    Camera.mouse_y = y;
  }
  else if (Camera.panning) {
    glfwGetWindowSize (window, &width, &height);

    if (width && height) {
      dx = x - Camera.mouse_x;
      dy = y - Camera.mouse_y;

      Camera.pan_x += 2.f*Camera.zoom*(float) dx/width;
      Camera.pan_y -= 2.f*Camera.zoom*(float) dy/height;

      Camera.mouse_x = x;
      Camera.mouse_y = y;
    }
  }
  else if (Camera.zooming) {
    dy = y - Camera.mouse_y;

    if (dy) {
      display_zoom (window, powf (1.1f, (float) dy/mouse_zoom_sensitivity));
    }

    Camera.mouse_y = y;
  }
}

/**
# DisplayFrame ownership. */

static void display_frame_free (DisplayFrame *frame)
{
  if (!frame)
    return;
  free (frame->position);
  free (frame->normal);
  free (frame->color);
  free (frame->index);
  free (frame->lindex);
  free (frame);
}

/**
The simulation thread publishes a pointer.

No geometry is copied here. The display module takes
ownership of the complete DisplayFrame. */

void display_publish (DisplayFrame *frame)
{
  DisplayFrame *old;

  pthread_mutex_lock (&Display.mutex);

  old = Display.ready;
  Display.ready = frame;

  pthread_mutex_unlock (&Display.mutex);

  /**
  Only the newest frame is retained. */
  
  display_frame_free (old);
}

/**
The display thread takes ownership of the newest frame. */

static DisplayFrame *display_take_frame (void)
{
  DisplayFrame *frame;

  pthread_mutex_lock (&Display.mutex);

  frame = Display.ready;
  Display.ready = NULL;

  pthread_mutex_unlock (&Display.mutex);

  return frame;
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
  static const char *vertex_source =
    "#version 430\n"
    "\n"
    "layout (location = 0) in vec3 position;\n"
    "layout (location = 1) in vec3 normal;\n"
    "layout (location = 2) in vec3 color;\n"
    "\n"
    "uniform mat4 modelview;\n"
    "uniform mat4 projection;\n"
    "uniform vec3 constant_normal;\n"
    "uniform bool use_constant_normal;\n"
    "uniform vec3 constant_color;\n"
    "uniform bool use_constant_color;\n"
    "uniform bool lighting;\n"
    "\n"
    "out vec3 vertex_normal;\n"
    "out vec3 vertex_color;\n"
    "\n"
    "void main ()\n"
    "{\n"
    "  gl_Position = projection * modelview * vec4 (position, 1.);\n"
    "  if (lighting)\n"
    "    vertex_normal = mat3 (modelview)*(use_constant_normal ? constant_normal : normal);\n"
    "  vertex_color = use_constant_color ? constant_color : color;\n"
    "}\n";

  static const char *fragment_source =
    "#version 430\n"
    "\n"
    "in vec3 vertex_normal;\n"
    "in vec3 vertex_color;\n"
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
    "    intensity = 0.25 + 0.75*diffuse;\n"
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

static void display_gl_init (void)
{
  program = create_program ();

  if (!program)
    exit (1);

  Camera.rotation[0] = 1.f;
  Camera.rotation[1] = 0.f;
  Camera.rotation[2] = 0.f;
  Camera.rotation[3] = 0.f;
  Camera.zoom = 1.f;
  Camera.pan_x = Camera.pan_y = 0.f;
  Camera.rotating = Camera.panning = Camera.zooming = false;

  matrix_identity (Camera.modelview);
  matrix_identity (Camera.projection);
  
  constant_color = glGetUniformLocation (program, "constant_color");
  use_constant_color = glGetUniformLocation (program, "use_constant_color");
  
  constant_normal = glGetUniformLocation (program, "constant_normal");
  use_constant_normal = glGetUniformLocation (program, "use_constant_normal");

  lighting = glGetUniformLocation (program, "lighting");
  light_direction = glGetUniformLocation (program, "light_direction");
  
  glGenVertexArrays (1, &vao);

  glGenBuffers (1, &position_vbo);
  glGenBuffers (1, &normal_vbo);
  glGenBuffers (1, &color_vbo);
  glGenBuffers (1, &index_vbo);
  glGenBuffers (1, &lindex_vbo);

  glBindVertexArray (vao);

  glBindBuffer (GL_ARRAY_BUFFER, position_vbo);
  glEnableVertexAttribArray (0);
  glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), NULL);

  glBindBuffer (GL_ARRAY_BUFFER, normal_vbo);
  glEnableVertexAttribArray (1);
  glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), NULL);

  glBindBuffer (GL_ARRAY_BUFFER, color_vbo);
  glEnableVertexAttribArray (2);
  glVertexAttribPointer (2, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), NULL);

  glBindVertexArray (0);

  glEnable (GL_DEPTH_TEST);

  matrix_identity (Camera.projection);
  Camera.zoom = 1.f;
  Camera.pan_x = Camera.pan_y = 0.f;
  Camera.panning = false;
  Camera.zooming = false;
}

/**
Upload a newly received frame */

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


/**
Render the current VBO */

static void display_render (GLFWwindow *window)
{
  int width, height;
  float aspect, sx, sy;

  glfwGetFramebufferSize (window, &width, &height);

  if (!height)
    return;

  glViewport (0, 0, width, height);

  aspect = (float) width/height;
  sx = Camera.zoom*aspect;
  sy = Camera.zoom;

  matrix_ortho (Camera.projection,
                -sx, sx, -sy, sy, -100.f, 100.f);

  display_quaternion_matrix (Camera.modelview, Camera.rotation);
  Camera.modelview[12] = Camera.pan_x;
  Camera.modelview[13] = Camera.pan_y;

  glClearColor (0.15f, 0.15f, 0.15f, 1.f);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glUseProgram (program);

  glUniformMatrix4fv (glGetUniformLocation (program, "modelview"),
                      1, GL_FALSE, Camera.modelview);
  glUniformMatrix4fv (glGetUniformLocation (program, "projection"),
                      1, GL_FALSE, Camera.projection);
  
  glUniform1i (lighting, Gview.lighting);
  glUniform3f (light_direction, 0.3f, 0.5f, 1.f);
 
  glUniform1i (use_constant_normal, constant_normal_enabled);

  if (constant_normal_enabled)
    glUniform3fv (constant_normal, 1, constant_normal_value);
  
  glUniform1i (use_constant_color, constant_color_enabled);
  
  if (constant_color_enabled)
    glUniform3fv (constant_color, 1, constant_color_value);
  
  glBindVertexArray (vao);

  if (index_count) {
    glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, index_vbo);
    glDrawElements (GL_TRIANGLES, (GLsizei) index_count,
                    GL_UNSIGNED_INT, NULL);
  }
  if (lindex_count) {
    glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, lindex_vbo);
    glDrawElements (GL_LINES, (GLsizei) lindex_count,
                    GL_UNSIGNED_INT, NULL);
  }

  glBindVertexArray (0);
}

/**
# GLFW/OpenGL thread */

#include "gui/ui.c"

static void *display_thread (void *arg)
{
  GLFWwindow *window;

  (void) arg;

  if (!glfwInit()) {
    fprintf (stderr, "display: glfwInit() failed\n");
    return NULL;
  }

  window = glfwCreateWindow (window_width, window_height,
                             "gview",
                             NULL, NULL);

  if (!window) {
    fprintf (stderr, "display: glfwCreateWindow() failed\n");
    glfwTerminate();
    return NULL;
  }

  glfwMakeContextCurrent (window);

  if (!gladLoadGLLoader ((GLADloadproc) glfwGetProcAddress)) {
    fprintf (stderr, "display: failed to initialize GLAD\n");
    glfwDestroyWindow (window);
    glfwTerminate ();
    return NULL;
  }
  
  glfwSetScrollCallback (window, display_scroll);
  glfwSetMouseButtonCallback (window, display_mouse_button);
  glfwSetCursorPosCallback (window, display_cursor_pos);
 
  /**
  The display thread runs independently of the simulation.*/
  
  glfwSwapInterval (1);

  display_gl_init();
  ui_init (window);
  
  while (!glfwWindowShouldClose (window)) {
    DisplayFrame *frame;

    glfwPollEvents();

    pthread_mutex_lock (&Display.mutex);
    bool stop = Display.stop;
    pthread_mutex_unlock (&Display.mutex);

    if (stop)
      break;

    frame = display_take_frame();

    if (frame) {
      display_upload (frame);
      display_frame_free (frame);
    }

    display_render (window);
    ui (window);
    
    glfwSwapBuffers (window);
  }

  ui_destroy();
  
  glDeleteBuffers (1, &position_vbo);
  glDeleteBuffers (1, &normal_vbo);
  glDeleteBuffers (1, &color_vbo);
  glDeleteBuffers (1, &index_vbo);
  glDeleteBuffers (1, &lindex_vbo);
  glDeleteVertexArrays (1, &vao);
  glDeleteProgram (program);

  glfwDestroyWindow (window);
  glfwTerminate ();

  return NULL;
}


int display_start (void)
{
  memset (&Display, 0, sizeof(Display));
  if (pthread_mutex_init (&Display.mutex, NULL))
    return -1;
  if (pthread_create (&Display.thread, NULL,
                      display_thread, NULL)) {
    pthread_mutex_destroy (&Display.mutex);
    return -1;
  }
  return 0;
}

void display_wait (void)
{
  if (!Display.joined) {
    pthread_join (Display.thread, NULL);
    Display.joined = true;
  }
}

void display_stop (void)
{
  pthread_mutex_lock (&Display.mutex);

  Display.stop = true;

  pthread_mutex_unlock (&Display.mutex);

  display_wait ();

  display_frame_free (Display.ready);
  Display.ready = NULL;

  pthread_mutex_destroy (&Display.mutex);
}
