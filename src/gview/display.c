#include "display.h"

#include "../grid/gpu/glad.h"
#include <GLFW/glfw3.h>

#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  pthread_t thread;
  pthread_mutex_t mutex;

  DisplayFrame *ready;

  bool stop;
} DisplayState;

static DisplayState Display;

/**
OpenGL state belongs exclusively to the display thread. */

static GLuint vao;
static GLuint position_vbo;
static GLuint normal_vbo;
static GLuint color_vbo;
static GLuint index_vbo;
static GLuint program;

static size_t index_count;
static int primitive;

/**
# Camera

This is deliberately independent of Basilisk. We will replace/extend
this with the jsview-style interaction later. */

typedef struct {
  float modelview[16];
  float projection[16];
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
    "\n"
    "out vec3 vertex_color;\n"
    "\n"
    "void main ()\n"
    "{\n"
    "  gl_Position = projection * modelview * vec4 (position, 1.);\n"
    "  vertex_color = color;\n"
    "}\n";

  static const char *fragment_source =
    "#version 430\n"
    "\n"
    "in vec3 vertex_color;\n"
    "out vec4 frag_color;\n"
    "\n"
    "void main ()\n"
    "{\n"
    "  frag_color = vec4 (vertex_color, 1.);\n"
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

    fprintf (stderr,
             "display: shader linking failed:\n%s\n",
             log);

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

  glGenVertexArrays (1, &vao);

  glGenBuffers (1, &position_vbo);
  glGenBuffers (1, &normal_vbo);
  glGenBuffers (1, &color_vbo);
  glGenBuffers (1, &index_vbo);

  glBindVertexArray (vao);

  glBindBuffer (GL_ARRAY_BUFFER, position_vbo);
  glEnableVertexAttribArray (0);
  glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE,
                         3*sizeof(float), NULL);

  glBindBuffer (GL_ARRAY_BUFFER, normal_vbo);
  glEnableVertexAttribArray (1);
  glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE,
                         3*sizeof(float), NULL);

  glBindBuffer (GL_ARRAY_BUFFER, color_vbo);
  glEnableVertexAttribArray (2);
  glVertexAttribPointer (2, 3, GL_FLOAT, GL_FALSE,
                         3*sizeof(float), NULL);

  glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, index_vbo);

  glBindVertexArray (0);

  glEnable (GL_DEPTH_TEST);

  matrix_identity (Camera.modelview);
  matrix_identity (Camera.projection);
}


/**
Upload a newly received frame */

static void display_upload (DisplayFrame *frame)
{
  glBindVertexArray (vao);

  glBindBuffer (GL_ARRAY_BUFFER, position_vbo);
  glBufferData (GL_ARRAY_BUFFER,
                3*frame->nvertex*sizeof(float),
                frame->position,
                GL_DYNAMIC_DRAW);

  glBindBuffer (GL_ARRAY_BUFFER, normal_vbo);
  glBufferData (GL_ARRAY_BUFFER,
                3*frame->nvertex*sizeof(float),
                frame->normal,
                GL_DYNAMIC_DRAW);

  glBindBuffer (GL_ARRAY_BUFFER, color_vbo);
  glBufferData (GL_ARRAY_BUFFER,
                3*frame->nvertex*sizeof(float),
                frame->color,
                GL_DYNAMIC_DRAW);

  glBindBuffer (GL_ELEMENT_ARRAY_BUFFER, index_vbo);
  glBufferData (GL_ELEMENT_ARRAY_BUFFER,
                frame->nindex*sizeof(unsigned int),
                frame->index,
                GL_DYNAMIC_DRAW);

  index_count = frame->nindex;
  primitive = frame->primitive;

  glBindVertexArray (0);
}


/**
Render the current VBO */

static void display_render (GLFWwindow *window)
{
  int width, height;

  glfwGetFramebufferSize (window, &width, &height);

  if (!height)
    return;

  glViewport (0, 0, width, height);

  /**
  Temporary projection. The proper camera will replace
  this once we add rotation/panning/zooming. */
  
  matrix_ortho (Camera.projection,
                -1.f, 1.f,
                -1.f, 1.f,
                -100.f, 100.f);

  glClearColor (0.15f, 0.15f, 0.15f, 1.f);

  glClear (GL_COLOR_BUFFER_BIT |
           GL_DEPTH_BUFFER_BIT);

  glUseProgram (program);

  glUniformMatrix4fv (
    glGetUniformLocation (program, "modelview"),
    1, GL_FALSE, Camera.modelview);

  glUniformMatrix4fv (
    glGetUniformLocation (program, "projection"),
    1, GL_FALSE, Camera.projection);

  glBindVertexArray (vao);

  if (index_count) {
    GLenum mode = primitive == 0 ?
      GL_LINES : GL_TRIANGLES;

    glDrawElements (mode,
                    (GLsizei) index_count,
                    GL_UNSIGNED_INT,
                    NULL);
  }

  glBindVertexArray (0);
}


/**
# GLFW/OpenGL thread */

static void *display_thread (void *arg)
{
  GLFWwindow *window;

  (void) arg;

  if (!glfwInit()) {
    fprintf (stderr, "display: glfwInit() failed\n");
    return NULL;
  }

  window = glfwCreateWindow (800, 800,
                             "Basilisk",
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
  
  /**
  The display thread runs independently of the simulation.*/
  
  glfwSwapInterval (1);

  display_gl_init ();

  while (!glfwWindowShouldClose (window)) {
    DisplayFrame *frame;

    glfwPollEvents ();

    pthread_mutex_lock (&Display.mutex);
    bool stop = Display.stop;
    pthread_mutex_unlock (&Display.mutex);

    if (stop)
      break;

    frame = display_take_frame ();

    if (frame) {
      display_upload (frame);
      display_frame_free (frame);
    }

    display_render (window);

    glfwSwapBuffers (window);
  }

  glDeleteBuffers (1, &position_vbo);
  glDeleteBuffers (1, &normal_vbo);
  glDeleteBuffers (1, &color_vbo);
  glDeleteBuffers (1, &index_vbo);
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


void display_stop (void)
{
  pthread_mutex_lock (&Display.mutex);

  Display.stop = true;

  pthread_mutex_unlock (&Display.mutex);

  pthread_join (Display.thread, NULL);

  display_frame_free (Display.ready);
  Display.ready = NULL;

  pthread_mutex_destroy (&Display.mutex);
}
