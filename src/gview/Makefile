.PHONY: all clean

CFLAGS += -I$(BASILISK) -I/usr/include/imgui -O2

all: libgview.a

# these are not Basilisk programs
EXCLUDE = display.c

libgview.a: display.o ../grid/gpu/glad.o dcimgui_impl_glfw.o dcimgui_impl_opengl3.o dcimgui.o
	ar cr libgview.a $^

display.o: display.h ../grid/gpu/glad.h \
	gui/ui.c \
	gui/squares.c

%.o: gui/bindings/%.cpp
	g++ $(CFLAGS) -I/usr/include/imgui -I/usr/include/imgui/backends -c -o $@ $<

include ../Makefile.defs
