.PHONY: all clean

CFLAGS += -I$(BASILISK) -O2

all: libgview.a

# these are not Basilisk programs
EXCLUDE = display.c

libgview.a: display.o ../grid/gpu/glad.o
	ar cr libgview.a ../grid/gpu/glad.o display.o

display.o: display.h ../grid/gpu/glad.h

include ../Makefile.defs
