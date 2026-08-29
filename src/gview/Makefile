CFLAGS += -I$(BASILISK) -O2

# these are not Basilisk programs
EXCLUDE = display.c transfer_ws.c wsclient.c

TOPTARGETS = all clean

SUBDIRS = gui

.PHONY: subdirs $(SUBDIRS) $(TOPTARGETS)

all: subdirs gview-nogui gview

subdirs: $(SUBDIRS)

$(TOPTARGETS): $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@ $(filter-out $(SUBDIRS),$(MAKECMDGOALS))

HEADERS = glad.h ../gl/trackball.h ../parse.h transfer.h display.h

display-nogui.o: display.c $(HEADERS)
	$(CC) -c $(CFLAGS) -o $@ $<

display.o: display.c $(HEADERS) gui/ui.h
	$(CC) -c $(CFLAGS) -DUI -o $@ $<

gview: display.o transfer_ws.o glad.o \
	gui/ui.o gui/dcimgui_impl_glfw.o gui/dcimgui_impl_opengl3.o gui/dcimgui.o \
	gui/tinyfiledialogs.o
	$(CC) $^ -o $@ -L../gl -lglutils -lglfw -limgui -lstb -lstdc++ -lm

gview-nogui: display-nogui.o transfer_ws.o glad.o
	$(CC) $^ -o $@ -L../gl -lglutils -lglfw -lm

wsclient: wsclient.o transfer_ws.o
	$(CC) $^ -o $@

wsclient.o: transfer.h

transfer_ws.o: transfer.h

include ../Makefile.defs
