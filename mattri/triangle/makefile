# makefile for Triangle adapted from original makefile
#
# Type "make" to compile Triangle

SYSTEM = $(shell uname)

SRC = ./
BIN = ./

# CC should be set to the name of your favorite C compiler.

CC = gcc

# CSWITCHES is a list of all switches passed to the C compiler.  I strongly
#   recommend using the best level of optimization.  I also strongly
#   recommend timing each level of optimization to see which is the
#   best.  For instance, when I had a DEC Alpha using DEC's optimizing
#   compiler, the -O2 switch generated a notably faster version of Triangle
#   than the -O3 switch.  Go figure.
#
# By default, Triangle and Show Me use double precision floating point
#   numbers.  If you prefer single precision, use the -DSINGLE switch.
#   Double precision uses more memory, but improves the resolution of
#   the meshes you can generate with Triangle.  It also reduces the
#   likelihood of a floating exception due to overflow.  Also, it is
#   much faster than single precision on many architectures.  I recommend
#   double precision unless you want to generate a mesh for which you do
#   not have enough memory to use double precision.
#
# If yours is not a Unix system, use the -DNO_TIMER switch to eliminate the
#   Unix-specific timer code.  Also, don't try to compile Show Me; it only
#   works with X Windows.
#
# To get the exact arithmetic to work right on an Intel processor, use the
#   -DCPU86 switch with Microsoft C, or the -DLINUX switch with gcc running
#   on Linux.  The floating-point arithmetic might not be robust otherwise.
#   Please see http://www.cs.cmu.edu/~quake/robust.pc.html for details.
#
# If you are modifying Triangle, I recommend using the -DSELF_CHECK switch
#   while you are debugging.  Defining the SELF_CHECK symbol causes
#   Triangle to include self-checking code.  Triangle will execute more
#   slowly, however, so be sure to remove this switch before compiling a
#   production version.
#
# If the size of the Triangle binary is important to you, you may wish to
#   generate a reduced version of Triangle.  The -DREDUCED switch gets rid
#   of all features that are primarily of research interest.  Specifically,
#   defining the REDUCED symbol eliminates the -i, -F, -s, and -C switches.
#   The -DCDT_ONLY switch gets rid of all meshing algorithms above and beyond
#   constrained Delaunay triangulation.  Specifically, defining the CDT_ONLY
#   symbol eliminates the -r, -q, -a, -u, -D, -S, and -s switches.  The
#   REDUCED and CDT_ONLY symbols may be particularly attractive when Triangle
#   is called by another program that does not need all of Triangle's
#   features; in this case, these switches should appear as part of
#   "TRILIBDEFS" below.
#
# On some systems, you may need to include -I/usr/local/include and/or
#   -L/usr/local/lib in the compiler options to ensure that the X include
#   files and libraries that Show Me needs are found.  If you get errors
#   like "Can't find include file X11/Xlib.h", you need the former switch.
#   Try compiling without them first; add them if that fails.
#
# An example CSWITCHES line is:
#
#   CSWITCHES = -O -DNO_TIMER -DLINUX -I/usr/X11R6/include -L/usr/X11R6/lib

ifeq ($(SYSTEM),Linux)
	CSWITCHES = -O3 -DLINUX
else
	CSWITCHES = -O3
endif

# RM should be set to the name of your favorite rm (file deletion program).

RM = /bin/rm

# The action starts here.

all: $(SRC)triangle.c
	echo $(SYSTEM)
	$(CC) $(CSWITCHES) -o $(BIN)triangle $(SRC)triangle.c -lm

distclean:
	$(RM) $(BIN)triangle $(BIN)triangle.o $(BIN)tricall $(BIN)showme
