# ---------------------------------------------------------------------------
#
#	Makefile for HydroPlayground 
#
#	Devon Powell 
#	December 2016
#
#	usage: make
#
# ---------------------------------------------------------------------------


# Source files
SOURCES = src/driver.c src/eos.c src/radiation.c src/geometry.c
COMMON = include/common.h include/eos.h include/geometry.h Makefile hydroPlayground.py 
LIBOUT = hydroplay.so 

# compiler options
CC = gcc
CFLAGS = -shared -fPIC -O3 #-Wall
INC = -I./include -I../r3d 
OBJ = $(SOURCES:.c=.o)
LDFLAGS += -L../r3d -lm -lr3d


# Makefile rules!

all: $(LIBOUT)

$(LIBOUT): $(COMMON) $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS) $(CFLAGS)

.c.o: $(COMMON)
	$(CC) -c -o $@ $(INC) $(CFLAGS) $<

clean:
	rm -rf src/*.o $(LIBOUT) *~ core
