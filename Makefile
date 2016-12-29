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
SOURCES = src/driver.c src/eos.c src/radiation.c
COMMON = include/common.h include/eos.h Makefile hydroPlayground.py 
LIBOUT = hydroplay.so 

# compiler options
CC = gcc
CFLAGS = -shared -fPIC -O3 -Wall
INC = -I./include 
OBJ = $(SOURCES:.c=.o)
LDFLAGS += -lm


# Makefile rules!

all: $(LIBOUT)

$(LIBOUT): $(COMMON) $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS) $(CFLAGS)

.c.o: $(COMMON)
	$(CC) -c -o $@ $(INC) $(CFLAGS) $<

clean:
	rm -rf src/*.o $(LIBOUT) *~ core
