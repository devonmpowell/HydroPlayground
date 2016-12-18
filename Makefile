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
SOURCES = driver.c eos.c
COMMON = common.h eos.h 
OBJ = $(SOURCES:.c=.o)
LIBOUT = lib/hydroplay.so 

# compiler options
CC = gcc
CFLAGS = -shared -fPIC -O3 -Wall
LDFLAGS += -lm

# enable MPI?
#ifeq ($(strip $(USE_MPI)), 1)
#DEF += -DUSE_MPI
#CC = mpicc
#endif

all: $(LIBOUT)

$(LIBOUT): $(COMMON) $(OBJ) dirs 
	$(CC) $(OBJ) -o $@ $(LDFLAGS) $(CFLAGS)

.c.o: $(COMMON)
	$(CC) -c -o $@ $(INC) $(CFLAGS) $<

dirs:
	@- if ! test -e obj; then mkdir obj; fi
	@- if ! test -e lib; then mkdir lib; fi

clean:
	rm -f $(OBJ) *~ core $(EXE)
