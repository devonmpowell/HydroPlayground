# ---------------------------------------------------------------------------
#
#	Makefile for d_euler_hydro
#
#	Devon Powell 
#	(with some useful bits from Jonathan Zrake)
#
#	Do not modify this file!
#	All user-set options live in Makefile.in
#
#	usage: make
#
# ---------------------------------------------------------------------------

# if there is no Makefile.in then use the template
ifneq ($(strip $(MAKEFILE_IN)),)
# use value of MAKEFILE_IN if provided on the command line
else ifeq ($(shell test -e Makefile.in && echo 1), 1)
MAKEFILE_IN = Makefile.in
else
MAKEFILE_IN = Makefile.in.template
endif
include $(MAKEFILE_IN)

# Source files
SOURCES = driver.c eos.c
COMMON = common.h eos.h 
OBJ = $(SOURCES:.c=.o)
EXE = d_euler_hydro

LDFLAGS += -lm

# Set up HDF5 dependencies
#INC += -I$(HDF5_HOME)/include
#LIB += -L$(HDF5_HOME)/lib
#LDFLAGS += -lhdf5

# enable MPI?
#ifeq ($(strip $(USE_MPI)), 1)
#DEF += -DUSE_MPI
#CC = mpicc
#endif

all: $(EXE)

$(EXE): $(COMMON) $(OBJ) 
	$(CC) $(LIB) $(OBJ) -o $@ $(LDFLAGS)

.c.o: $(COMMON)
	$(CC) $(INC) $(DEF) $(CFLAGS) $(OPT) -c $< -o $@

clean:
	rm -f $(OBJ) *~ core $(EXE)
