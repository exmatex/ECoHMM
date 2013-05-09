# This is a skeleton makefile in case cmake is having problems..
# All the paths need to be adjusted to your system. 
# This file provides the base SRCS that you need.

# Uncomment the appropriate lines to enable build options
# Precision option

#PREC_FLAG = -DDOUBLE
PREC_FLAG = -DSINGLE

#CHEBY_FLAG = -DUSE_CHEBY 

OPT_FLAG = -O3
FLAGS = $(PREC_FLAG) $(CHEBY_FLAG) ${OPT_FLAG} -DGZIPSUPPORT



#INCS = {path to opencl.h}
#LIBS = {path to libopencl.so}

LIBS = -framework OpenCL -lz

HEAD = 	 -I. -I../src-flat/include -I./include -DUSEOCL

SRCS =  pmdOCL.c \
	helpers.c \
	cl_utils.c \
	docompute.c\
	../src-flat/cmdLineParser.c \
        ../src-flat/domains.c \
        ../src-flat/eam.c \
        ../src-flat/ljForce.c \
        ../src-flat/ic_fcc.c \
        ../src-flat/cheby.c \
        ../src-flat/mycommand.c \
        ../src-flat/read.c \
        ../src-flat/write.c \
        ../src-flat/utility.c \
        ../src-flat/timestep.c

OBJS = $(SRCS:.c=.o)

DEPS = ./include/*.h ../src-flat/include/*.h
CC = gcc

all: CoMDOCL

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(INCS) $(HEAD) $(FLAGS)

CoMDOCL : $(OBJS)

	$(CC) -o ../CoMDOCL $(OBJS) $(HEAD)  $(INCS) $(LIBS)

clean:
	rm *.o ../src-flat/*.o

