# This is a skeleton makefile in case cmake is having problems..
# All the paths need to be adjusted to your system. 
# This file provides the base SRCS that you need.

# Uncomment the appropriate lines to enable build options
#
# Precision option

#PREC_FLAG = -DDOUBLE
PREC_FLAG = -DSINGLE

#Chebyshev coefficients for EAM

#CHEBY_FLAG = -DUSE_CHEBY

OPT_FLAG = -O3

FLAGS = $(PREC_FLAG) $(CHEBY_FLAG) $(OPT_FLAG) -DGZIPSUPPORT

# gzip library
#
LIBS = -lz

HEAD = 	 -I. -I./include 

SRCS =  pmd.c \
	docompute.c\
	cmdLineParser.c \
        domains.c \
        eam.c \
        ljForce.c \
        ic_fcc.c \
        cheby.c \
        mycommand.c \
        read.c \
        write.c \
        utility.c \
        timestep.c \
	yamlOutput.c

OBJS = $(SRCS:.c=.o)

DEPS = ./include/*.h
CC = gcc

all: generate_info CoMD

%.o: %.c $(DEPS)
	$(CC)  -c -o $@ $< $(INCS) $(HEAD) $(FLAGS)

CoMD : $(OBJS)

	$(CC) -o ../CoMD $(OBJS) $(HEAD) $(INCS) $(LIBS)

generate_info : 
	./generate_info_header "$(CC)" "$(FLAGS)"

clean:
	rm *.o include/CoMD_info.h


