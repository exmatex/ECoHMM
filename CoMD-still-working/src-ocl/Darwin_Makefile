val = $(shell nvidia-smi -a 2>/dev/null)

ifeq ($(val),)
#@echo AMD
//STREAMSDK = /home/stream/AMD-APP-SDK-v2.4-lnx64
//STREAMSDK = /home/stream/ati-stream-sdk-v2.3-lnx64
STREAMSDK = /home/opt/AMDAPP_v2.6
INCS = -I$(STREAMSDK)/include 
STREAMSDKLIB = $(STREAMSDK)/lib/x86_64
LIBS = -L$(STREAMSDKLIB) -lOpenCL -Wl,-rpath,$(STREAMSDKLIB) -lm
else
#@echo nVidia
STREAMSDK =  /home/cuda/sdk4.1.28/OpenCL/common/
INCS = -I$(STREAMSDK)/inc 
LIBS = -L/usr/lib64/ -lOpenCL -lm
endif

ifeq ($HOSTNAME , 'cn1' )
STREAMSDK =  /usr
INCS = -I$(STREAMSDK)/include
LIBS = -L/usr/lib64/ -lOpenCL -lm
endif

HEAD = -I. -I../src-flat/include
SRCS = pmdOCL.c \
	   helpers.c  \
	   cl_utils.c  \
	   docompute.c\
	   ../src-flat/cmdLineParser.c \
	   ../src-flat/domains.c \
	   ../src-flat/eam.c \
	   ../src-flat/ljForce.c \
	   ../src-flat/ic_fcc.c \
	   ../src-flat/mycommand.c \
	   ../src-flat/read.c \
	   ../src-flat/write.c \
	   ../src-flat/utility.c \
	   ../src-flat/timestep.c

OBJS = $(SRCS:.c=.o)

DEPS = *.h
CC = gcc

all: cruft

%.o: %.c $(DEPS)
$(CC) -c -o $@ $< $(INCS) $(HEAD)

cruft: $(OBJS)
$(CC) -o ../cruftOCL -O3 $(OBJS) $(HEAD)  $(INCS) $(LIBS)

clean:
rm *.o ../src-flat/*.o

