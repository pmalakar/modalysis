ifdef NERSC_HOST
CC=cc
CXX=CC
CFLAGS=-g -O3 -DNERSC_HOST #-pg
FFTW_DIR=/opt/cray/pe/fftw/3.3.8.2/x86_64
else 
ifdef ALCF_HOST
CC=mpixlc
CXX=mpixlcxx
CFLAGS=-g -O3 -DALCF_HOST #-pg
else
CC=mpicc
CXX=mpicxx
CFLAGS=-g -O3 #-pg
FFTW_DIR=         #wherever FFTW is installed
endif
endif

#CFLAGS+=-DDEBUG

SRCS = 	vacf.cpp \
		msd.cpp \
		histo.cpp \
		fft1d.cpp \
		modalysis.cpp \
		postprocess.cpp \
		coanalysis.cpp \
		driver.cpp
		
OBJS = 	$(SRCS:.cpp=.o)

TARGET = modalysis

ifdef ALCF_HOST
INC	 += -I/soft/libraries/alcf/current/xl/FFTW3/include
LIBS += -L/soft/perftools/hpctw/lib -lmpihpm -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm -lrt -lstdc++ 
LIBS += -L/soft/libraries/alcf/current/xl/FFTW3/lib -lfftw3_mpi -lfftw3 -lm #-lfftw3f
else
INC  += -I$(FFTW_DIR)/include
LIBS += -L$(FFTW_DIR)/lib -lfftw3
endif

all:    $(TARGET)
		@echo Compilation done.

%.o:%.cpp
		$(CXX) $(CFLAGS) $(INC) -c $< -o $@

$(TARGET): $(OBJS) 
		$(CXX) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

clean:
		$(RM) *.o *~ $(TARGET)

