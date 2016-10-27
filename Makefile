ifdef NERSC_HOST
CC=cc
CXX=CC
CFLAGS=-g -O3 -DNERSC_HOST #-pg
else
CC=mpixlc
CXX=mpixlcxx
CFLAGS=-g -O3 #-pg
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

ifdef NERSC_HOST 
LIBS += -lfftw3
else
INC	 += -I/soft/libraries/alcf/current/xl/FFTW3/include
LIBS += -L/soft/perftools/hpctw/lib -lmpihpm -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm -lrt -lstdc++ 
LIBS += -L/soft/libraries/alcf/current/xl/FFTW3/lib -lfftw3_mpi -lfftw3 -lm #-lfftw3f
endif

all:    $(TARGET)
		@echo Compilation done.

%.o:%.cpp
		$(CXX) $(CFLAGS) $(INC) -c $< -o $@

$(TARGET): $(OBJS) 
		$(CXX) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

clean:
		$(RM) *.o *~ $(TARGET)

