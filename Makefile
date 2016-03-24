CC=mpixlc
CXX=mpixlcxx

CFLAGS=-O3 -g #-pg

SRCS = 	vacf.cpp \
		msd.cpp \
		histo.cpp \
		fft.cpp \
		modalysis.cpp \
		driver.cpp
		
OBJS = 	$(SRCS:.cpp=.o)

TARGET = modalysis

INC	 += -I/soft/libraries/alcf/current/xl/FFTW3/include
LIBS += -L/soft/perftools/hpctw/lib -lmpihpm -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm -lrt -lstdc++ 
LIBS += -L/soft/libraries/alcf/current/xl/FFTW3/lib -lfftw3

all:    $(TARGET)
		@echo Compilation done.

%.o:%.cpp
		$(CXX) $(CFLAGS) $(INC) -c $< -o $@

$(TARGET): $(OBJS) 
		$(CXX) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

clean:
		$(RM) *.o *~ $(TARGET)

