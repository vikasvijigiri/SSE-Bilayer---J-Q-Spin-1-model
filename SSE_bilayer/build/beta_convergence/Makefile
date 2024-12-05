CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall  -pedantic-errors -g -msse2 -march=native -funroll-loops -ffast-math -fomit-frame-pointer -fstrict-aliasing  
#-O3 -g -Wall -msse2 -march=native -funroll-loops -ffast-math -fomit-frame-pointer -fstrict-aliasing

#VAR = ${CURDIR}
#VAR = $(wildcard ../*cpp)
#_DEPS = *.cpp
#SRCS = $(patsubst %,$(SDIR)/%,$(_DEPS))
#$(info    VAR is $(VAR))
SDIR = $(wildcard ../../../../../../src)

SRCS = $(SDIR)/*.cpp
#OBJD = /home/vikas/SSE_QQ_Bi_H/build/BetaVp40/L4/J_H_0.0/
OBJS = ${SRCS:%.cpp=${CURDIR}/*.o}
OBJD = ${SRCS:%.cpp=%.o}


MAIN = main

all: ${MAIN}
	echo   Compilation Successful.................

${MAIN}: ${OBJD}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${MAIN}

.cpp.o:
	${CXX} ${CXXFLAGS} -c ${SRCS}

clean:
	${RM} ${PROGS} ${OBJS} *.o *~


	

