CC=gcc
CXX=g++
RM=rm -f

CCFLAGS=-g -O2 -Wall
CXXFLAGS=-g -O2 -I../ --std=c++11 -Wall
LDFLAGS=
LDLIBS=-lm



CXX_HDRS=timer.h mtw_discrete.h
C_HDRS=timer.h

CC_SRCS=test_rng.cc
CC_OBJS=$(CC_SRCS:.cc=.o)

C_SRCS=kfunc.c
C_OBJS=$(C_SRCS:.c=.o)

.cc.o: $(CXX_HDRS)
	$(CXX) -c $(CXXFLAGS) $*.cc

.c.o: $(C_HDRS)
	$(CXX) -c $(CCFLAGS) $*.c


all: test_rng

test_rng: $(CC_OBJS) $(C_OBJS)
	$(CXX) $(LDFLAGS) -o test_rng $(CC_OBJS) $(C_OBJS) $(LDLIBS) 

test_rng.o: test_rng.cc timer.h


clean:
	$(RM) $(CC_OBJS) $(C_OBJS)

dist-clean: clean
	$(RM) test_rng
