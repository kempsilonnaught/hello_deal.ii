IFLAGS += -I/usr/include

LFLAGS += -L/usr/lib/x86_64-linux-gnu

LIBS += -ldeal.ii

CXXFLAGS += -Wall
CXXFLAGS += -O0 -ggdb
CXXFLAGS += -std=c++11

step2exe : step2.o
	g++ -o step2exe $(CXXFLAGS) $(LFLAGS) step2.o $(LIBS)

step2.o : step2.cpp
	g++ -o step2.o $(CXXFLAGS) $(IFLAGS) -c step2.cpp