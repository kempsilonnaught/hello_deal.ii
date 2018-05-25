IFLAGS += -I/usr/include

LFLAGS += -L/usr/lib/x86_64-linux-gnu

LIBS += -ldeal.ii

CXXFLAGS += -Wall
CXXFLAGS += -O0 -ggdb
CXXFLAGS += -std=c++11

step1exe : step1.o
	g++ -o step1exe $(CXXFLAGS) $(LFLAGS)  step1.o $(LIBS)  

step1.o : step1.cpp
	g++ -o step1.o $(CXXFLAGS) $(IFLAGS) -c step1.cpp  