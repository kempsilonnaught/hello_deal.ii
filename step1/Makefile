IFLAGS += -I/usr/include

LFLAGS += -L/usr/lib/x86_64-linux-gnu

LIBS += -ldeal.ii

CXXFLAGS += -Wall
CXXFLAGS += -O0 -ggdb

step1exe : step1.o
	g++ $(CXXFLAGS) $(LFLAGS) $(LIBS) step1.o -o step1exe 

step1 : step1.cpp
	g++ $(CXXFLAGS) $(IFLAGS) -o step1.cpp -c step1