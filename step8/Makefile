IFLAGS += -I/usr/include

LFLAGS += -L/usr/lib/x86_64-linux-gnu

LIBS += -ldeal.ii
LIBS += -ltbb

CXX = g++
CXXFLAGS += -Wall
CXXFLAGS += -O0 -ggdb
CXXFLAGS += -std=c++11

step8exe : step8.o
	$(CXX) -o step8exe $(CXXFLAGS) $(LFLAGS) step8.o $(LIBS)

step8.o : step8.cpp
	$(CXX) -o step8.o $(CXXFLAGS) $(IFLAGS) -c step8.cpp