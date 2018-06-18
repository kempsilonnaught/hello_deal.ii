IFLAGS += -I/usr/include

LFLAGS += -L/usr/lib/x86_64-linux-gnu

LIBS += -ldeal.ii
LIBS += -ltbb

CXX = g++
CXXFLAGS += -Wall
CXXFLAGS += -O0 -ggdb
CXXFLAGS += -std=c++11

step10exe : step10.o 
	$(CXX) -o step10exe $(CXXFLAGS) $(LFLAGS) step10.o $(LIBS)

step10.o : step10.cpp
	$(CXX) -o step10.o $(CXXFLAGS) $(IFLAGS) -c step10.cpp