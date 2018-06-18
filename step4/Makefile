IFLAGS += -I/usr/include

LFLAGS += -L/usr/lib/x86_64-linux-gnu

LIBS += -ldeal.ii
LIBS += -ltbb

CXX = g++
CXXFLAGS += -Wall
CXXFLAGS += -O0 -ggdb
CXXFLAGS += -std=c++11

step4exe : step4.o
	$(CXX)	-o step4exe $(CXXFLAGS) $(LFLAGS) step4.o $(LIBS)

step4.o : step4.cpp
	$(CXX)	-o step4.o $(CXXFLAGS) $(IFLAGS) -c step4.cpp
