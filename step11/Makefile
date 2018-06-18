IFLAGS += -I/usr/include

LFLAGS += -L/usr/lib/x86_64-linux-gnu

LIBS += -ldeal.ii
LIBS += -ltbb

CXX = g++
CXXFLAGS += -Wall 
CXXFLAGS += -O0 -ggdb
CXXFLAGS += -std=c++11

step11exe : step11.o
	$(CXX) -o step11exe $(CXXFLAGS) $(LFLAGS) step11.o $(LIBS)

step11.o : step11.cpp
	$(CXX) -o step11.o $(CXXFLAGS) $(IFLAGS) -c step11.cpp