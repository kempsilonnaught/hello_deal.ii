IFLAGS += -I/usr/include

LFLAGS += -L/usr/lib/x86_64-linux-gnu

LIBS += -ldeal.ii
LIBS += -ltbb

CXX = g++
CXXFLAGS += -Wall
CXXFLAGS += -O0 -ggdb
CXXFLAGS += -std=c++11

step5exe : step5.o 
	$(CXX) -o step5exe $(CXXFLAGS) $(LFLAGS) step5.o $(LIBS)

step5.o : step5.cpp
	$(CXX) -o step5.o $(CXXFLAGS) $(IFLAGS) -c step5.cpp 