IFLAGS += -I/usr/include

LFLAGS += -L/usr/lib/x86_64-linux-gnu

LIBS += -ldeal.ii
LIBS += -ltbb

CXX = g++
CXXFLAGS += -Wall
CXXFLAGS += -O0 -ggdb
CXXFLAGS += -std=c++11

step5exe : step6.o 
	$(CXX) -o step6exe $(LFLAGS) $(CXXFLAGS) step6.o $(LIBS)

step5.o : step6.cpp
	$(CXX) -o step6.o $(IFLAGS) $(CXXFLAGS) -c step6.cpp