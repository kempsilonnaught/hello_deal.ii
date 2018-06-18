IFLAGS += -I/usr/include

LFLAGS += -L/usr/lib/x86_64-linux-gnu

LIBS += -ldeal.ii
LIBS += -ltbb

CXX = g++
CXXFLAGS += -Wall
CXXFLAGS += -O0 -ggdb
CXXFLAGS += -std=c++11

# CXX = clang
# CXXFLAGS += -std=c++11

step3exe : step3.o 
	$(CXX) -o step3exe $(CXXFLAGS) $(LFLAGS) step3.o $(LIBS)

step3.o : step3.cpp
	$(CXX) -o step3.o $(CXXFLAGS) $(IFLAGS) -c step3.cpp