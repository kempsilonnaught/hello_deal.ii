# Oops
IFLAGS += -I/usr/include

LFLAGS += -L/usr/lib/x86_64-linux-gnu

LIBS += -ldeal.ii

CXXFLAGS += -Wall
CXXFLAGS += -O0 -ggdb

progname : objects
	g++ $(CXXFLAGS) -o progname objects $(LFLAGS) $(LIBS)

object : source
	g++ $(CXXFLAGS) $(IFLAGS) -c -o object source 