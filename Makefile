#############################################################################
# Makefile for building: fpa
#############################################################################

MAKEFILE      = Makefile

####### Compiler, tools and options

CC            = gcc
CXX           = g++
CFLAGS        = -m64 -pipe -O2 -Wall -std=c++11
CXXFLAGS      = -m64 -pipe -O2 -Wall -std=c++11
INCPATH       = -I. 
LINK          = g++
LFLAGS        = -m64 -Wl,-O1
LIBS          = -lboost_program_options $(SUBLIBS)
DEL_FILE      = rm -f

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = fpa.cpp 
OBJECTS       = fpa.o

DESTDIR       = 
TARGET        = fpa


first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: Makefile $(TARGET)

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

clean: compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core

compiler_clean: 

####### Compile

fpa.o: fpa.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fpa.o fpa.cpp

####### Install

install:  FORCE

uninstall:  FORCE

FORCE:

