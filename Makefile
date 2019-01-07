.EXPORT_ALL_VARIABLES:

.SUFFIXES:

.PHONY: clean all

# := is only evaluated once

SHELL = /bin/sh

NAME = SortData

LIB_DIR = $(HOME)/lib

# ROOT specfic stuff
#ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)
ROOTINC    := -I$(shell root-config --incdir)

# GRSIsort specfic stuff
#GRSICFLAGS := $(shell $(GRSISYS)/bin/grsi-config --cflags)
GRSILIBS   := $(shell $(GRSISYS)/bin/grsi-config --all-libs)
GRSIINC    := -I$(GRSISYS)/include

# define any directories containing header files other than /usr/include
INCLUDES = -I.

# define the C compiler to use
CC = gcc

# define the C++ compiler to use
CXX = g++

# define any compile-time flags
CPPFLAGS = $(ROOTINC) $(GRSIINC) $(INCLUDES) -fPIC
CXXFLAGS = -std=gnu++0x -pedantic -Wall -Wno-long-long -g -O3
LDFLAGS  = -g -fPIC

# define library paths in addition to /usr/lib
# if I wanted to include libraries not in /usr/lib I'd specify
# their path using -Lpath, something like:
LDLIBS    = -L$(LIB_DIR) -Wl,-rpath,/opt/gcc/lib64 $(ROOTLIBS) $(GRSILIBS) -lm
LOADLIBES  = \
	ReadConfig.o

# -------------------- implicit rules --------------------
# n.o made from n.c by 		$(CC) -c $(CPPFLAGS) $(CFLAGS)
# n.o made from n.cc by 	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS)
# n made from n.o by 		$(CC) $(LDFLAGS) n.o $(LOADLIBES) $(LDLIBS)

# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above
#

# -------------------- rules --------------------

all: $(NAME)
	@echo Done

# -------------------- pattern rules --------------------
# this rule sets the name of the .cc file at the beginning of the line (easier to find)
# it is a suffix replacement rule for building .o's from .cc's uses automatic variables 
# $<: the name of the prerequisite of the rule (a .cc file) and $@: the name of the target 
# of the rule (a .o file) (see the gnu make manual section about automatic variables)

%.o: %.cc %.hh
	$(CXX) $< -c $(CPPFLAGS) $(CXXFLAGS) -o $@

# -------------------- default rule for executables --------------------

%: %.cc $(LOADLIBES)
	$(CXX) $< $(CXXFLAGS) $(CPPFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

# -------------------- clean --------------------

clean:
	rm -f $(NAME) *.o *~
