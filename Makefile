.EXPORT_ALL_VARIABLES:

.PHONY: clean all

BIN_DIR = $(HOME)/bin
LIB_DIR = $(HOME)/lib

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTINC      := -I$(shell root-config --incdir)

COMMON_DIR = $(HOME)/common
BASELIBS  = -lm $(ROOTLIBS) $(ROOTGLIBS) -L$(LIB_DIR)
ALLIBS  = $(BASELIBS) -lCommandLineInterface

CPP             = g++
CFLAGS		= -Wall -Wno-long-long -g -O3 $(ROOTCFLAGS) -fPIC

INCLUDES        = -I./ -I$(COMMON_DIR) 
LFLAGS		= -g -fPIC
LIBS 		= $(ALLIBS)

CFLAGS += -Wl,--no-as-needed
LFLAGS += -Wl,--no-as-needed
DFLAGS += -Wl,--no-as-needed

O_FILES = Range.o

all: LineShape

LineShape: LineShape.cc $(O_FILES)
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) $< $(LIBS) $(O_FILES) -o $(BIN_DIR)/$@ 

%.o: %.cc %.hh
	@echo "Compiling $@"
	@$(CPP) $(CFLAGS) $(INCLUDES) -c $< -o $@ 

clean:
	@echo "Cleaning up"
	@rm -f *~ *.o

