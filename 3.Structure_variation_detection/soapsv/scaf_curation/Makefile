SHELL         = /bin/sh
PROG          = scaf_restrictor
DEBUG         = NO
PROFILE       = NO
PTHREADS      = YES
STATIC        = YES
CC            = g++
CXX           = g++
DEBUG_FLAGS   = -g3 -Wall -O2
PROFILE_FLAGS = -fprofile-arcs -ftest-coverage -pg 
RELEASE_FLAGS = -O3 -fomit-frame-pointer
STATIC_FLAGS  = -static
DFLAGS        = -DMAKE_TIME=\""`date`"\"
LIBS          = -lz -lboost_regex
INC           = -I/FULLPATH/boost/include/
LIB           = -L/FULLPATH/boost/stage/lib/
#TARBALL_EXCLUDE = "*.(o,gz,zip)"
#ZIP_EXCLUDE     = *.o *.gz *.zip

ifeq (YES, $(DEBUG))
        CFLAGS   = $(DEBUG_FLAGS)
        DFLAGS  += -DDEBUG
#        PTHREADS = NO
else
        CFLAGS   = $(RELEASE_FLAGS)
endif

ifeq (YES, $(PTHREADS))
        LIBS   +=  -pthread
        #DFLAGS +=  -DPTHREADS
endif

ifeq (YES, $(STATIC))
        CFLAGS += $(STATIC_FLAGS)
endif

ifeq (YES, $(PROFILE))
       DFLAGS += $(PROFILE_FLAGS)
endif

OBJ = main.o func.o func2.o pileup.o general.o gzstream.o threadmanager.o
.SUFFIX:
.SUFFIX: .cpp .o

.cpp.o:
	$(CXX) -c $(CFLAGS) $(DFLAGS) $(INC) $(LIB) $< -o $@

all: $(PROG)

$(PROG): $(OBJ)
	$(CXX) $(CFLAGS) $(DFLAGS) $(OBJ) $(INC) $(LIB) -o $@ $(LIBS)

pileup.o:pileup.h
main.o:main.h general.h func.h func2.h
func.o:general.h func.h
func2.o:general.h func.h func2.h main.h threadmanager.h
general.o: general.h
gzstream.o: gzstream.h
threadmanager.o:threadmanager.h

clean:
	rm -f *.o $(PROG)
