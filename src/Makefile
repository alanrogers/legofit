destination := $(HOME)/bin
global_destination := $(HOME)/group/bin
opt := -DNDEBUG -O3  -finline-functions  # For full optimization
#opt :=  -O0 -fno-inline-functions -rdynamic -DDEBUG     # For debugging
#opt := -O3  -finline-functions  # Optimization + debugging
#prof := -pg -rdynamic                    # For profiling
prof :=
incl := -I/usr/local/include -I/opt/local/include

targets := legosim legofit tabpat daf numcores
pytargets := diverg.py bootci.py flatfile.py
tests := xzeroin xbinary

CC := gcc

# Flags to determine the warning messages issued by the compiler
warn := \
 -Wall \
 -Wcast-align \
 -Wcast-qual \
 -Wmissing-declarations \
 -Wmissing-prototypes \
 -Wnested-externs \
 -Wpointer-arith \
 -Wstrict-prototypes \
 -Wno-unused-parameter \
 -Wno-unused-function \
 -Wshadow \
 -Wundef \
 -Wwrite-strings

CFLAGS := -g -std=gnu99 $(warn) $(incl) $(opt) $(prof) $(osargs)

lib := -L/usr/local/lib -lgsl -lgslcblas -lpthread -lm

.c.o:
	$(CC) $(CFLAGS) $(incl) -c -o ${@F}  $<

all : $(targets)

test : $(tests)
	-./xbinary
	-./xzeroin
	@echo "ALL UNIT TESTS WERE COMPLETED."

LEGOSIM := legosim.o patprob.o gptree.o binary.o jobqueue.o misc.o parse.o \
  branchtab.o popnodetab.o lblndx.o tokenizer.o parstore.o parkeyval.o \
  popnode.o gene.o dprintf.o rngseed.o dtnorm.o
legosim : $(LEGOSIM)
	$(CC) $(CFLAGS) -o $@ $(LEGOSIM) $(lib)

LEGOFIT := legofit.o patprob.o gptree.o binary.o jobqueue.o misc.o \
  parse.o branchtab.o popnodetab.o lblndx.o tokenizer.o parstore.o \
  parkeyval.o popnode.o gene.o cost.o diffev.o dprintf.o rngseed.o \
  simsched.o dtnorm.o
legofit : $(LEGOFIT)
	$(CC) $(CFLAGS) -o $@ $(LEGOFIT) $(lib)

TABPAT := tabpat.o misc.o binary.o lblndx.o parkeyval.o dafreader.o \
  tokenizer.o strint.o boot.o
tabpat : $(TABPAT)
	$(CC) $(CFLAGS) -o $@ $(TABPAT) $(lib)

DAF := daf.o misc.o
daf : $(DAF)
	$(CC) $(CFLAGS) -o $@ $(DAF) $(lib)

NUMCORES := numcores.o misc.o
numcores : $(NUMCORES)
	$(CC) $(CFLAGS) -o $@ $(NUMCORES) $(lib)

zeroin.o : zeroin.h

TRY := try.o
try : $(TRY)
	$(CC) $(CFLAGS) -o $@ $(TRY) $(lib)

# Make dependencies file
depend : *.c *.h
	echo '#Automatically generated dependency info' > depend
	$(CC) -MM $(incl) *.c >> depend

clean :
	rm -f *.a *.o *~ gmon.out *.tmp $(targets) $(tests) core.* vgcore.*

install : $(targets) $(pytargets)
	cp -p $(pytargets) $(destination)
	cp $(targets) $(destination)

ginstall : $(targets) $(pytargets)
	cp -p $(pytargets) $(global_destination)
	cp $(targets) $(global_destination)

include depend

.SUFFIXES:
.SUFFIXES: .c .o
.PHONY: clean

