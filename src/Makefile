destination := $(HOME)/bin
#opt := -DNDEBUG -O3  -finline-functions  # For full optimization
#opt := -O3  -finline-functions  # Optimization + debugging
opt :=  -O0 -fno-inline-functions -rdynamic -DDEBUG     # For debugging
#prof := -pg -rdynamic                    # For profiling
prof :=
incl := -I/usr/local/include -I/opt/local/include

targets := lego legofit
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

LEGO := lego.o patprob.o gptree.o binary.o jobqueue.o misc.o parse.o \
  branchtab.o hashtab.o lblndx.o tokenizer.o parstore.o parkeyval.o \
  popnode.o gene.o dprintf.o
lego : $(LEGO)
	$(CC) $(CFLAGS) -o $@ $(LEGO) $(lib)

LEGOFIT := legofit.o patprob.o gptree.o binary.o jobqueue.o misc.o \
  parse.o branchtab.o hashtab.o lblndx.o tokenizer.o parstore.o \
  parkeyval.o popnode.o gene.o cost.o diffev.o dprintf.o
legofit : $(LEGOFIT)
	$(CC) $(CFLAGS) -o $@ $(LEGOFIT) $(lib)

TABPAT := tabpat.o misc.o binary.o lblndx.o parkeyval.o vcfreader.o \
  tokenizer.o
tabpat : $(TABPAT)
	$(CC) $(CFLAGS) -o $@ $(TABPAT) $(lib)

zeroin.o : zeroin.h

# Make dependencies file
depend : *.c *.h
	echo '#Automatically generated dependency info' > depend
	$(CC) -MM $(incl) *.c >> depend

clean :
	rm -f *.a *.o *~ gmon.out *.tmp $(targets) $(tests) core.* vgcore.*

install : $(targets)
	cp diverg.py $(destination)
	cp $(targets) $(destination)

include depend

.SUFFIXES:
.SUFFIXES: .c .o
.PHONY: clean

