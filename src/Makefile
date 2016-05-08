destination := $(HOME)/bin
#opt := -DNDEBUG -O3  -finline-functions  # For full optimization
#opt := -O3  -finline-functions  # Optimization + debugging
opt :=  -O0 -fno-inline-functions -DDEBUG     # For debugging
#prof := -pg -rdynamic                    # For profiling
prof :=
incl := -I/usr/local/include -I/opt/local/include

target := lego
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

all : $(target)

test : $(tests)
	-./xbinary
	-./xzeroin
	@echo "ALL UNIT TESTS WERE COMPLETED."

LEGO := lego.o gptree.o binary.o jobqueue.o misc.o branchtab.o \
  hashtab.o sampndx.o parse.o tokenizer.o parstore.o parkeyval.o
lego : $(LEGO)
	$(CC) $(CFLAGS) -o $@ $(LEGO) $(lib)

LEGO2 := lego2.o patprob.o gptree.o binary.o jobqueue.o misc.o parse.o \
  branchtab.o hashtab.o sampndx.o tokenizer.o parstore.o parkeyval.o
lego2 : $(LEGO2)
	$(CC) $(CFLAGS) -o $@ $(LEGO2) $(lib)

zeroin.o : zeroin.h

# Make dependencies file
depend : *.c *.h
	echo '#Automatically generated dependency info' > depend
	$(CC) -MM $(incl) *.c >> depend

clean :
	rm -f *.a *.o *~ gmon.out *.tmp $(targets) $(tests) core.* vgcore.*

install : $(target)
	cp diverg.py $(destination)
	cp $(target) $(destination)

include depend

.SUFFIXES:
.SUFFIXES: .c .o
.PHONY: clean

