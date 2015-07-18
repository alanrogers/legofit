vpath %.c ../src

#opt := -DNDEBUG -O3  -finline-functions  # For full optimization
opt :=  -O0 -fno-inline-functions -DDEBUG     # For debugging
#prof := -pg -rdynamic                    # For profiling
prof :=
incl := -I/usr/local/include -I/opt/local/include -I../src

tests := xzeroin xbinary xhashtab xmktree

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
	$(CC) $(CFLAGS) -c -o ${@F}  $<

all : $(targets)

test : $(tests)
	-./xbinary
	-./xzeroin
	@echo "ALL UNIT TESTS WERE COMPLETED."

zeroin.o : zeroin.h

xzeroin.o : zeroin.c zeroin.h
	$(CC) $(CFLAGS) -c -DTEST -o $@ zeroin.c

XZEROIN := xzeroin.o
xzeroin : $(XZEROIN)
	$(CC) $(CFLAGS) -o $@ $(XZEROIN) $(lib)

XBINARY := xbinary.o binary.o
xbinary : $(XBINARY)
	$(CC) $(CFLAGS) -o $@ $(XBINARY) $(lib)

XHASHTAB := xhashtab.o hashtab.o misc.o
xhashtab : $(XHASHTAB)
	$(CC) $(CFLAGS) -o $@ $(XHASHTAB) $(lib)

xmktree.o : mktree.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/mktree.c

XMKTREE := xmktree.o hashtab.o misc.o tokenizer.o gptree.o
xmktree : $(XMKTREE)
	$(CC) $(CFLAGS) -o $@ $(XMKTREE) $(lib)

# Make dependencies file
depend : *.c 
	echo '#Automatically generated dependency info' > depend
	$(CC) -MM $(incl) *.c >> depend

clean :
	rm -f *.a *.o *~ gmon.out *.tmp $(targets) $(tests) core.* vgcore.*

include depend

.SUFFIXES:
.SUFFIXES: .c .o
.PHONY: clean

