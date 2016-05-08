vpath %.c ../src

#opt := -DNDEBUG -O3  -finline-functions  # For full optimization
opt :=  -O0 -fno-inline-functions -DDEBUG     # For debugging
#prof := -pg -rdynamic                    # For profiling
prof :=
incl := -I/usr/local/include -I/opt/local/include -I../src

#tests := xzeroin xbinary xhashtab xparse xgptree xbranchtab
tests := xbinary xhashtab xparse xgptree xbranchtab xparstore \
  xparkeyval xlblndx

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

all : $(tests)

test : $(tests)
	-./xbinary
	-./xhashtab
	-./xparse
	-./xgptree
	-./xparkeyval
	-./xparstore
	-./xlblndx
#	-./xzeroin
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

XHASHTAB := xhashtab.o hashtab.o misc.o gptree.o branchtab.o
xhashtab : $(XHASHTAB)
	$(CC) $(CFLAGS) -o $@ $(XHASHTAB) $(lib)

xparse.o : parse.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/parse.c

XPARSE := xparse.o hashtab.o misc.o tokenizer.o gptree.o lblndx.o \
       branchtab.o parstore.o parkeyval.o
xparse : $(XPARSE)
	$(CC) $(CFLAGS) -o $@ $(XPARSE) $(lib)

xgptree.o : gptree.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/gptree.c

XGPTREE := xgptree.o misc.o branchtab.o
xgptree : $(XGPTREE)
	$(CC) $(CFLAGS) -o $@ $(XGPTREE) $(lib)

XPARSTORE := xparstore.o misc.o parstore.o parkeyval.o
xparstore : $(XPARSTORE)
	$(CC) $(CFLAGS) -o $@ $(XPARSTORE) $(lib)

XPARKEYVAL := xparkeyval.o parkeyval.o misc.o
xparkeyval : $(XPARKEYVAL)
	$(CC) $(CFLAGS) -o $@ $(XPARKEYVAL) $(lib)

xlblndx.o : lblndx.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/lblndx.c

XLBLNDX := xlblndx.o misc.o
xlblndx : $(XLBLNDX)
	$(CC) $(CFLAGS) -o $@ $(XLBLNDX) $(lib)

xbranchtab.o : branchtab.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/branchtab.c

XBRANCHTAB := xbranchtab.o gptree.o misc.o binary.o
xbranchtab : $(XBRANCHTAB)
	$(CC) $(CFLAGS) -o $@ $(XBRANCHTAB) $(lib)

Q := q.o misc.o gptree.o hashtab.o branchtab.o
q : $(Q)
	$(CC) $(CFLAGS) -o $@ $(Q) $(lib)

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

