vpath %.c ../src

#opt := -DNDEBUG -O3  -finline-functions  # For full optimization
opt :=  -O0 -fno-inline-functions -DDEBUG     # For debugging
#prof := -pg -rdynamic                    # For profiling
prof :=
incl := -I/usr/local/include -I/opt/local/include -I../src

tests := xbinary xhashtab xparse xgene xgptree xbranchtab xparstore \
  xparkeyval xlblndx xpopnode xjobqueue xdiffev xdafreader xstrint \
  xboot xsimsched

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
	-./xboot
	-./xhashtab
	-./xgene
	-./xparse
	-./xpopnode
	-./xgptree
	-./xparkeyval
	-./xparstore
	-./xsimsched
	-./xstrint
	-./xlblndx
	-./xdiffev
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

XHASHTAB := xhashtab.o hashtab.o misc.o gptree.o branchtab.o parstore.o \
   parkeyval.o tokenizer.o binary.o lblndx.o popnode.o gene.o parse.o
xhashtab : $(XHASHTAB)
	$(CC) $(CFLAGS) -o $@ $(XHASHTAB) $(lib)

XJOBQUEUE := xjobqueue.o jobqueue.o
xjobqueue : $(XJOBQUEUE)
	$(CC) $(CFLAGS) -o $@ $(XJOBQUEUE) $(lib)

XDIFFEV := xdiffev.o diffev.o misc.o binary.o lblndx.o jobqueue.o parkeyval.o
xdiffev : $(XDIFFEV)
	$(CC) $(CFLAGS) -o $@ $(XDIFFEV) $(lib)

XDAFREADER := xdafreader.o dafreader.o misc.o lblndx.o tokenizer.o binary.o \
   parkeyval.o strint.o
xdafreader : $(XDAFREADER)
	$(CC) $(CFLAGS) -o $@ $(XDAFREADER) $(lib)

xparse.o : parse.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/parse.c

XPARSE := xparse.o hashtab.o misc.o tokenizer.o gptree.o lblndx.o \
       branchtab.o parstore.o parkeyval.o popnode.o binary.o gene.o
xparse : $(XPARSE)
	$(CC) $(CFLAGS) -o $@ $(XPARSE) $(lib)

xgene.o : gene.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/gene.c

XGENE := xgene.o branchtab.o misc.o binary.o tokenizer.o lblndx.o parkeyval.o
xgene : $(XGENE)
	$(CC) $(CFLAGS) -o $@ $(XGENE) $(lib)

XBOOT := xboot.o misc.o boot.o binary.o
xboot : $(XBOOT)
	$(CC) $(CFLAGS) -o $@ $(XBOOT) $(lib)

xpopnode.o : popnode.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/popnode.c

XPOPNODE := xpopnode.o misc.o gene.o branchtab.o binary.o lblndx.o tokenizer.o \
   parkeyval.o
xpopnode : $(XPOPNODE)
	$(CC) $(CFLAGS) -o $@ $(XPOPNODE) $(lib)

XSIMSCHED := xsimsched.o simsched.o misc.o
xsimsched : $(XSIMSCHED)
	$(CC) $(CFLAGS) -o $@ $(XSIMSCHED) $(lib)

xgptree.o : gptree.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/gptree.c

XGPTREE := xgptree.o misc.o branchtab.o parstore.o parse.o lblndx.o \
        parkeyval.o tokenizer.o hashtab.o gene.o popnode.o binary.o
xgptree : $(XGPTREE)
	$(CC) $(CFLAGS) -o $@ $(XGPTREE) $(lib)

XPARSTORE := xparstore.o misc.o parstore.o parkeyval.o binary.o lblndx.o
xparstore : $(XPARSTORE)
	$(CC) $(CFLAGS) -o $@ $(XPARSTORE) $(lib)

XPARKEYVAL := xparkeyval.o parkeyval.o misc.o binary.o lblndx.o
xparkeyval : $(XPARKEYVAL)
	$(CC) $(CFLAGS) -o $@ $(XPARKEYVAL) $(lib)

xlblndx.o : lblndx.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/lblndx.c

XLBLNDX := xlblndx.o misc.o parkeyval.o binary.o
xlblndx : $(XLBLNDX)
	$(CC) $(CFLAGS) -o $@ $(XLBLNDX) $(lib)

xbranchtab.o : branchtab.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/branchtab.c

XBRANCHTAB := xbranchtab.o gptree.o misc.o binary.o parstore.o popnode.o gene.o \
   lblndx.o parse.o parkeyval.o tokenizer.o hashtab.o
xbranchtab : $(XBRANCHTAB)
	$(CC) $(CFLAGS) -o $@ $(XBRANCHTAB) $(lib)

xstrint.o : strint.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/strint.c

XSTRINT := xstrint.o misc.o binary.o lblndx.o parkeyval.o
xstrint : $(XSTRINT)
	$(CC) $(CFLAGS) -o $@ $(XSTRINT) $(lib)

Q := q.o misc.o gptree.o hashtab.o branchtab.o binary.o lblndx.o parstore.o \
   popnode.o parse.o gene.o tokenizer.o parkeyval.o
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

