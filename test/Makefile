vpath %.c ../src

opt :=  -O0 -fno-inline-functions -DDEBUG
prof :=
incl := -I/usr/local/include -I/opt/local/include -I../src
tests := xbinary xbma xboot xbranchtab xclic xdafreader xrafreader xdiffev \
  xgene xgptree xhessian xpopnodetab xjobqueue xlblndx xllrbtree xparkeyval \
  xparse xparstore xpopnode xscrmreader xsimsched xstrdblqueue xstrint \
  xstrdblmap xdtnorm xmisc xerror xstate xwraparound xpointbuff

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
	-./xbma
	-./xboot
	-./xbranchtab
	-./xclic
	-./xdafreader
	-./xrafreader
	-./xdiffev
	-./xdtnorm
	-./xerror
	-./xgene
	-./xgptree
	-./xhessian
	-./xjobqueue
	-./xlblndx
	-./xllrbtree
	-./xmisc
	-./xparkeyval
	-./xparse
	-./xparstore
	-./xpointbuff
	-./xpopnode
	-./xpopnodetab
	-./xscrmreader
	-./xsimsched
	-./xstate
	-./xstrint
	-./xstrdblmap
	-./xstrdblqueue
	-./xwraparound
	@echo "ALL UNIT TESTS WERE COMPLETED."

XBINARY := xbinary.o binary.o
xbinary : $(XBINARY)
	$(CC) $(CFLAGS) -o $@ $(XBINARY) $(lib)

XMISC := xmisc.o misc.o
xmisc : $(XMISC)
	$(CC) $(CFLAGS) -o $@ $(XMISC) $(lib)

XWRAPAROUND := xwraparound.o wraparound.o
xwraparound : $(XWRAPAROUND)
	$(CC) $(CFLAGS) -o $@ $(XWRAPAROUND) $(lib)

XSTATE := xstate.o state.o misc.o
xstate : $(XSTATE)
	$(CC) $(CFLAGS) -o $@ $(XSTATE) $(lib)

XPOPNODETAB := xpopnodetab.o popnodetab.o misc.o popnode.o gene.o \
   branchtab.o lblndx.o tokenizer.o dtnorm.o binary.o \
   parkeyval.o parstore.o tinyexpr.o
xpopnodetab : $(XPOPNODETAB)
	$(CC) $(CFLAGS) -o $@ $(XPOPNODETAB) $(lib)

XDTNORM := xdtnorm.o dtnorm.o
xdtnorm : $(XDTNORM)
	$(CC) $(CFLAGS) -o $@ $(XDTNORM) $(lib)

XERROR := xerror.o error.o misc.o
xerror : $(XERROR)
	$(CC) $(CFLAGS) -o $@ $(XERROR) $(lib)

XHESSIAN := xhessian.o misc.o hessian.o
xhessian : $(XHESSIAN)
	$(CC) $(CFLAGS) -o $@ $(XHESSIAN) $(lib)

XJOBQUEUE := xjobqueue.o jobqueue.o misc.o
xjobqueue : $(XJOBQUEUE)
	$(CC) $(CFLAGS) -o $@ $(XJOBQUEUE) $(lib)

XPOINTBUFF := xpointbuff.o misc.o pointbuff.o
xpointbuff : $(XPOINTBUFF)
	$(CC) $(CFLAGS) -o $@ $(XPOINTBUFF) $(lib)

XDIFFEV := xdiffev.o diffev.o binary.o lblndx.o jobqueue.o parkeyval.o \
  simsched.o state.o misc.o
xdiffev : $(XDIFFEV)
	$(CC) $(CFLAGS) -o $@ $(XDIFFEV) $(lib)

XDAFREADER := xdafreader.o dafreader.o misc.o lblndx.o tokenizer.o binary.o \
   parkeyval.o strint.o
xdafreader : $(XDAFREADER)
	$(CC) $(CFLAGS) -o $@ $(XDAFREADER) $(lib)

XRAFREADER := xrafreader.o rafreader.o misc.o lblndx.o tokenizer.o binary.o \
   parkeyval.o strint.o error.o
xrafreader : $(XRAFREADER)
	$(CC) $(CFLAGS) -o $@ $(XRAFREADER) $(lib)

XCLIC := xclic.o misc.o hessian.o strdblqueue.o
xclic : $(XCLIC)
	$(CC) $(CFLAGS) -o $@ $(XCLIC) $(lib)

xbma.o : bma.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/bma.c

XBMA := xbma.o misc.o strdblqueue.o strint.o tokenizer.o
xbma : $(XBMA)
	$(CC) $(CFLAGS) -o $@ $(XBMA) $(lib)

xscrmreader.o : scrmreader.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/scrmreader.c

XSCRMREADER := xscrmreader.o misc.o tokenizer.o
xscrmreader : $(XSCRMREADER)
	$(CC) $(CFLAGS) -o $@ $(XSCRMREADER) $(lib)

xparse.o : parse.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/parse.c

XPARSE := xparse.o popnodetab.o misc.o tokenizer.o gptree.o lblndx.o \
       branchtab.o parstore.o parkeyval.o popnode.o binary.o gene.o \
       dprintf.o dtnorm.o tinyexpr.o
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

XPOPNODE := xpopnode.o misc.o gene.o branchtab.o binary.o lblndx.o \
   tokenizer.o parkeyval.o dtnorm.o parstore.o tinyexpr.o
xpopnode : $(XPOPNODE)
	$(CC) $(CFLAGS) -o $@ $(XPOPNODE) $(lib)

XSIMSCHED := xsimsched.o simsched.o misc.o
xsimsched : $(XSIMSCHED)
	$(CC) $(CFLAGS) -o $@ $(XSIMSCHED) $(lib)

xgptree.o : gptree.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/gptree.c

XGPTREE := xgptree.o misc.o branchtab.o parstore.o lblndx.o \
        parkeyval.o tokenizer.o popnodetab.o gene.o popnode.o binary.o \
        dprintf.o dtnorm.o parse.o tinyexpr.o
xgptree : $(XGPTREE)
	$(CC) $(CFLAGS) -o $@ $(XGPTREE) $(lib)

XPARSTORE := xparstore.o misc.o parstore.o parkeyval.o binary.o lblndx.o \
        dtnorm.o tinyexpr.o
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

xllrbtree.o : llrbtree.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/llrbtree.c

XLLRBTREE := xllrbtree.o
xllrbtree : $(XLLRBTREE)
	$(CC) $(CFLAGS) -o $@ $(XLLRBTREE) $(lib)

xbranchtab.o : branchtab.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/branchtab.c

XBRANCHTAB := xbranchtab.o gptree.o misc.o binary.o parstore.o popnode.o \
   gene.o lblndx.o parse.o parkeyval.o tokenizer.o popnodetab.o \
   dprintf.o dtnorm.o tinyexpr.o
xbranchtab : $(XBRANCHTAB)
	$(CC) $(CFLAGS) -o $@ $(XBRANCHTAB) $(lib)

xstrint.o : strint.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/strint.c

XSTRINT := xstrint.o misc.o binary.o lblndx.o parkeyval.o
xstrint : $(XSTRINT)
	$(CC) $(CFLAGS) -o $@ $(XSTRINT) $(lib)

xstrdblmap.o : strdblmap.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/strdblmap.c

XSTRDBLMAP := xstrdblmap.o misc.o
xstrdblmap : $(XSTRDBLMAP)
	$(CC) $(CFLAGS) -o $@ $(XSTRDBLMAP) $(lib)

XSTRDBLQUEUE := xstrdblqueue.o strdblqueue.o misc.o
xstrdblqueue : $(XSTRDBLQUEUE)
	$(CC) $(CFLAGS) -o $@ $(XSTRDBLQUEUE) $(lib)

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
