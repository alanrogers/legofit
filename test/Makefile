vpath %.c ../src

opt :=  -O0 -fno-inline-functions -DDEBUG
prof :=
incl := -I/usr/local/include -I/opt/local/include -I../src
tests := xaddrparmap xbinary xbooma xboot xbranchtab xclic xcomb xdafreader \
  xdiffev \
  xdtnorm xerror xgene xgptree xhessian xidset xintpart \
  xjobqueue xlblndx xlinereader xmigoutcome xmisc xmpfrmatcoal xnodestore \
  xsimreader \
  xparam xparse xparseopf xparstore xpointbuff xpopnode xptrlst xptrptrmap \
  xptrqueue xptru32map \
  xptrvec \
  xstrptrmap xrafreader xrational xscrmreader xsampndx xsegment xsetpart \
  xsimsched xstate xstrdblqueue xstrint xstrparmap xtinyexpr xtokenizer \
  xu64u64map xuintqueue xwraparound xmatcoal

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
 -Wwrite-strings  \
 -Wno-nullability-completeness

# Abort on integer overflow or underflow
ifeq ($(shell uname -s), Darwin)
    # clang syntax
    #ubsan := -fsanitize=undefined,integer -fno-sanitize-recover=undefined,integer
    ubsan := -fsanitize=undefined,integer
else
    # gcc syntax
    #ubsan := -fsanitize=undefined -fno-sanitize-recover=undefined
    ubsan := -fsanitize=undefined
endif

CFLAGS := -g -std=gnu99 $(warn) $(ubsan) $(incl) $(opt) $(prof) $(osargs)

lib := -L/usr/local/lib -lgsl -lgslcblas -lpthread -lm

.c.o:
	$(CC) $(CFLAGS) -c -o ${@F}  $<

all : $(tests)

test : $(tests)
	-./xaddrparmap
	-./xbinary
	-./xbooma
	-./xboot
	-./xbranchtab
	-./xclic
	-./xcomb
	-./xdafreader
	-./xdiffev
	-./xdtnorm
	-./xerror
	-./xgene
	-./xgptree
	-./xhessian
	-./xidset
	-./xintpart
	-./xjobqueue
	-./xlblndx
	-./xlinereader
	-./xmatcoal
	-./xmigoutcome
	-./xmisc
	-./xmpfrmatcoal
	-./xnodestore
	-./xsimreader
	-./xparam
	-./xparse
	-./xparseopf
	-./xparstore
	-./xpointbuff
	-./xpopnode
	-./xptrlst
	-./xptrptrmap
	-./xptrqueue
	-./xptru32map
	-./xptrvec
	-./xrafreader
	-./xrational
	-./xsampndx
	-./xscrmreader
	-./xsegment
	-./xsetpart
	-./xsimsched
	-./xstate
	-./xstrdblqueue
	-./xstrint
	-./xstrparmap
	-./xstrptrmap
	-./xtinyexpr
	-./xtokenizer
	-./xu64u64map
	-./xuintqueue
	-./xwraparound
	@echo "ALL UNIT TESTS WERE COMPLETED."

XMATCOAL := xmatcoal.o matcoal.o misc.o rational.o mpfrmatcoal.o
xmatcoal : $(XMATCOAL)
	$(CC) $(CFLAGS) -o $@ $(XMATCOAL) $(lib) -lmpfr

XMPFRMATCOAL := xmpfrmatcoal.o mpfrmatcoal.o misc.o rational.o
xmpfrmatcoal : $(XMPFRMATCOAL)
	$(CC) $(CFLAGS) -o $@ $(XMPFRMATCOAL) $(lib) -lmpfr

XRATIONAL := xrational.o rational.o misc.o
xrational : $(XRATIONAL)
	$(CC) $(CFLAGS) -o $@ $(XRATIONAL) $(lib)

XBINARY := xbinary.o binary.o misc.o
xbinary : $(XBINARY)
	$(CC) $(CFLAGS) -o $@ $(XBINARY) $(lib)

XCOMB := xcomb.o comb.o misc.o
xcomb : $(XCOMB)
	$(CC) $(CFLAGS) -o $@ $(XCOMB) $(lib)

XSEGMENT := xsegment.o misc.o nodestore.o segment.o idset.o binary.o error.o \
    branchtab.o lblndx.o tokenizer.o matcoal.o comb.o rational.o setpart.o
xsegment : $(XSEGMENT)
	$(CC) $(CFLAGS) -o $@ $(XSEGMENT) $(lib)

XINTPART := xintpart.o intpart.o misc.o u64u64map.o binary.o
xintpart : $(XINTPART)
	$(CC) $(CFLAGS) -o $@ $(XINTPART) $(lib)

XMISC := xmisc.o misc.o
xmisc : $(XMISC)
	$(CC) $(CFLAGS) -o $@ $(XMISC) $(lib)

XPTRLST := xptrlst.o ptrlst.o misc.o
xptrlst : $(XPTRLST)
	$(CC) $(CFLAGS) -o $@ $(XPTRLST) $(lib)

XPTRQUEUE := xptrqueue.o ptrqueue.o misc.o
xptrqueue : $(XPTRQUEUE)
	$(CC) $(CFLAGS) -o $@ $(XPTRQUEUE) $(lib)

XPTRVEC := xptrvec.o ptrvec.o misc.o
xptrvec : $(XPTRVEC)
	$(CC) $(CFLAGS) -o $@ $(XPTRVEC) $(lib)

XMIGOUTCOME := xmigoutcome.o migoutcome.o misc.o
xmigoutcome : $(XMIGOUTCOME)
	$(CC) $(CFLAGS) -o $@ $(XMIGOUTCOME) $(lib)

XSETPART := xsetpart.o setpart.o misc.o u64u64map.o binary.o
xsetpart : $(XSETPART)
	$(CC) $(CFLAGS) -o $@ $(XSETPART) $(lib)

XLINEREADER := xlinereader.o linereader.o misc.o
xlinereader : $(XLINEREADER)
	$(CC) $(CFLAGS) -o $@ $(XLINEREADER) $(lib)

XWRAPAROUND := xwraparound.o wraparound.o
xwraparound : $(XWRAPAROUND)
	$(CC) $(CFLAGS) -o $@ $(XWRAPAROUND) $(lib)

XSTATE := xstate.o state.o misc.o
xstate : $(XSTATE)
	$(CC) $(CFLAGS) -o $@ $(XSTATE) $(lib)

XSTRPTRMAP := xstrptrmap.o strptrmap.o misc.o dprintf.o
xstrptrmap : $(XSTRPTRMAP)
	$(CC) $(CFLAGS) -o $@ $(XSTRPTRMAP) $(lib)

XDTNORM := xdtnorm.o dtnorm.o
xdtnorm : $(XDTNORM)
	$(CC) $(CFLAGS) -o $@ $(XDTNORM) $(lib)

XERROR := xerror.o error.o misc.o
xerror : $(XERROR)
	$(CC) $(CFLAGS) -o $@ $(XERROR) $(lib)

XPARAM := xparam.o misc.o param.o tinyexpr.o dtnorm.o 
xparam : $(XPARAM)
	$(CC) $(CFLAGS) -o $@ $(XPARAM) $(lib)

XTINYEXPR := xtinyexpr.o tinyexpr.o misc.o
xtinyexpr : $(XTINYEXPR)
	$(CC) $(CFLAGS) -o $@ $(XTINYEXPR) $(lib)

XTOKENIZER := xtokenizer.o tokenizer.o misc.o
xtokenizer : $(XTOKENIZER)
	$(CC) $(CFLAGS) -o $@ $(XTOKENIZER) $(lib)

XHESSIAN := xhessian.o misc.o hessian.o
xhessian : $(XHESSIAN)
	$(CC) $(CFLAGS) -o $@ $(XHESSIAN) $(lib)

XJOBQUEUE := xjobqueue.o jobqueue.o misc.o
xjobqueue : $(XJOBQUEUE)
	$(CC) $(CFLAGS) -o $@ $(XJOBQUEUE) $(lib)

XPOINTBUFF := xpointbuff.o misc.o pointbuff.o
xpointbuff : $(XPOINTBUFF)
	$(CC) $(CFLAGS) -o $@ $(XPOINTBUFF) $(lib)

XDIFFEV := xdiffev.o diffev.o binary.o lblndx.o jobqueue.o \
  simsched.o state.o pointbuff.o misc.o
xdiffev : $(XDIFFEV)
	$(CC) $(CFLAGS) -o $@ $(XDIFFEV) $(lib)

XDAFREADER := xdafreader.o dafreader.o misc.o lblndx.o tokenizer.o binary.o \
   strint.o
xdafreader : $(XDAFREADER)
	$(CC) $(CFLAGS) -o $@ $(XDAFREADER) $(lib)

XRAFREADER := xrafreader.o rafreader.o misc.o lblndx.o tokenizer.o binary.o \
   strint.o error.o
xrafreader : $(XRAFREADER)
	$(CC) $(CFLAGS) -o $@ $(XRAFREADER) $(lib)

xsampndx.o : sampndx.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/sampndx.c
XSAMPNDX := xsampndx.o misc.o
xsampndx : $(XSAMPNDX)
	$(CC) $(CFLAGS) -o $@ $(XSAMPNDX) $(lib)

xnodestore.o : nodestore.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/nodestore.c

XNODESTORE := xnodestore.o misc.o
xnodestore : $(XNODESTORE)
	$(CC) $(CFLAGS) -o $@ $(XNODESTORE) $(lib)

xbooma.o : booma.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/booma.c

XBOOMA := xbooma.o misc.o strdblqueue.o strint.o tokenizer.o
xbooma : $(XBOOMA)
	$(CC) $(CFLAGS) -o $@ $(XBOOMA) $(lib)

xscrmreader.o : scrmreader.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/scrmreader.c

XSCRMREADER := xscrmreader.o misc.o tokenizer.o uintqueue.o
xscrmreader : $(XSCRMREADER)
	$(CC) $(CFLAGS) -o $@ $(XSCRMREADER) $(lib)

xsimreader.o : simreader.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/simreader.c

XSIMREADER := xsimreader.o misc.o tokenizer.o
xsimreader : $(XSIMREADER)
	$(CC) $(CFLAGS) -o $@ $(XSIMREADER) $(lib)

XUINTQUEUE := xuintqueue.o uintqueue.o misc.o
xuintqueue : $(XUINTQUEUE)
	$(CC) $(CFLAGS) -o $@ $(XUINTQUEUE) $(lib)

xparse.o : parse.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/parse.c

XPARSE := xparse.o misc.o tokenizer.o gptree.o lblndx.o \
       branchtab.o network.o parstore.o popnode.o binary.o gene.o \
       dprintf.o dtnorm.o tinyexpr.o param.o strparmap.o addrparmap.o \
       sampndx.o strptrmap.o nodestore.o
xparse : $(XPARSE)
	$(CC) $(CFLAGS) -o $@ $(XPARSE) $(lib)

xboot.o : boot.c
	echo making xboot.o
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/boot.c

XBOOT := xboot.o misc.o binary.o
xboot : $(XBOOT)
	$(CC) $(CFLAGS) -o $@ $(XBOOT) $(lib)

xgene.o : gene.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/gene.c

XGENE := xgene.o misc.o binary.o branchtab.o tokenizer.o lblndx.o
xgene : $(XGENE)
	$(CC) $(CFLAGS) -o $@ $(XGENE) $(lib)

xpopnode.o : ../src/popnode.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/popnode.c

XPOPNODE := xpopnode.o dtnorm.o misc.o gene.o parstore.o tinyexpr.o \
   param.o branchtab.o nodestore.o ptrptrmap.o ptrqueue.o strint.o \
   binary.o
xpopnode : $(XPOPNODE)
	$(CC) $(CFLAGS) -o $@ $(XPOPNODE) $(lib)

XSIMSCHED := xsimsched.o simsched.o misc.o
xsimsched : $(XSIMSCHED)
	$(CC) $(CFLAGS) -o $@ $(XSIMSCHED) $(lib)

xgptree.o : gptree.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/gptree.c

XGPTREE := xgptree.o addrparmap.o misc.o branchtab.o parstore.o lblndx.o \
        tokenizer.o strptrmap.o gene.o param.o popnode.o binary.o \
        dprintf.o dtnorm.o parse.o strparmap.o tinyexpr.o sampndx.o \
        nodestore.o network.o ptrptrmap.o ptrqueue.o strint.o error.o
xgptree : $(XGPTREE)
	$(CC) $(CFLAGS) -o $@ $(XGPTREE) $(lib)

XPARSTORE := xparstore.o misc.o parstore.o binary.o ptrptrmap.o ptrqueue.o \
        tinyexpr.o param.o strint.o \
	dtnorm.o 
xparstore : $(XPARSTORE)
	$(CC) $(CFLAGS) -o $@ $(XPARSTORE) $(lib)

# On the Mac, "make" fails to compile this if I name the object
# file "xbranchtab". I have no clue why.
xbrtab.o : branchtab.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/branchtab.c

XBRANCHTAB := xbrtab.o misc.o binary.o 
xbranchtab : $(XBRANCHTAB)
	$(CC) $(CFLAGS) -o $@ $(XBRANCHTAB) $(lib)

xlblndx.o : lblndx.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/lblndx.c

XLBLNDX := xlblndx.o misc.o binary.o
xlblndx : $(XLBLNDX)
	$(CC) $(CFLAGS) -o $@ $(XLBLNDX) $(lib)

XPARSEOPF := xparseopf.o parseopf.o branchtab.o misc.o lblndx.o tokenizer.o \
   binary.o
xparseopf : $(XPARSEOPF)
	$(CC) $(CFLAGS) -o $@ $(XPARSEOPF) $(lib)

xstrint.o : strint.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/strint.c

XSTRINT := xstrint.o misc.o binary.o lblndx.o
xstrint : $(XSTRINT)
	$(CC) $(CFLAGS) -o $@ $(XSTRINT) $(lib)

xstrparmap.o : strparmap.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/strparmap.c

XSTRPARMAP := xstrparmap.o misc.o param.o tinyexpr.o dtnorm.o gptree.o \
   parstore.o gene.o lblndx.o popnode.o branchtab.o tokenizer.o \
   parse.o binary.o strptrmap.o dprintf.o addrparmap.o sampndx.o
xstrparmap : $(XSTRPARMAP)
	$(CC) $(CFLAGS) -o $@ $(XSTRPARMAP) $(lib)

xaddrparmap.o : addrparmap.c
	$(CC) $(CFLAGS) -c -DTEST -o $@ ../src/addrparmap.c

XADDRPARMAP := xaddrparmap.o misc.o param.o tinyexpr.o dtnorm.o gptree.o \
   parstore.o gene.o lblndx.o popnode.o branchtab.o strparmap.o \
   tokenizer.o parse.o binary.o strptrmap.o dprintf.o sampndx.o \
   nodestore.o network.o
xaddrparmap : $(XADDRPARMAP)
	$(CC) $(CFLAGS) -o $@ $(XADDRPARMAP) $(lib)

u64u64map.o : ../src/u64u64map.c ../src/u64u64map.h ../src/hashmap.h \
   ../src/hashmap.c

XU64U64MAP := xu64u64map.o misc.o binary.o u64u64map.o
xu64u64map : $(XU64U64MAP)
	$(CC) $(CFLAGS) -o $@ $(XU64U64MAP) $(lib)

ptru32map.o : ../src/ptru32map.c ../src/ptru32map.h ../src/hashmap.h \
   ../src/hashmap.c

XPTRU32MAP := xptru32map.o misc.o binary.o ptru32map.o
xptru32map : $(XPTRU32MAP)
	$(CC) $(CFLAGS) -o $@ $(XPTRU32MAP) $(lib)

ptrptrmap.o : ../src/ptrptrmap.c ../src/ptrptrmap.h ../src/hashmap.h \
   ../src/hashmap.c

XPTRPTRMAP := xptrptrmap.o misc.o binary.o ptrptrmap.o
xptrptrmap : $(XPTRPTRMAP)
	$(CC) $(CFLAGS) -o $@ $(XPTRPTRMAP) $(lib)

XSTRDBLQUEUE := xstrdblqueue.o strdblqueue.o misc.o
xstrdblqueue : $(XSTRDBLQUEUE)
	$(CC) $(CFLAGS) -o $@ $(XSTRDBLQUEUE) $(lib)

XIDSET := xidset.o idset.o misc.o binary.o migoutcome.o
xidset : $(XIDSET)
	$(CC) $(CFLAGS) -o $@ $(XIDSET) $(lib)

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
