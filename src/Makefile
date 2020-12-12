destination := $(HOME)/bin
git_version := "$(shell git describe --dirty --always --tags)"
global_destination := $(HOME)/group/bin
opt := -DNDEBUG -O2  -finline-functions  # For full optimization
#opt := -O0 -g -fno-inline-functions -rdynamic -DDEBUG     # For debugging
#opt := -O2  -g -finline-functions  # Optimization + debugging
#prof := -pg -g -rdynamic                    # For profiling under gcc
#prof := --fprofile-generate -g              # For profiling under clang
prof :=
incl := -I/usr/local/include -I/opt/local/include

targets := legosim legofit tabpat sitepat scrmpat daf raf numcores \
     joinraf raf2daf bepe booma pclgo resid simpat ms2sim clic rafdaf
pytargets := diverg.py bootci.py flatfile.py axt2raf.py flatpat.py legotree.py
tests := xzeroin xbinary

CC := gcc

# Flags to determine the warning messages issued by the compiler
warn := \
 -Wall \
 -Wcast-align \
 -Wcast-qual \
 -Wmissing-declarations \
 -Wmissing-prototypes \
 -Wpointer-arith \
 -Wstrict-prototypes \
 -Wno-unused-parameter \
 -Wno-unused-function \
 -Wshadow \
 -Wwrite-strings

# MacOS only
ifeq ($(shell uname -s), Darwin)
  warn +=  -Wno-nullability-completeness
endif

# Abort on integer overflow or underflow
ubsan =
#ifeq ($(shell uname -s), Darwin)
#    # clang syntax
#    ubsan += -fsanitize=undefined,integer \
#              -fno-sanitize-recover=undefined,integer
#else
#    # gcc syntax
#    ubsan += -fsanitize=undefined -fno-sanitize-recover=undefined
#endif

CFLAGS := -std=gnu99 $(ubsan) $(warn) $(incl) $(opt) $(prof) \
 -DGIT_VERSION=\"$(git_version)\"

lib := -L/usr/local/lib -lgsl -lgslcblas -lpthread -lm

.c.o:
	$(CC) $(CFLAGS) $(incl) -c -o ${@F}  $<

all : $(targets)

test : $(tests)
	-./xbinary
	-./xzeroin
	@echo "ALL UNIT TESTS WERE COMPLETED."

LEGOSIM := legosim.o binary.o branchtab.o dprintf.o dtnorm.o \
   error.o externs.o eventlst.o gene.o gptree.o jobqueue.o lblndx.o network.o \
   misc.o param.o parse.o parstore.o patprob.o popnode.o rngseed.o sampndx.o \
   strptrmap.o tinyexpr.o tokenizer.o ptrptrmap.o ptrqueue.o strint.o \
   segment.o idset.o idsetset.o matcoal.o ptrlst.o comb.o \
   setpart.o rational.o u64u64map.o mctree.o u64i64map.o
legosim : $(LEGOSIM)
	$(CC) $(CFLAGS) -o $@ $(LEGOSIM) $(lib)

LEGOFIT := legofit.o binary.o branchtab.o cost.o diffev.o \
   dprintf.o dtnorm.o error.o externs.o gene.o gptree.o jobqueue.o lblndx.o \
   misc.o network.o param.o parse.o parseopf.o parstore.o patprob.o \
   pointbuff.o popnode.o rngseed.o sampndx.o simsched.o state.o \
   strptrmap.o tinyexpr.o tokenizer.o ptrptrmap.o ptrqueue.o strint.o \
   mctree.o matcoal.o segment.o eventlst.o idset.o idsetset.o ptrlst.o \
   rational.o comb.o u64i64map.o setpart.o u64u64map.o
legofit : $(LEGOFIT)
	$(CC) $(CFLAGS) -o $@ $(LEGOFIT) $(lib)

PCLGO := pclgo.o misc.o strdblqueue.o gptree.o parstore.o param.o popnode.o \
  gene.o branchtab.o lblndx.o dtnorm.o dprintf.o externs.o \
  tinyexpr.o tokenizer.o parse.o binary.o strptrmap.o sampndx.o network.o \
  error.o ptrptrmap.o ptrqueue.o strint.o mctree.o segment.o eventlst.o \
  idsetset.o matcoal.o ptrlst.o comb.o idset.o u64i64map.o rational.o \
  setpart.o u64u64map.o
pclgo : $(PCLGO)
	$(CC) $(CFLAGS) -o $@ $(PCLGO) $(lib)

BEPE := bepe.o misc.o strdblqueue.o
bepe : $(BEPE)
	$(CC) $(CFLAGS) -o $@ $(BEPE) $(lib)

CLIC := clic.o hessian.o strdblqueue.o misc.o
clic : $(CLIC)
	$(CC) $(CFLAGS) -o $@ $(CLIC) $(lib)

RESID := resid.o misc.o strdblqueue.o branchtab.o lblndx.o binary.o \
  tokenizer.o
resid : $(RESID)
	$(CC) $(CFLAGS) -o $@ $(RESID) $(lib)

BOOMA := booma.o misc.o strdblqueue.o strint.o tokenizer.o
booma : $(BOOMA)
	$(CC) $(CFLAGS) -o $@ $(BOOMA) $(lib)

TABPAT := tabpat.o misc.o binary.o lblndx.o dafreader.o \
  tokenizer.o strint.o boot.o error.o
tabpat : $(TABPAT)
	$(CC) $(CFLAGS) -o $@ $(TABPAT) $(lib)

SITEPAT := sitepat.o misc.o binary.o lblndx.o rafreader.o \
  tokenizer.o strint.o boot.o error.o
sitepat : $(SITEPAT)
	$(CC) $(CFLAGS) -o $@ $(SITEPAT) $(lib)

RAFDAF := rafdaf.o misc.o binary.o lblndx.o rafreader.o \
  tokenizer.o error.o strint.o
rafdaf : $(RAFDAF)
	$(CC) $(CFLAGS) -o $@ $(RAFDAF) $(lib)

JOINRAF := joinraf.o rafreader.o error.o tokenizer.o misc.o binary.o
joinraf : $(JOINRAF)
	$(CC) $(CFLAGS) -o $@ $(JOINRAF) $(lib)

RAF2DAF := raf2daf.o rafreader.o error.o tokenizer.o misc.o binary.o
raf2daf : $(RAF2DAF)
	$(CC) $(CFLAGS) -o $@ $(RAF2DAF) $(lib)

SCRMPAT := scrmpat.o misc.o binary.o lblndx.o scrmreader.o \
  tokenizer.o error.o uintqueue.o
scrmpat : $(SCRMPAT)
	$(CC) $(CFLAGS) -o $@ $(SCRMPAT) $(lib)

SIMPAT := simpat.o misc.o binary.o lblndx.o simreader.o \
  tokenizer.o error.o
simpat : $(SIMPAT)
	$(CC) $(CFLAGS) -o $@ $(SIMPAT) $(lib)

MS2SIM := ms2sim.o binary.o misc.o lblndx.o linereader.o tokenizer.o
ms2sim : $(MS2SIM)
	$(CC) $(CFLAGS) -o $@ $(MS2SIM) $(lib)

DAF := daf.o misc.o
daf : $(DAF)
	$(CC) $(CFLAGS) -o $@ $(DAF) $(lib)

RAF := raf.o misc.o
raf : $(RAF)
	$(CC) $(CFLAGS) -o $@ $(RAF) $(lib)

NUMCORES := numcores.o misc.o
numcores : $(NUMCORES)
	$(CC) $(CFLAGS) -o $@ $(NUMCORES) $(lib)

zeroin.o : zeroin.h

TRY := try.o binary.o
try : $(TRY)
	$(CC) $(CFLAGS) -o $@ $(TRY) $(lib)

ptrptrmap.o : ptrptrmap.c ptrptrmap.h hashmap.h hashmap.src

ptru32map.o : ptru32map.c ptru32map.h hashmap.h hashmap.src

u64u64map.o : u64u64map.c u64u64map.h hashmap.h hashmap.src

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
