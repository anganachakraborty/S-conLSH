CC=g++
CFLAGS=-O3 -c -Wall
DFLAGS=-g -Wall
LDFLAGS= -lz -lm -lpthread

KSWSOURCE=ksw.c 
#ALGNSOURCES=main.cpp aligner.cpp  readfl.cpp hash.cpp form.cpp btree.cpp conLSH.cpp
INDSOURCES=whash.cpp genhash.cpp formh.cpp conLSH.cpp readfl.cpp aligner.cpp

INDOBJECTS=$(INDSOURCES:.cpp=.o) $(KSWSOURCE:.c=.o)
#ALGNOBJECTS=$(ALGNSOURCES:.cpp=.o) $(KSWSOURCE:.c=.o)

ALIGNER=S-conLSH
#ALIGNER=S-conLSH-aligner

all: $(INDSOURCES) $(KSWSOURCE) $(ALIGNER)
	
#$(ALIGNER): $(ALGNOBJECTS)
	#$(CC)  $(ALGNOBJECTS) -o $@ $(LDFLAGS)

$(ALIGNER): $(INDOBJECTS)
	$(CC)  $(INDOBJECTS) -o $@ $(LDFLAGS)

debug:
	$(CC) $(DFLAGS) $(KSWSOURCE) $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
.c.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm -rf *.o $(ALIGNER) a.out

 
