# My 
SYSTEM=WORKING
# If your
# SYSTEM=WINDOWS/DOS
# it should also work (well, apart from the OS)
# If you want to use this on an android tablet or smartphone your
# SYSTEM=WORKING/ARM

# Compile options:
# Debugging
# CFLAGS= -g -pedantic -Wall -Wno-unknown-pragmas
# use optimizations
# CFLAGS= -O3 -pedantic -Wall -Wno-unknown-pragmas
# optimize more
CFLAGS= -O3 -ffast-math -pedantic -Wall -Wno-unknown-pragmas
# parallelization using openmp
# CFLAGS= -O3 -pedantic -Wall -fopenmp
# full optimization and openmp
# CFLAGS= -O3 -ffast-math -pedantic -Wall -fopenmp 

# Linker Options:
# Note to self: when static linking the order of object specification matters
# The symbols are parsed from the right. Make sure libm.a comes last!
# LFLAGS= -static -lm
LFLAGS= -lm
# LFLAGS= -lm -fopenmp

VERSION = 1.2

SRC=sourcefield.c sourcewarp.c main.c parse.c utils.c mod_bessel_01.c scd2.c
HDR=sourcefield.h sourcewarp.h parsedef.h parse.h main.h utils.h mod_bessel_01.h scd2.h
OBJ=sourcefield.o sourcewarp.o parse.o main.o utils.o mod_bessel_01.o scd2.o

ifeq ($(SYSTEM), WINDOWS/DOS)
# mingw cross-compiler
CC = i686-pc-mingw32-gcc
TARGET = SourceField.exe 
else
ifeq ($(SYSTEM), WORKING/ARM)
CC=arm-gp2x-linux-gcc
TARGET = SourceField
else
CC=gcc
TARGET = SourceField
endif
endif

test: SourceField
	./test.sh $(TARGET)
SourceField: $(OBJ) keywords.h Makefile
	$(CC)  $(OBJ) $(LFLAGS) -o $(TARGET)
keywords.h: keywordtable.txt
	./KeywordTable2h.sh > keywords.h
utils.o: utils.c main.h Makefile
	$(CC) $(CFLAGS) -c -DVERSION=\"$(VERSION)\" utils.c
main.o: main.h parse.h Makefile main.c defaults.h keywords.h
parse.o: main.h utils.h sourcefield.h parsedef.h parse.h Makefile parse.c
sourcefield.o: main.h utils.h sourcefield.h Makefile sourcefield.c mod_bessel_01.h
sourcewarp.o: main.h utils.h sourcewarp.h Makefile sourcewarp.c scd2.h
mod_bessel_01.o: mod_bessel_01.c
scd2.o: scd2.c
clean:
	rm *.o $(TARGET) keywords.h
