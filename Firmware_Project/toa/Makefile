INCLUDE_DIRS =
LIB_DIRS =
CC=gcc

CDEFS=
CFLAGS1= -c -lrt -lm
CFLAGS2= -g -DUNITTESTS -lrt -lm

HFILES=toa.h
CFILES=toa.c

SRCS = ${HFILES} ${CFILES}

all: toa

clean:
	-rm -f *.o circ_buff

toa: toa.c toa.h
	$(CC) $(CFLAGS1)  toa.c

unit: toa.c toa.h
	$(CC) $(CFLAGS2)  toa.c

