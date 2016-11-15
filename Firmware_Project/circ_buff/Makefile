INCLUDE_DIRS =
LIB_DIRS =
CC=gcc

CDEFS=
CFLAGS1= -c
CFLAGS2= -g -DUNITTESTS

HFILES=circ_buff.h
CFILES=circ_buff.c

SRCS = ${HFILES} ${CFILES}

all: circ_buff

clean:
	-rm -f *.o circ_buff

circ_buff: circ_buff.c
	$(CC) $(CFLAGS1)  circ_buff.c

unit: circ_buff.c
	$(CC) $(CFLAGS2)  circ_buff.c

