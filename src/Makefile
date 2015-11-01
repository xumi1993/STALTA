CC = gcc
LIBS = -lsac -lsacio -L/usr/local/sac/lib
INCLUDE = -I/usr/local/sac/include

all: stalta envelope

BIN = ../bin

stalta:  stalta.c sac.h sacio.c smooth.c
	${CC} -o ${BIN}/$@ $^ ${INCLUDE} ${LIBS}

envelope: envelope.c sacio.c sac.h
	${CC} -o ${BIN}/$@ $^ ${INCLUDE} ${LIBS}

clean:
	rm ${BIN}/* *.o
