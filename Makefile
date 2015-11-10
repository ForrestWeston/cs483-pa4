CC=mpicc
CFLAGS=-g -lrt -lm -Wall -Werror
EXE=-o run

all:
		$(CC) $(CFLAGS) $(EXE) main.c matrix.c prn.c
