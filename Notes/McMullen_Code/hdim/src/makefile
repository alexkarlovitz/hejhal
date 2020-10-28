
CC     = gcc
CFLAGS = -g
BIN	= ../bin

hdim:	hdim.o julia.o klein.o io.o eigen.o \
	markov.o cx.o ps.o alloc.o hdim.h
	if(test ! -x ${BIN}) then (mkdir ${BIN}) fi
	gcc $(CFLAGS) hdim.o julia.o klein.o io.o \
	eigen.o markov.o cx.o ps.o alloc.o -lm -o ${BIN}/hdim

hdim.o:		hdim.h
io.o:		hdim.h
eigen.o:	hdim.h
markov.o:	hdim.h
klein.o:	hdim.h
julia.o:	hdim.h
ps.o:		hdim.h
alloc.o:	hdim.h	
