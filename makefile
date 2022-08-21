IDIR =./include
BDIR = ./bin
SRC = src
CC=g++ -Wall -g -std=c++1z
CFLAGS=-I$(IDIR)

ODIR=./obj
LDIR =./lib

_OBJ = FactorGraph.o SPSolver.o main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all:$(BDIR)/main

$(ODIR)/%.o: $(SRC)/%.cc $(IDIR)/*.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(BDIR)/main: $(OBJ)
	$(CC) -o $@ $^ 


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(BDIR)/* 