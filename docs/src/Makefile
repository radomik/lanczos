CC=/home/darek/bin/mex
LDFLAGS=$(shell pkg-config --libs lapack)

all:
	$(CC) -v mrrr.c $(LDFLAGS)
