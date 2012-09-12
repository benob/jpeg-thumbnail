CFLAGS=-O3
LDFLAGS=-ljpeg -lm
CC=gcc

jpeg-thumbnail:jpeg-thumbnail.o

clean:
	rm -rf *.o jpeg-thumbnail
