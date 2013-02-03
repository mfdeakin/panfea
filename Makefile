
CC=gcc
CFLAGS=-g -Wall -std=c99
LIBS=-lGL -lglut

panfea: panfea.c
	$(CC) $(CFLAGS) $(LIBS) -o panfea panfea.c
