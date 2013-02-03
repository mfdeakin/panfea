
CC=gcc
CFLAGS=-g -Wall -std=c99
LIBS=

panfea: panfea.c
	$(CC) $(CFLAGS) $(LIBS) -o panfea panfea.c
