# Public domain.

d2src = digest2.c d2model.c d2math.c d2mpc.c

digest2: $(d2src) digest2.h d2model.h
	gcc -o digest2 -std=c11 -pthread $(d2src) -lm -static
