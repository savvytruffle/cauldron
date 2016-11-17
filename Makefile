# -*- makefile -*-

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
EXT = linux
GCC_FLAGS1_gsl = -I/astro/users/aleezah/Documents/cauldron/src/include -fPIC -Wl,-Bsymbolic-functions -c -O3
GCC_FLAGS2_gsl = -L/astro/users/aleezah/Documents/cauldron/src/lib -lgsl -lgslcblas -lm -shared -O3 -Wl,-Bsymbolic-functions,-soname,helpers_linux.so
endif

GCC = gcc

.PHONY: all
.SILENT: all

all:
	echo "Compiling C source code for helper funcs..."
	${GCC} ${GCC_FLAGS1_gsl} helpers.c
	echo "Generating shared library for helper funcs..."
	gcc ${GCC_FLAGS2_gsl} -o helpers_${EXT}.so helpers.o -lc
	rm helpers.o
	echo "Install successful."
	
