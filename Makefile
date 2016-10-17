# -*- makefile -*-

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
EXT = linux
GCC_FLAGS1_gsl = -I$(GSL_PATH)/include -fPIC -Wl,-Bsymbolic-functions -c -O3
GCC_FLAGS2_gsl = -L$(GSL_PATH)/lib -lgsl -lgslcblas -lm -shared -O3 -Wl,-Bsymbolic-functions,-soname,helpers_linux.so
GCC_FLAGS1 = -fPIC -Wl,-Bsymbolic-functions -c -O3
GCC_FLAGS2 = -lm -shared -O3 -Wl,-Bsymbolic-functions,-soname,gravitation_linux.so
endif
ifeq ($(UNAME_S),Darwin)
EXT = mac
GCC_FLAGS1_gsl = -I$(GSL_PATH)/include -fPIC -c
GCC_FLAGS2_gsl = -L$(GSL_PATH)/lib -lgsl -lgslcblas -shared -Wl,-install_name,helpers_mac.so
GCC_FLAGS1 = -fPIC -c
GCC_FLAGS2 = -lm -shared -Wl,-install_name,gravitation_mac.so
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
	
