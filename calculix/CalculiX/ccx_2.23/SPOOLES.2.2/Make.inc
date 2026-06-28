# SPOOLES.2.2 Make.inc for Windows (gcc via MSYS2)
CC = gcc
OPTLEVEL = -O2
CFLAGS = $(OPTLEVEL) -D_POSIX_C_SOURCE=199506L -fcommon -Wno-int-conversion -Wno-implicit-function-declaration -Wno-incompatible-pointer-types -Wno-return-type
LDFLAGS =
THREAD_LIBS = -lpthread
PURIFY =
PURIFY_GCC_VERSION =
AR = ar
ARFLAGS = rv
RANLIB = ranlib

.c.o :
	$(PURIFY) $(CC) -c $(CFLAGS) $<

.c.a :
	$(PURIFY) $(CC) -c $(CFLAGS) $<
	$(AR) $(ARFLAGS) $@ $*.o
	rm -f $*.o

MPI_INSTALL_DIR = /usr/local/mpich-1.0.13
MPI_LIB_PATH = -L$(MPI_INSTALL_DIR)/lib/solaris/ch_p4
MPI_LIBS = $(MPI_LIB_PATH) -D_REENTRANT -lmpi -lsocket -lnsl -lthread
MPI_INCLUDE_DIR = -I$(MPI_INSTALL_DIR)/include
