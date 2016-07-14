# Definicion del directorio donde estan los ejecutables
#
XDIR	=	$(HOME)/pvm3/bin/$(PVM_ARCH)
#
# Compilador de fortran, opciones y librerias
F77     =       gfortran
FFLAGS  =       -O -w
FLIBS	=	-lfpvm3 -lpvm3
GLIB    =       -lgpvm3

# Path absoluto donde se encuentran las librerias anteriores
LFLAGS	=	-L$(PVM_ROOT)/lib/$(PVM_ARCH)


# Path donde se encuentran los archivos include
INCL    =       -I/usr/include


all: norma jacobi

norma: norma.f
	$(F77) -c $(FFLAGS) $(INCL) norma.f

jacobi: jacobi.f norma.o
	$(F77) $(FFLAGS) -o $@ $@.f norma.o $(LFLAGS) $(FLIBS) $(GLIB) $(INCL)
	mv $@ $(XDIR)
