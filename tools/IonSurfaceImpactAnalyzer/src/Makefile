FC = gfortran
CFLAGS = -O2 -Wall 

# flags
DBG := -O0 -g -fbacktrace -Wall -Wextra #-pedantic -fcheck=all -finit-real=snan -ffpe-trap=invalid -std=f2018
FFLAGS := $(DBG)
#
OBJS = share.o utilities.o parametersIO.o pusher.o boundaries.o
#
.PHONY: all clean
.SUFFIXES: .f90 .o

all: IonSurfaceImpactAnalysis

IonSurfaceImpactAnalysis: main.f90 $(OBJS)
	$(FC) $(FCFLAGS) $< $(OBJS) -o $@

.f90.o:
	$(FC) -c $(FCFLAGS) $<

%.o: %.mod

clean:
	$(RM) impact *.o *.mod
