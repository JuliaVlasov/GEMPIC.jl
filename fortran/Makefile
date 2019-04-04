PROG =	VM.exe

SRCS =	zone.f90 initialisation.f90 quietstart.f90 particules.f90 \
	villasenor.f90 maxwell.f90 diagno.f90 main.f90

OBJS =	zone.o initialisation.o quietstart.o particules.o \
	villasenor.o maxwell.o diagno.o main.o

LIBS =	-llapack

F90      = gfortran
F90FLAGS =  -O3
LDFLAGS  = ${F90FLAGS}

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod fort.* 

debug:	
	make F90FLAGS=-g FFLAGS=-g

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $*.f90
.mod.o:
	$(F90) $(F90FLAGS) -c $*.f90
.f90.ps:
	a2ps $*.f90 -o $*.ps

main.o : main.f90 initialisation.o particules.o villasenor.o \
	maxwell.o diagno.o zone.o 
initialisation.o : initialisation.f90 zone.o
particules.o : particules.f90 quietstart.o zone.o
quietstart.o : quietstart.f90 zone.o
villasenor.o : villasenor.f90 zone.o
maxwell.o : maxwell.f90 zone.o
diagno.o : diagno.f90 zone.o
zone.o : zone.f90
