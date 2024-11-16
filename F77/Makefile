PROGRAM=	main

OBJECTS=	main.o \
		air.o \
		matm62.o \
		peak.o \
		scalkt.o \
		shock.o \
		speed.o \
		terror.o \
		timar1.o \
		well.o \
		wfdrmt.o \
		wfdzr.o \
		wfpkod.o \
		wfpkop.o \
		wfpkv.o \
		wfpr.o \
		wfprmt.o \
		wfvrmt.o \
		wfvzr.o \
		wfzr.o

F77=	pgf77
FLAGS=	-Wall -Msave -Ktrap=divz,denorm,fp,ovf -C -g
LIBS=	
INCLUDES=	

F77=	gfortran
FLAGS=	-O3 -Wall -Wunused-parameter -fno-range-check -std=legacy -ffpe-trap=zero,denormal,invalid,overflow -C -g
LIBS=	
INCLUDES=	

$(PROGRAM):	$(OBJECTS)
	$(F77) $(INCLUDES) $(FLAGS) -o $(PROGRAM) $(OBJECTS) $(LIBS)

main.o:	main.for\
	air.f \
	matm62.f \
	peak.f \
	scalkt.f \
	shock.f \
	speed.f \
	terror.f \
	timar1.f \
	well.f \
	wfdrmt.f \
	wfdzr.f \
	wfpkod.f \
	wfpkop.f \
	wfpkv.f \
	wfpr.f \
	wfprmt.f \
	wfvrmt.f \
	wfvzr.f \
	wfzr.f
	$(F77) $(INCLUDES) $(FLAGS) -c $<

air.o:	air.f
	$(F77) $(INCLUDES) $(FLAGS) -c $<

matm62.o:	matm62.f
	$(F77) $(INCLUDES) $(FLAGS) -c $<

peak.o:	peak.f\
	wfrt.h 
	$(F77) $(INCLUDES) $(FLAGS) -c $<

scalkt.o:	scalkt.f
	$(F77) $(INCLUDES) $(FLAGS) -c $<

shock.o:	shock.f
	$(F77) $(INCLUDES) $(FLAGS) -c $<

speed.o:	speed.f
	$(F77) $(INCLUDES) $(FLAGS) -c $<

terror.o:	terror.f
	$(F77) $(INCLUDES) $(FLAGS) -c $<

timar1.o:	timar1.f
	$(F77) $(INCLUDES) $(FLAGS) -c $<

well.o:	well.f\
	wfrt.h 
	$(F77) $(INCLUDES) $(FLAGS) -c $<

wfdrmt.o:	wfdrmt.f\
	wfrt.h 
	$(F77) $(INCLUDES) $(FLAGS) -c $<

wfdzr.o:	wfdzr.f
	$(F77) $(INCLUDES) $(FLAGS) -c $<

wfpkod.o:	wfpkod.f\
	wfrt.h 
	$(F77) $(INCLUDES) $(FLAGS) -c $<

wfpkop.o:	wfpkop.f\
	wfrt.h 
	$(F77) $(INCLUDES) $(FLAGS) -c $<

wfpkv.o:	wfpkv.f\
	wfrt.h 
	$(F77) $(INCLUDES) $(FLAGS) -c $<

wfpr.o:	wfpr.f
	$(F77) $(INCLUDES) $(FLAGS) -c $<

wfprmt.o:	wfprmt.f\
	wfrt.h 
	$(F77) $(INCLUDES) $(FLAGS) -c $<

wfvrmt.o:	wfvrmt.f\
	wfrt.h 
	$(F77) $(INCLUDES) $(FLAGS) -c $<

wfvzr.o:	wfvzr.f
	$(F77) $(INCLUDES) $(FLAGS) -c $<

wfzr.o:	wfzr.f
	$(F77) $(INCLUDES) $(FLAGS) -c $<


clean:
	rm *.o $(PROGRAM)
