PROGRAM=	main

OBJECTS=	main.o \
		AFWL.o

CC=	g++
FLAGS=	-Wall -g -C

LIBS=	-lm

INCLUDES=	

$(PROGRAM):	$(OBJECTS)
	$(CC) $(INCLUDES) $(FLAGS) -o $(PROGRAM) $(OBJECTS) $(LIBS)

main.o:	main.cpp\
	AFWL.cpp
	$(CC) $(INCLUDES) $(FLAGS) -c $<

AFWL.o:	AFWL.cpp
	$(CC) $(INCLUDES) $(FLAGS) -c $<


clean:
	rm *.o $(PROGRAM)
