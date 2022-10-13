# compiler
CC = g++

# compiler flags
#  -g : adds debugging information to the executable file
#  -Wall : turn on compiler warnings

# Debug
CFLAGS = -Wall -g

# Release
# CFLAGS = -Wall

# include directory
INC=-I../inc

# dependencies
objects = main.o\
	model.o\
	io.o\
	matrix.o\
	paraview.o\
	utils.o
	
# make the program
MPM-Particle-Generator: $(objects)
	$(CC) $(CFLAGS) $(INC) -o MPM-Particle-Generator $(objects)

# compile dependencies
main.o : ../src/main.cpp ../inc/io.h ../inc/model.h
	$(CC) $(CFLAGS) $(INC) -c ../src/main.cpp

model.o : ../src/model.cpp ../inc/utils.h ../inc/matrix.h
	$(CC) $(CFLAGS) $(INC) -c ../src/model.cpp

io.o : ../src/io.cpp ../inc/model.h ../inc/paraview.h
	$(CC) $(CFLAGS) $(INC) -c ../src/io.cpp

matrix.o: ../src/matrix.cpp
	$(CC) $(CFLAGS) $(INC) -c ../src/matrix.cpp

paraview.o: ../src/paraview.cpp ../inc/io.h ../inc/model.h
	$(CC) $(CFLAGS) $(INC) -c ../src/paraview.cpp

utils.o: ../src/utils.cpp ../inc/matrix.h
	$(CC) $(CFLAGS) $(INC) -c ../src/utils.cpp

# clean
clean :
	rm MPM-Particle-Generator $(objects)
