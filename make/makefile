# compiler to be used
CC = g++

# compiler flags
#  -g : adds debugging information to the executable file
#  -Wall : turn on compiler warnings

# Debug
# CFLAGS = -Wall -g

# Release
# CFLAGS = -Wall

# include directory
INC = -I../inc

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
main.o : ../src/main.cpp
	$(CC) $(CFLAGS) $(INC) -c ../src/main.cpp

model.o : ../src/model.cpp
	$(CC) $(CFLAGS) $(INC) -c ../src/model.cpp

io.o : ../src/io.cpp
	$(CC) $(CFLAGS) $(INC) -c ../src/io.cpp

matrix.o: ../src/matrix.cpp
	$(CC) $(CFLAGS) $(INC) -c ../src/matrix.cpp

paraview.o: ../src/paraview.cpp
	$(CC) $(CFLAGS) $(INC) -c ../src/paraview.cpp

utils.o: ../src/utils.cpp
	$(CC) $(CFLAGS) $(INC) -c ../src/utils.cpp

# this target is not associated with any file
.PHONY: clean

# clean
clean :
	rm MPM-Particle-Generator $(objects)
