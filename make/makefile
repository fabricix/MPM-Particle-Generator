# compiler to be used
CC = g++

# compiler flags
#  -g : adds debugging information to the executable file
#  -Wall : turns on compiler warnings
CFLAGS = -Wall -g

# include directory
INC = -I../inc

# dependencies
OBJ = main.o model.o io.o matrix.o paraview.o utils.o

# make the program
MPM-Particle-Generator: $(OBJ)
	$(CC) $(CFLAGS) $(INC) -o MPM-Particle-Generator $(OBJ)

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

# this target is not associated with any file
.PHONY: clean

# clean
clean :
	rm MPM-Particle-Generator $(OBJ)
