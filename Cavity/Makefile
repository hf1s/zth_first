CC = g++ -O2 
LINK = -lz 
TARGET = cavity
all: $(TARGET)
$(TARGET): main.o atom.o grid.o parameter.o pocket.o protein.o utils.o ligand.o
	$(CC) -o $(TARGET) main.o atom.o grid.o parameter.o pocket.o protein.o utils.o ligand.o $(LINK)
main.o: main.cpp utils.h parameter.h protein.h pocket.h ligand.h cmdline.h
	$(CC) -c main.cpp -o main.o
pocket.o : pocket.cpp pocket.h grid.h protein.h parameter.h ligand.h
	$(CC) -c pocket.cpp -o pocket.o
grid.o : grid.cpp grid.h parameter.h protein.h
	$(CC) -c grid.cpp -o grid.o
protein.o : protein.cpp protein.h utils.h parameter.h
	$(CC) -c protein.cpp -o protein.o
parameter.o : parameter.cpp parameter.h utils.h cmdline.h
	$(CC) -c parameter.cpp -o parameter.o
ligand.o : ligand.cpp ligand.h utils.h atom.h
	$(CC) -c ligand.cpp -o ligand.o
atom.o : atom.cpp atom.h utils.h
	$(CC) -c atom.cpp -o atom.o
utils.o : utils.cpp utils.h
	$(CC) -c utils.cpp -o utils.o

install:
	install $(TARGET) ../bin/$(TARGET)	

.PHONY : clean
clean :
	rm *.o
