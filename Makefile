MPICC = mpic++

NP = -n 1

FILE = ./graph/web-BerkStan.txt

OUT = paths.txt

SOURCE = 0

#MACHINE = -machinefile machines
MACHINE =

all: clean proj run

MyQueue.o: MyQueue.h MyQueue.cpp
	$(MPICC) -c MyQueue.cpp

short.o: GetPath.cpp MyQueue.h
	$(MPICC) -c GetPath.cpp -o short.o

proj: short.o MyQueue.o
	$(MPICC) short.o MyQueue.o -o proj

run:
	mpirun $(MACHINE) $(NP) ./proj $(FILE) $(OUT) $(SOURCE)

clean:
	rm -f *o proj paths.txt
