
CC=mpicc

TARGET=src

all:
	$(CC) -o $(TARGET) src.c -lm
	g++ ~/UGP/allocator/src/allocator_improved.cpp -o ~/UGP/allocator/src/allocator.out

clean:
	rm $(TARGET)	
	rm ~/UGP/allocator/src/allocator.out
