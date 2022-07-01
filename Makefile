CC=g++

runner.o: MatrixImpl.cpp
	$(CC)  Runner.cpp -o Runner.o

clean:
	rm -r Runner.o 
