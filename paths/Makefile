all: paths

paths: main.o graph.o parser.o
	g++ main.o graph.o parser.o -o paths

main.o: main.cpp config.h
	g++ -c -W main.cpp 

paths.o: graph.cpp graph.h
	g++ -c -W graph.cpp

parser.o: parser.cpp parser.h
	g++ -c -W parser.cpp

clean:
	rm -rf *o paths
