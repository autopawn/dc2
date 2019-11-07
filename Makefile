compile:
	rm -rf bin || true
	mkdir bin
	gcc -g -Wall src/*.c -lpthread -D THREADS=4 -o bin/dc