compile:
	rm -rf bin || true
	mkdir bin
	gcc -g -Wall src/*.c -D THREADS=0 -o bin/dc