compile:
	rm -rf bin || true
	mkdir bin
	gcc -g -O2 -Wall src/*.c -lpthread -o bin/dc