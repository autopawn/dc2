compile:
	rm -rf bin || true
	mkdir bin
	gcc -g -Wall src/*.c -lpthread -o bin/dc