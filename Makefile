CC = gcc

CFLAGS = -g3 -O2 -Wall -Wextra -Werror

.PHONY: all

all: seir_d.out

seir_d.out: seir_d.c
	$(CC) $(CFLAGS) -o $@ $^ -lm

.PHONY: clean
clean:
	find . -type f -name "*.so" -exec rm -f {} \;
	find . -type f -name "*.a" -exec rm -f {} \;
	find . -type f -name "*.out" -exec rm -f {} \;
	find . -type f -name "*.o" -exec rm -f {} \;
	find . -type f -name "*.d" -exec rm -f {} \;

rebuild: clean all

run:
	./seir_d.out
