.PHONY: bench-% %.o

CFLAGS := $(CFLAGS) -g

all: rtree-2.o

bench-%: rtree-%.o
	$(CC) $(CFLAGS) bench.c rtree-$*/$^ -o $@

%.o: */%.c
	$(CC) $(CFLAGS) -c $^ -o $*/$@

clean:
	$(RM) rtree-*/rtree-*.o bench-[0-9]*
