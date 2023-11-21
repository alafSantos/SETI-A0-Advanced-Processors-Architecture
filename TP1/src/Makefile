CC = gcc
EXE = tp1
CFILES = tp1.c
CFLAGS = -DN=500 -DTYPE=float -O2

all:
	$(CC) $(CFILES) $(CFLAGS) -o $(EXE)

o0:
	$(CC) -O0 $(CFILES) -o $(EXE)

o1:
	$(CC) -O1 $(CFILES) -o $(EXE)

o2:
	$(CC) -O2 $(CFILES) -o $(EXE)

o3:
	$(CC) -O3 $(CFILES) -o $(EXE)

run:
	taskset -c 1 ./$(EXE)

clean:
	rm -rf $(EXE)