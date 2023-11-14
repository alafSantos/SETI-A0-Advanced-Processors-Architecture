CC = gcc
EXE = tp1
CFILES = tp1.c

o0:
	$(CC) -O0 $(CFILES) -o $(EXE)

o1:
	$(CC) -O1 $(CFILES) -o $(EXE)

o2:
	$(CC) -O2 $(CFILES) -o $(EXE)

o3:
	$(CC) -O3 $(CFILES) -o $(EXE)

run:(EXE)
	./$(EXE)

clean:
	rm -rf $(EXE)