CC = g++
CFLAGS = -c -std=c++11 -g
OBJ = main.o particle.o collisions.o field.o
LIBS = -lsfml-graphics -lsfml-window -lsfml-system 
EXE = sims

all: $(OBJ)
	$(CC) $(OBJ) -o $(EXE) $(LIBS)

main.o: main.cxx $(DEPS)
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm $(OBJ) $(EXE)
