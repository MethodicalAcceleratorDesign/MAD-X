all: name_mangler_c

name_mangler_c: name_mangler_c.o name_mangler.o
	g++ name_mangler_c.o  name_mangler.o -o name_mangler_c

name_mangler_c.o: name_mangler_c.c
	gcc name_mangler_c.c -c

name_mangler.o: name_mangler.cc
	g++ -c -o name_mangler.o name_mangler.cc

clean:
	rm -f name_mangler_c.o name_mangler.o name_mangler_c
