#Make file for inlapardiso interface
OBJECTS = interface.o ../libpardiso500-GNU461-X86-64.so 

#gcc mypardiso_sym.c -g -o mypardiso_sym.exe -L/home/litvina/PARDISO libpardiso500-GNU461-X86-64.so -llapack -#lblas  -lgfortran -fopenmp -lpthread -lm


prefix = /home/litvina/PARDISO/inlapardiso

exec_prefixH = ${prefixH}
#CFLAGS = -p -O2 -Wall -pedantic
CFLAGS = -pg
LDFLAGS =  
LIBS =  -llapack -lblas -lgfortran -lm -fopenmp -lpthread
CC = gcc


projects: inlapardiso

clean:
	rm -f inlapardiso
	rm -f *~ *.ps--
	rm -f *.o
	
interface.o: interface.h
interface.o: interface.c
	$(CC) $(CFLAGS) -c interface.c -o $@


inlapardiso: $(OBJECTS)
	$(CC) $(CFLAGS)  $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@
