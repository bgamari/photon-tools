CXXFLAGS = -ggdb

all : autocorr dump_pt2

autocorr : pt2.o
dump_pt2 : pt2.o

clean : 
	rm autocorr *.o

