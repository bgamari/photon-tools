CXXFLAGS = -ggdb

all : autocorr

autocorr : pt2.o

clean : 
	rm autocorr *.o

