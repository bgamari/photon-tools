CXXFLAGS = -ggdb
PROGS = autocorr_pt2 dump_pt2

all : ${PROGS}

autocorr_pt2 : pt2.o
dump_pt2 : pt2.o

clean : 
	rm ${PROGS} *.o

