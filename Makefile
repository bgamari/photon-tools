CXXFLAGS = -ggdb -std=c++0x
PROGS = autocorr dump_pt2 extract_pt2_timestamps

all : ${PROGS}

dump_pt2 : pt2.o
extract_pt2_timestamps : pt2.o

clean : 
	rm ${PROGS} *.o

