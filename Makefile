PROGS = autocorr dump_pt2 extract_pt2_timestamps
CXXFLAGS = -ggdb -std=c++0x ${INCLUDE}
LDFLAGS = ${LIBS}
CC=g++

all : ${PROGS} bin_photons

dump_pt2 : pt2.o
extract_pt2_timestamps : pt2.o

timetag_parse.so : timetag_parse.o
	gcc ${LDFLAGS} -shared -o $@ $+
timetag_parse.o : timetag_parse.c
	gcc -c -I/usr/include/python2.6 -fPIC ${CFLAGS} -o $@ $<

%.c : %.pyx
	cython $<

test : ${PROGS}
	./test_bin_photons.py

clean : 
	rm -f ${CLEAN} ${PROGS} *.o *.pyc

# For automatic header dependencies
.deps/%.d : %
	@mkdir -p .deps
	@makedepend  ${INCLUDES} -f - $< 2>/dev/null | sed 's,\($*\.o\)[ :]*,\1 $@ : ,g' >$@
SOURCES = $(wildcard *.cpp) $(wildcard *.c)
-include $(addprefix .deps/,$(addsuffix .d,$(SOURCES)))

