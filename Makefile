PROGS = autocorr dump_pt2 extract_pt2_timestamps bin_photons
CXXFLAGS = -ggdb -std=c++0x ${INCLUDE}
LDFLAGS = ${LIBS}
CC=g++

all : ${PROGS}

dump_pt2 : pt2.o
extract_pt2_timestamps : pt2.o
bin_photons : _bin_photons.so
_bin_photons.so : _bin_photons.o
	gcc ${LDFLAGS} -shared -o $@ $+
_bin_photons.o : _bin_photons.c
	gcc -c -I/usr/include/python2.6 -fPIC ${CFLAGS} -o $@ $<

%.c : %.pyx
	cython $<

test : ${PROGS}
	./test_bin_photons.py

clean : 
	rm -f ${PROGS} *.o

# For automatic header dependencies
.deps/%.d : %
	@mkdir -p .deps
	@makedepend  ${INCLUDES} -f - $< 2>/dev/null | sed 's,\($*\.o\)[ :]*,\1 $@ : ,g' >$@
SOURCES = $(wildcard *.cpp) $(wildcard *.c)
-include $(addprefix .deps/,$(addsuffix .d,$(SOURCES)))

