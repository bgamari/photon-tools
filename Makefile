PROGS = bin_photons
CXXFLAGS = -ggdb -std=c++0x ${INCLUDE}
LDFLAGS = ${LIBS}
CC=g++
DESTDIR ?= /usr/local

include pt2/Makefile

all : ${PROGS}

timetag_parse.so : timetag_parse.o
	gcc ${LDFLAGS} -shared -o $@ $+
timetag_parse.o : timetag_parse.c
	gcc -c -I/usr/include/python2.6 -fPIC ${CFLAGS} -o $@ $<

%.c : %.pyx
	cython $<

test : ${PROGS}
	./test_bin_photons.py

.PHONY : install
install : ${PROGS}
	cp ${PROGS} ${DESTDIR}/bin

.PHONY : clean
clean : 
	rm -f ${CLEAN} ${PROGS} *.o *.pyc

# For automatic header dependencies
.deps/%.d : %
	@mkdir -p .deps
	@makedepend  ${INCLUDES} -f - $< 2>/dev/null | sed 's,\($*\.o\)[ :]*,\1 $@ : ,g' >$@
SOURCES = $(wildcard *.cpp) $(wildcard *.c)
-include $(addprefix .deps/,$(addsuffix .d,$(SOURCES)))

