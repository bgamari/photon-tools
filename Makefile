PROGS=
DESTDIR?=/usr/local

CXXFLAGS += -std=c++11

all : programs

include picoharp/Makefile

programs : ${PROGS}

.PHONY : install
install : ${PROGS}
	cp ${PROGS} ${DESTDIR}/bin

.PHONY : clean
clean : 
	rm -f ${CLEAN} ${PROGS} *.o *.pyc

readme.pdf : readme.mkd
	pandoc $< -o $@
