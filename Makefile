include pt2/Makefile

.PHONY : install
install : ${PROGS}
	cp ${PROGS} ${DESTDIR}/bin

.PHONY : clean
clean : 
	rm -f ${CLEAN} ${PROGS} *.o *.pyc


