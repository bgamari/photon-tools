include picoharp/Makefile

.PHONY : install
install : ${PROGS}
	cp ${PROGS} ${DESTDIR}/bin

.PHONY : clean
clean : 
	rm -f ${CLEAN} ${PROGS} *.o *.pyc

readme.pdf : readme.mkd
	pandoc $< -o $@
