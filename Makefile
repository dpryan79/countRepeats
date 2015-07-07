include=-I../libGTF -I../htslib
libs=../libGTF/libGTF.a ../htslib/libhts.a 

.PHONY: all htslib libGTF

all: countRepeats

htslib:
	$(MAKE) -C ../htslib

libGTF:
	$(MAKE) -C ../libGTF

countRepeats: htslib libGTF
	gcc -g $(include) -o $@ cntHash.c EM.c countRepeats.c $(libs) -lz -lpthread -lpcre -lgsl -lm -lblas

clean:
	rm -f countRepeats
