include=-IlibGTF -IlibGTF/htslib
libs=libGTF/libGTF.a libGTF/htslib/libhts.a 

.PHONY: all htslib libGTF

all: countRepeats

htslib:
	$(MAKE) -C libGTF/htslib

libGTF:
	$(MAKE) -C libGTF

countRepeats: htslib libGTF
	gcc -g $(include) -o $@ cntHash.c EM.c countRepeats.c $(libs) -lz -lpthread -lpcre -lgsl -lm -lblas

clean:
	rm -f countRepeats
