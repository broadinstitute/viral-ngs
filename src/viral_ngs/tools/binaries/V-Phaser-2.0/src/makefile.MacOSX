PROG=../MacOSX/variant_caller
BAMTOOLSPATH=../bamtools
BOOSTPATH=/usr/local/include/boost
COMPILER=/usr/local/bin/g++-4.9

all:
	$(COMPILER) -fopenmp -O3 -I$(BAMTOOLSPATH)/include/ -I$(BOOSTPATH)/ $(BAMTOOLSPATH)/lib/libjsoncpp.a $(BAMTOOLSPATH)/lib/libbamtools-utils.a $(BAMTOOLSPATH)/lib/libbamtools.a /usr/lib/libz.dylib *.cpp -o $(PROG) -lpthread
	
clean:
	rm -rf $(PROG)