all: rotate

rotate: rotate.cc
	g++ -I /usr/include/netpbm -o rotate rotate.cc -Wall -lnetpbm

clean:
	rm -f *~

clean-all: clean
	rm -f rotate
