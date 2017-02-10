all: gConvu gConvd

gConvu: gConv.h gConv.cpp
	g++ -g -DHAL_GRID gConv.cpp -fopenmp -o gConvu

gConvd: gConv.h gConv.cpp
	g++ -g gConv.cpp -fopenmp -o gConvd

clean:
	rm gConvd gConvu
