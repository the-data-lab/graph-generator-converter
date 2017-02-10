exe		= graph500_csr
cxx		= "$(shell which g++)"
cc		= "$(shell which gcc)"

objs	= $(patsubst %.c,%.o,$(wildcard ./generator/*.c)) \
		  $(patsubst %.cpp,%.o,$(wildcard *.cpp))

deps	= $(wildcard ./generator/*.h) \
		  Makefile

cflags	= -std=c99 -Wall -Drestrict=__restrict__ -O3 -DNDEBUG -ffast-math -DGRAPH_GENERATOR_OMP # -g -pg
ldflags	= -O3 -I./generator/ -I. -Wall -fopenmp

ifeq ($(dump), 1)
	ldflags += -DDUMP_RES
endif

%.o:%.c $(deps)
	${cc} -c $< -o $@ $(cflags) $(ldflags) -fopenmp

%.o:%.cpp $(deps)
	${cxx} -c $< -o $@ $(ldflags)
	echo $(objs)

$(exe):$(objs)
	${cxx} $(objs) -o $(exe) $(ldflags) -fopenmp

clean:
	rm -rf $(objs) $(exe) 
