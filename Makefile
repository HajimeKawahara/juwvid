FC=/usr/local/bin/gfortran
FFLAGS=-fPIC
OBJS=nufft_src/nufft1df90.o nufft_src/dirft1d.o nufft_src/dfftpack.o nufft_src/next235.o

all: libnufft.so

libnufft.so: $(OBJS)	
	$(FC) -shared $(FFLAGS) $(OBJS) ionufft.f90 -o libnufft.so


clean: 
	\rm -rf *.o *.mod libnufft.so
