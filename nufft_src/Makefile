FC=/usr/local/bin/gfortran
FFLAGS=-fPIC
OBJS=nufft1df90.o dirft1d.o dfftpack.o next235.o

all: libnufft.so sample

sample: sample.o $(OBJS) 
	$(FC) $(FFLAGS) -o sample sample.o $(OBJS) 

libnufft.so: $(OBJS)	
	$(FC) -shared $(FFLAGS) $(OBJS) ionufft.f90 -o libnufft.so 

clean: 
	\rm -rf *.o *.mod
