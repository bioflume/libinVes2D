#location of mex compiler and mex libraries (you can override using env vars, i.e. env MEX=/my/mex make)
MATLABDIR ?= /opt/matlab/
MEX       ?= ${MATLABDIR}/bin/mex
MEXLIBDIR ?= ${MATLABDIR}/bin/bin/glnxa64/

FORTRAN = ifort
FFLAGS  = -c -fpic -openmp -O3
MEXLIB  = -lifcore -lirc -lsvml -limf -lgfortran -liomp5 

OBJS=\
lfmm2drouts.o\
cfmm2dpart.o\
rfmm2dpart.o\
d2mtreeplot.o\
d2tstrcr_omp.o\
l2dterms.o\
laprouts2d.o\
prini.o\
stokesSLP.o

all: cfmm2dpart d2mtreeplot d2tstrcr_omp l2dterms laprouts2d lfmm2drouts prini rfmm2dpart stokesSLP mex

mex: stokesSLPfmm.F $(OBJS) 
	$(MEX) -O stokesSLPfmm.F $(OBJS) -L $(MEXLIBDIR) $(MEXLIB)

cfmm2dpart: cfmm2dpart.f
	$(FORTRAN) $(FFLAGS) cfmm2dpart.f

d2mtreeplot: d2mtreeplot.f
	$(FORTRAN) $(FFLAGS) d2mtreeplot.f

d2tstrcr_omp: d2tstrcr_omp.f
	$(FORTRAN) $(FFLAGS) d2tstrcr_omp.f

l2dterms: l2dterms.f
	$(FORTRAN) $(FFLAGS) l2dterms.f

laprouts2d: laprouts2d.f
	$(FORTRAN) $(FFLAGS) laprouts2d.f

lfmm2drouts: lfmm2drouts.f
	$(FORTRAN) $(FFLAGS) lfmm2drouts.f

prini: prini.f
	$(FORTRAN) $(FFLAGS) prini.f

rfmm2dpart: rfmm2dpart.f
	$(FORTRAN) $(FFLAGS) rfmm2dpart.f

stokesSLP: stokesSLP.f
	$(FORTRAN) $(FFLAGS) stokesSLP.f

clean: 
	rm *.o
	cp *.mex* ../src
	rm *.mex*
