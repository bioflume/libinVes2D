
## Add any compiler flags you like to this
CXXFLAGS := -Wall -Wextra -Wno-unused -Woverloaded-virtual -pipe -fPIC -ftree-vectorize -funroll-loops -ffast-math -shared -fno-common -fpermissive -O3 -g -DGL_GLEXT_PROTOTYPES

## Set up include paths and libraries used
#CPPFLAGS := -Iacm -Iacm/eigen3 -Iacm/vcglib -I/home/shravan/pkgs/glew/include
include iagm/iagm/make.inc
CPPFLAGS := -Iiagm -I${EIGEN_INC} -Iiagm/iagm/inc -Iiagm/iagm/src
LIBS 	 := -Liagm/iagm/lib -liagm

MEXSUFFIX  = mexglx
MATLABHOME = /Applications/MATLAB_R2015b.app
MEX        = /Applications/MATLAB_R2015b.app/bin/mex
CXX        = g++-5
CFLAGS     = -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -DMX_COMPAT_32 \
               -DMATLAB_MEX_FILE
#CFLAGS     = -ansi -D_GNU_SOURCE -fexceptions -fPIC -fno-omit-frame-pointer -DMATLAB_MEX_FILE

MLIBS      = -lm
MATINCLUDE   = -I${MATLABHOME}/extern/include -I${MATLABHOME}/simulink/include
MEXFLAGS  = -cxx CC='${CXX}' CXX='${CXX}' LD='${CXX}' 

#MEX_LIB = -Wl,-rpath-link,/opt/matlab/bin/glnx86 -L/opt/matlab/bin/glnx86 -lmx -lmex -lmat -lm
MEX_LIB = -L${MATLABHOME}/bin/glnxa64 -lmx -lmex -lmat -lm

## Places for make to search for source files
VPATH := ./iagm 
#./iagm/iagm ./iagm/iagm/src

## Everything to build
TARGETS := getCollision.${MEXSUFFIX}

## Normal sources for the code library
SOURCES := getCollision.cpp $(wildcard ./iagm/*.cpp)

OBJDIR := build/

## Objects for the code library
OBJS    := $(addprefix $(OBJDIR), $(addsuffix .o, $(basename $(notdir $(SOURCES)))))

## Build everything
all: $(TARGETS)

## Library of code
getCollision.${MEXSUFFIX}: $(OBJS)
	${MEX} ${MEXFLAGS} ${MEX_LIB} -output getCollision $(LIBS) $^

$(OBJDIR)%.o: %.cpp
	$(shell mkdir -p $(OBJDIR))
	$(CXX) ${CFLAGS} ${MATINCLUDE} $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

depend:
	$(CXX) ${MATINCLUDE} $(CPPFLAGS) -MM $(SOURCES) > .depend

ifeq (.depend, $(wildcard .depend))
  include .depend
endif

clean: 
	-rm -f $(TARGETS) $(OBJDIR)*.o .depend

distclean: clean

rebuild: distclean depend all

