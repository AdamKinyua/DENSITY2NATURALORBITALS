mqcroot = "$(HOME)/.mqc"
FC = gfortran
ifeq ($(FC),gfortran)
	FCFLAGS = -std=f2008 -fdefault-real-8 -fdefault-integer-8 -fopenmp
	MQCOBJS = -I$(mqcroot)/GNU/mod
	LIBS = -llapack -lblas $(mqcroot)/GNU/lib/libmqc.a
else ifeq ($(FC),pgfortran)
	FCFLAGS = -Mallocatable=03 -r8 -i8 -mp
	MQCOBJS = -module $(mqcroot)/PGI/mod
	LIBS = -llapack -lblas $(mqcroot)/PGI/lib/libmqc.a
endif

all: compile

compile: density2NatOrbs.f03
	$(FC) $(FCFLAGS) -g $(MQCOBJS) -o density2NatOrbs.exe density2NatOrbs.f03 $(LIBS) 

clean:
	rm -f density2NatOrbs.exe density2NatOrbs.o

