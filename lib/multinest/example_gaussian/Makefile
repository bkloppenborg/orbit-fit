NESTLIBDIR = ..
LIBS = -L$(NESTLIBDIR) -lmultinest $(LAPACKLIB)
 
OBJFILES = params.o like.o nestwrap.o main.o

all: gaussian 

%.o: %.f90
	$(FC) $(FFLAGS) -I$(NESTLIBDIR) -c $*.f90
 
gaussian: $(OBJFILES)
	$(FC) -o ../gaussian $(OBJFILES) \
	$(FFLAGS) $(LIBS)

clean:
	rm -f *.o *.mod ../gaussian
