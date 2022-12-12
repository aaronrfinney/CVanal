# Makefile for CVanal (ARF 07/17)

EXE = ./CVanal

OBJ1 = variables.o channels.o error_close.o util.o

OBJ2 = readinp.o toxyz.o readcon.o readfld.o  \
	cvatdef.o readhist.o images.o coordinate.o \
	dstrb_smooth.o

OBJ3 = CVanal.o

#========= For testing and debugging  ================
deb: 
	$(MAKE) CF="gfortran -Wall -pedantic -fbacktrace -fbounds-check -c" \
	LD="gfortran -Wall -pedantic -fbacktrace -fbounds-check" \
	LDFLAGS="-L/usr/local/lib " \
	EXE=$(EXE).deb LIBS=""  makit

#========= For production runs  ======================
prod: 
	$(MAKE) CF="gfortran -c -O3" \
	LD="gfortran -O3" \
	LDFLAGS="-L/usr/local/lib " \
	EXE=$(EXE) LIBS=""  makit

#=====================================================
makit:  $(OBJ1) $(OBJ2) $(OBJ3)
	$(LD) -o $(EXE) $(LDFLAGS) $(OBJ1) $(OBJ2) $(OBJ3) $(LIBS) 

clean:
	\rm *.o *.T *.tmp.f *.mod $(EXE)

# declare dependencies
%.o:%.f90
	$(CF) $*.f90
.f.o:
	$(CF) $*.f

#.F.o:
#       $(GPP) $(GPPFLAGS) $*.F $*.tmp.f
#       $(CF) $(CFLAG) $*.tmp.f
#       mv $*.tmp.o $*.o


