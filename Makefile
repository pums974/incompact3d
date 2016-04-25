include Make.Gandalf

.SUFFIXES:
.SUFFIXES:.o .f90

OBJS=\
body.o \
convdiff.o \
convection.o \
derivevitesse.o \
derivetemperature.o \
fft.o \
filtre.o \
forces.o \
geom_complex.o \
incompact3d.o \
lorenz.o \
module_param.o \
navier.o \
parametre.o \
poisson.o \
scalar.o \
schemas.o \
slfft2d.o \
slfft2d_shift.o \
slfft3d_shift.o \
stats.o \
tools.o \
waves.o

SRC=*.f90

all : $(EXENAME)

$(EXENAME):     $(OBJS)
	$(FF) -o $@ $(LDFLAGS) $(OBJS)

.f90.o:
	@echo
	@echo "-----------Compilation ($<)-----------"
	@echo
	$(FF) $(FFLAGS) -c $<

depend : $(SRC)
	makedepf90 $(SRC) > dependance.dep

clean:
	@echo
	@echo "-----------Nettoyage -----------"
	@echo
	rm *.o *.mod *.*~ 


include dependance.dep
