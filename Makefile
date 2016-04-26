include Make.Bardwell

.SUFFIXES:
.SUFFIXES:.o .f90

OBJS=\
tools.o \
module_param.o \
derivevitesse.o \
scalar.o \
filtre.o \
geom_complex.o \
convdiff.o \
derivetemperature.o \
fft.o \
forces.o \
lorenz.o \
parametre.o \
slfft2d.o \
slfft2d_shift.o \
slfft3d_shift.o \
poisson.o \
schemas.o \
stats.o \
waves.o \
navier.o \
convection.o \
body.o \
incompact3d.o

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
	rm *.o *.mod *.*~  *__genmod.f90


include dependance.dep
