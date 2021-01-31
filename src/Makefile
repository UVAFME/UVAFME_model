PROG =	UVAFME_v1_NABoreal

SRCS =	Climate.f90 Constants.f90 csv_file.f90 GenusGroups.f90 Input.f90       \
	IO_Utils.f90 lists.f90 Model.f90 Output.f90 Parameters.f90 Plot.f90        \
	Random.f90 Site.f90 Soil.f90 Species.f90 strings.f90 Tree.f90              \
	Utilities.f90 UVAFME.f90 FileUtils.F90

OBJS =	Climate.o Constants.o csv_file.o GenusGroups.o Input.o IO_Utils.o      \
	lists.o Model.o Output.o Parameters.o Plot.o Random.o Site.o Soil.o        \
	Species.o strings.o Tree.o Utilities.o UVAFME.o FileUtils.o                \

LIBS =

F90 =  ifort

#DBG = -CB -g -traceback -fp-model strict
#OPT = -O2 -axcore-avx2 -xavx -align array64byte -ipo

F90FLAGS = $(DBG) $(OPT)
LDFLAGS =
all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

.PHONY: clean
clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90 .F90 .f95 .f

%.o: %.mod

.f90.o .f95.o .F90.o:
	$(F90) $(F90FLAGS) -c $<

Climate.o: Constants.o Random.o
csv_file.o: csv_file_1d.f90 csv_file_2d.f90
Input.o: Constants.o IO_Utils.o Parameters.o Species.o FileUtils.o
IO_Utils.o: Constants.o FileUtils.o Parameters.o csv_file.o
GenusGroups.o: Constants.o Site.o
lists.o: Species.o
Model.o: Climate.o Constants.o Input.o Parameters.o Random.o Site.o Soil.o     \
	Species.o Tree.o
Output.o: Constants.o IO_Utils.o Tree.o csv_file.o
Parameters.o: Constants.o
Plot.o: Constants.o Species.o Tree.o Soil.o Input.o
Site.o: Constants.o Input.o Plot.o Soil.o Species.o Utilities.o lists.o
Soil.o: csv_file.o Constants.o FileUtils.o Parameters.o
Species.o: Constants.o
Tree.o: Constants.o Random.o Species.o
Utilities.o: Constants.o
UVAFME.o: Constants.o GenusGroups.o Input.o Model.o Output.o Parameters.o      \
	Site.o Species.o
FileUtils.o: Constants.o Parameters.o strings.o
