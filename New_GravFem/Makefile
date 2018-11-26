fem_objects = ElmIntCode.o Grav_main.o
fem_modules = global_grav.o

f90comp = ifort -qopenmp

FFLAGS = -O2
FFLAGS = -g
#FFLAGS = -nostandard-realloc-lhs
FFLAGS = -check bounds
FFLAGS = -debug
#FFLAGS = -threads specified

%.o : %.f90
	$(f90comp) $(FFLAGS) -c  $<


FemGrav : $(fem_objects) $(fem_modules) 
	$(f90comp) $(FFLAGS) -o FemGrav $(fem_modules) $(fem_objects) 

$(fem_objects) : $(fem_modules)

clean:
	rm $(fem_objects) $(fem_modules) 

