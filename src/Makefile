# Find all source files, create a list of corresponding object files
SRCS=$(wildcard *.F90)
OBJS=$(patsubst %.F90,%.o,$(SRCS))

# Ditto for mods (They will be in both lists)
MODS=$(wildcard mod*.F90)
MOD_OBJS=$(patsubst %.F90,%.o,$(MODS))

# Compiler/Linker settings
FC = ifort
FLFLAGS = -g -L/home/wei.wang/CGNS-3.3.1/build/install/lib -lcgns
FCFLAGS = -g -c -I/home/wei.wang/CGNS-3.3.1/build/install/include
PROGRAM = gensurf
PRG_OBJ = $(PROGRAM).o

# make without parameters will make first target found.
default : $(PROGRAM)

# Compiler steps for all objects
$(OBJS) : %.o : %.F90
	$(FC) $(FCFLAGS) -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^

# If something doesn't work right, have a 'make debug' to 
# show what each variable contains.
debug:
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "MODS = $(MODS)"
	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "PRG_OBJ = $(PRG_OBJ)"

clean:
	rm -rf $(OBJS) $(PROGRAM) $(patsubst %.o,%.mod,$(MOD_OBJS))

.PHONY: debug default clean

# Dependencies

# Main program depends on all modules
$(PRG_OBJ) : $(MOD_OBJS)

# Blocks and allocations depends on certain mods
mod_func.o : mod_constant.o
