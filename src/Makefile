##Makefile 

FC = f95

src = $(wildcard *.F90)
obj = $(src:.F90=.o)
dep = $(obj:.o=.d)

DEBUG = -g
CPPFLAG = -I/root/CGNS/src
LDFLAG = -L/root/CGNS/src/lib
LIB = -lcgns

ALL:gensurf


%.o:%.F90
	$(FC) -c $< -o $@ $(CPPFLAG) $(DEBUG)

gensurf: $(obj)
	$(FC) -o $@ $^ $(LDFLAG) $(LIB)

.PHONY: ALL

.PHONY: clean

clean:
	rm -f $(obj) gensurf *.mod

