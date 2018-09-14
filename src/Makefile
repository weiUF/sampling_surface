##Makefile 

FC = f95

src = $(wildcard *.F90)
obj = $(src:.F90=.o)

DEBUG = -g
CPPFLAG = -I/root/temp/CGNS-3.3.1/src
LDFLAG = -L/root/temp/CGNS-3.3.1/src/lib
LIB = -lcgns

gensurf: $(obj)
	$(FC) -o $@ $^ $(LDFLAG) $(LIB)

%.o:%.F90
	$(FC) -c $^ -o $@ $(CPPFLAG) $(DEBUG)

.PHONY: clean

clean:
	rm -f $(obj) gensurf

