## Makefile for FrTMalign
## 

PROG =	frtmalign

SRCS =	main.f  main_dp.f  main_fragment.f  main_initial.f  main_tmsc.f

OBJS =	main.o  main_dp.o  main_fragment.o  main_initial.o  main_tmsc.o

## choose appropriate compiler, which is present on your OS

FC = gfortran
FFLAGS = -ffixed-line-length-none 

all: $(PROG)

$(PROG):$(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) 

clean: 

	rm  ${PROG} ${OBJS}

# End of the makefile
