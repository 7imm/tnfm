FC=gfortran
FLAGS=-Wall -fbounds-check -O3 -fbackslash

tfm_tools.o: tfm_tools.f95 tfm_essentials.o 
	$(FC) $(FLAGS) -c -I/usr/include tfm_tools.f95

tfm_essentials.o: tfm_essentials.f95 tfm_constants.o
	$(FC) $(FLAGS) -c tfm_essentials.f95

tfm_density.o: tfm_density.f95 tfm_essentials.o tfm_constants.o
	$(FC) $(FLAGS) -c tfm_density.f95

tfm_temperature.o: tfm_temperature.f95 tfm_essentials.o
	$(FC) $(FLAGS) -c tfm_temperature.f95

tfm_grain.o: tfm_grain.f95 tfm_essentials.o tfm_constants.o
	$(FC) $(FLAGS) -c tfm_grain.f95

tfm_liquid.o: tfm_liquid.f95 tfm_essentials.o tfm_constants.o
	$(FC) $(FLAGS) -c tfm_liquid.f95

tfm_num.o: tfm_num.f95 tfm_essentials.o tfm_constants.o tfm_liquid.o tfm_temperature.o tfm_density.o tfm_grain.o
	$(FC) $(FLAGS) -c tfm_num.f95

tfm_preprocessing.o: tfm_preprocessing.f95 tfm_essentials.o tfm_constants.o
	$(FC) $(FLAGS) -c tfm_preprocessing.f95

tfm_constants.o: tfm_constants.f95
	$(FC) $(FLAGS) -c tfm_constants.f95

tfm: tfm.f95 tfm_essentials.o tfm_num.o tfm_constants.o tfm_liquid.o tfm_temperature.o tfm_density.o tfm_tools.o
	$(FC) $(FLAGS) -o tfm tfm.f95 tfm_essentials.o tfm_num.o tfm_constants.o tfm_liquid.o tfm_temperature.o tfm_grain.o tfm_density.o tfm_tools.o -L/usr/lib -lnetcdff -I/usr/include

clean:
	rm *.o
	rm *.mod
