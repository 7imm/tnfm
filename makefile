FC=gfortran
FLAGS=-Wall -fbounds-check -O3 -fbackslash

test: test.f95 settings.o medley.o merra2.o tfm_tools.o tfm_density.o tfm_temperature.o tfm_liquid.o tfm_num.o tfm_preprocessing.o
	$(FC) $(FLAGS) -o test test.f95 settings.o medley.o merra2.o tfm_tools.o tfm_density.o tfm_temperature.o tfm_liquid.o tfm_num.o tfm_preprocessing.o tfm_constants.o -L/usr/lib -lnetcdff

merra2.o: merra2.f95 settings.o settings.o tfm_tools.o
	$(FC) $(FLAGS) -c -I/usr/include merra2.f95

tfm_tools.o: tfm_tools.f95 settings.o settings.o
	$(FC) $(FLAGS) -c -I/usr/include tfm_tools.f95

settings.o: settings.f95
	$(FC) $(FLAGS) -c settings.f95

medley.o: medley.f95 settings.o merra2.o tfm_density.o tfm_temperature.o tfm_liquid.o tfm_num.o tfm_preprocessing.o tfm_constants.o
	$(FC) $(FLAGS) -c medley.f95

tfm_density.o: tfm_density.f95 settings.o tfm_constants.o
	$(FC) $(FLAGS) -c tfm_density.f95

tfm_temperature.o: tfm_temperature.f95 settings.o
	$(FC) $(FLAGS) -c tfm_temperature.f95

tfm_liquid.o: tfm_liquid.f95 settings.o tfm_constants.o
	$(FC) $(FLAGS) -c tfm_liquid.f95

tfm_num.o: tfm_num.f95 settings.o tfm_constants.o tfm_liquid.o tfm_temperature.o tfm_density.o
	$(FC) $(FLAGS) -c tfm_num.f95

tfm_preprocessing.o: tfm_preprocessing.f95 settings.o tfm_constants.o
	$(FC) $(FLAGS) -c tfm_preprocessing.f95

tfm_constants.o: tfm_constants.f95
	$(FC) $(FLAGS) -c tfm_constants.f95

melt_test: melt_test.f95 settings.o medley.o merra2.o tfm_tools.o tfm_density.o tfm_temperature.o tfm_liquid.o tfm_num.o tfm_preprocessing.o
	$(FC) $(FLAGS) -o melt_test melt_test.f95 settings.o medley.o merra2.o tfm_tools.o tfm_density.o tfm_temperature.o tfm_liquid.o tfm_num.o tfm_preprocessing.o tfm_constants.o -L/usr/lib -lnetcdff

retmip: retmip.f95 settings.o tfm_num.o tfm_temperature.o tfm_temperature.o tfm_liquid.o tfm_tools.o
	$(FC) $(FLAGS) -o retmip retmip.f95 settings.o tfm_num.o tfm_temperature.o tfm_density.o tfm_liquid.o tfm_tools.o -L/usr/lib -lnetcdff -I/usr/include

clean:
	rm *.o
	rm *.mod
