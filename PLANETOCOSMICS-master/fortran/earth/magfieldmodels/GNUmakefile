OBJECTS = Tsy89c_mod.o  Tsy96_magcos.o  igrf_to_c++.o  \
	  Tsyg_01_scaled_to_B0_magcos.o TS04c.o

SRC = Tsy89c_mod.for  Tsy96_magcos.for  igrf_to_c++.for  \
	  Tsyg_01_scaled_to_B0_magcos.for TS04c.for	  

C = gfortran -O  

lib/libFortranEarthMagField.a : $(SRC)
	@echo "compiling ..."
	@$(C) -c $(SRC)
	@if [ ! -d lib ]; then mkdir lib; fi
	@echo "creating  lib/libFortranEarthMagField.a ..."
	@ar -r lib/libFortranEarthMagField.a *.o
	@rm *.o

clean  :
	@if [ -d lib ];then echo "cleaning ..."; rm -r  lib; else  echo "Nothing to clean"; fi	
