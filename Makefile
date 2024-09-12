FC = gfortran

TARGET = SU2_4d.exe

BIN = bin
SRC = src

#FFLAGS = -Wall -Wextra -fcheck=all -O0 -std=f2008
FFLAGS = -O3

SOURCE = statistics.f90 precision.f90 parameters.f90 number2string_mod.f90 check_files_directories_mod.f90 create_files.f90 pbc.f90 datatypes.f90 starts.f90 arrays.f90 wilson_loops.f90 smooth_configurations.f90 observables.f90 local_update_algorithms.f90 dynamics.f90 main.f90


OBJECT = $(patsubst %, $(BIN)/%, $(SOURCE:.f90=.o) )

compile: $(OBJECT)
	$(FC) -o $(BIN)/$(TARGET) $^ -llapack

$(BIN)/%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -I$(BIN) -J$(BIN) -c $< -o $@ 

.PHONY: clean run install load_module

install:
	chmod +x install.sh
	./install.sh

load_module:
	module load lamod/gcc/12.2

clean:
	rm -f $(BIN)/*

run:
	{ echo 'input/input_parameters.par'; echo 'input/beta.in'; } | time $(BIN)/$(TARGET)


