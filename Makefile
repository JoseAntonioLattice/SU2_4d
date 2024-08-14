FC = gfortran

TARGET = SU2_4d.exe

BIN = bin
SRC = src

FFLAGS = -Wall -Wextra -fcheck=all -O0 -std=f2008

SOURCE = data_types.f90 parameters.f90 arrays.f90 dynamics.f90 main.f90
OBJECT = $(patsubst %, $(BIN)/%, $(SOURCE:.f90=.o) )

$(BIN)/$(TARGET): $(OBJECT)
	$(FC) -o $@ $^

$(BIN)/%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -I$(BIN) -J$(BIN) -c $< -o $@

.PHONY: clean run

clean:
	rm -f $(BIN)/*

run:
	@echo input/input_parameters.par | time $(BIN)/$(TARGET)


