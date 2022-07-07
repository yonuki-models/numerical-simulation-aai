# module load odyssey mpi-fftw
MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
SRC_DIR = $(MAKEFILE_DIR)../src/
BUILD_DIR = $(MAKEFILE_DIR)../build/
RUN_DIR = $(MAKEFILE_DIR)../run/

FC = mpifrtpx
OPTS = -Kfast -Kparallel -Kopenmp -SSL2MPI -SSL2BLAMP
LIB_FLOQUET = -SSL2

FLOQUET_MODS = floquet_eigen_problem.o
FLOQUET_MODS += floquet_solver.o
FLOQUET_MODS += floquet_governing_equation.o
FLOQUET_MODS += floquet_control.o
FLOQUET_MODS += parallel.o
FLOQUET_MODS += misk.o

FLOQUET_SWEEP_TARGET = floquet_sweep.exe
FLOQUET_SWEEP_MAIN_OBJ = floquet_sweep.o

FLOQUET_SIMULATION_TARGET = floquet_simulation.exe
FLOQUET_SIMULATION_MAIN_OBJ = floquet_simulation.o

all: $(FLOQUET_SWEEP_TARGET)
	cp $(BUILD_DIR)*.mod $(SRC_DIR)
	cp $(FLOQUET_SWEEP_TARGET) $(RUN_DIR)

$(FLOQUET_SWEEP_TARGET): $(FLOQUET_SWEEP_MAIN_OBJ) $(FLOQUET_MODS)
	$(FC) $(OPTS) -o $@ $^ $(LIB_FLOQUET)

simulation: $(FLOQUET_SIMULATION_TARGET)
	cp $(BUILD_DIR)*.mod $(SRC_DIR)
	cp $(FLOQUET_SIMULATION_TARGET) $(RUN_DIR)

$(FLOQUET_SIMULATION_TARGET): $(FLOQUET_SIMULATION_MAIN_OBJ) $(FLOQUET_MODS)
	$(FC) $(OPTS) -o $@ $^ $(LIB_FLOQUET)

# make module file
%.o : $(SRC_DIR)%.f90
	$(FC) $(OPTS) -o $@ -c $< $(LIB)

# dependency of object files
floquet_sweep.o: parallel.o floquet_control.o floquet_solver.o

floquet.o: parallel.o floquet_control.o floquet_solver.o

floquet_solver.o: floquet_eigen_problem.o floquet_governing_equation.o misk.o

floquet_simulation.o: floquet_control.o floquet_governing_equation.o

parallel.o: misk.o

# clean
clean :
	rm -f $(SRC_DIR)*.mod
	rm -f $(BUILD_DIR)*.exe
	rm -f $(BUILD_DIR)*.mod
	rm -f $(BUILD_DIR)*.o
	rm -f $(RUN_DIR)*.exe
