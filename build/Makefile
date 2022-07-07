# module load odyssey mpi-fftw
MAKEFILE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
SRC_DIR = $(MAKEFILE_DIR)../src/
BUILD_DIR = $(MAKEFILE_DIR)../build/
RUN_DIR = $(MAKEFILE_DIR)../run/

FC = mpifrtpx
OPTS = -Kfast -Kparallel -Kopenmp -SSL2MPI -SSL2BLAMP
LIB_FLOQUET = -SSL2

MODS = analysis.o
MODS += control.o
MODS += dynamics.o
MODS += fftw_interface.o
MODS += geometry.o
MODS += external_field.o
MODS += initial_condition.o
MODS += input_output.o
MODS += integration.o
MODS += log.o
MODS += misk.o
MODS += parallel.o
MODS += parameters.o
MODS += simulation.o
MODS += time_module.o

TARGET = main.exe
MAIN_OBJ = main.o

all: $(TARGET)
	cp $(BUILD_DIR)*.mod $(SRC_DIR)
	cp $(TARGET) $(RUN_DIR)

$(TARGET): $(MAIN_OBJ) $(MODS)
	$(FC) $(OPTS) -o $@ $^ $(LIB) $(LDFLAGS)

# make module file
%.o : $(SRC_DIR)%.f90
	$(FC) $(OPTS) -o $@ -c $< $(LIB)

# dependency of object files
main.o: time_module.o simulation.o

simulation.o: parallel.o input_output.o geometry.o time_module.o log.o  \
 control.o initial_condition.o external_field.o integration.o \
 parameters.o

integration.o: geometry.o time_module.o external_field.o dynamics.o \
 analysis.o

input_output.o: parallel.o time_module.o geometry.o external_field.o  \
 analysis.o dynamics.o

dynamics.o: geometry.o parameters.o analysis.o

external_field.o: parameters.o time_module.o geometry.o

initial_condition.o: geometry.o parallel.o parameters.o external_field.o \
 time_module.o

log.o: input_output.o geometry.o control.o time_module.o parallel.o  \
 parameters.o external_field.o analysis.o

analysis.o: geometry.o external_field.o misk.o

control.o: input_output.o parameters.o geometry.o time_module.o  \
 external_field.o

geometry.o: fftw_interface.o parallel.o misk.o

fftw_interface.o: parallel.o misk.o

parallel.o: misk.o

# clean
clean :
	rm -f $(SRC_DIR)*.mod
	rm -f $(BUILD_DIR)*.exe
	rm -f $(BUILD_DIR)*.mod
	rm -f $(BUILD_DIR)*.o
	rm -f $(RUN_DIR)*.exe
