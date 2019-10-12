################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/solvers/Fvm_tvd_implicit.cpp \
../src/solvers/MatrixSolver.cpp \
../src/solvers/Method.cpp \
../src/solvers/SolveEilerEq.cpp \
../src/solvers/Solver.cpp 

CU_SRCS += \
../src/solvers/SolveHeatEq.cu 

CU_DEPS += \
./src/solvers/SolveHeatEq.d 

OBJS += \
./src/solvers/Fvm_tvd_implicit.o \
./src/solvers/MatrixSolver.o \
./src/solvers/Method.o \
./src/solvers/SolveEilerEq.o \
./src/solvers/SolveHeatEq.o \
./src/solvers/Solver.o 

CPP_DEPS += \
./src/solvers/Fvm_tvd_implicit.d \
./src/solvers/MatrixSolver.d \
./src/solvers/Method.d \
./src/solvers/SolveEilerEq.d \
./src/solvers/Solver.d 


# Each subdirectory must supply rules for building sources it contributes
src/solvers/%.o: ../src/solvers/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-10.0/bin/nvcc -G -g -O0 -std=c++11 -gencode arch=compute_30,code=sm_30  -odir "src/solvers" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-10.0/bin/nvcc -G -g -O0 -std=c++11 --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/solvers/%.o: ../src/solvers/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-10.0/bin/nvcc -G -g -O0 -std=c++11 -gencode arch=compute_30,code=sm_30  -odir "src/solvers" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-10.0/bin/nvcc -G -g -O0 -std=c++11 --compile --relocatable-device-code=false -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_30  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


