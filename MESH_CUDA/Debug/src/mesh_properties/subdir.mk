################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/mesh_properties/BndFaceTemperature.cpp \
../src/mesh_properties/CellFluidDynamicsProps.cpp \
../src/mesh_properties/CellTemperature.cpp 

OBJS += \
./src/mesh_properties/BndFaceTemperature.o \
./src/mesh_properties/CellFluidDynamicsProps.o \
./src/mesh_properties/CellTemperature.o 

CPP_DEPS += \
./src/mesh_properties/BndFaceTemperature.d \
./src/mesh_properties/CellFluidDynamicsProps.d \
./src/mesh_properties/CellTemperature.d 


# Each subdirectory must supply rules for building sources it contributes
src/mesh_properties/%.o: ../src/mesh_properties/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-10.0/bin/nvcc -G -g -O0 -std=c++11 -gencode arch=compute_30,code=sm_30  -odir "src/mesh_properties" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-10.0/bin/nvcc -G -g -O0 -std=c++11 --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


