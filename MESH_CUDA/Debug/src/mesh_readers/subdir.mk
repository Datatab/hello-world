################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/mesh_readers/MeshReaderUnv.cpp 

OBJS += \
./src/mesh_readers/MeshReaderUnv.o 

CPP_DEPS += \
./src/mesh_readers/MeshReaderUnv.d 


# Each subdirectory must supply rules for building sources it contributes
src/mesh_readers/%.o: ../src/mesh_readers/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-10.0/bin/nvcc -G -g -O0 -std=c++11 -gencode arch=compute_30,code=sm_30  -odir "src/mesh_readers" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-10.0/bin/nvcc -G -g -O0 -std=c++11 --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


