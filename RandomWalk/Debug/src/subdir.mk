################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Copy\ of\ Main.cpp \
../src/Main.cpp 

OBJS += \
./src/Copy\ of\ Main.o \
./src/Main.o 

CPP_DEPS += \
./src/Copy\ of\ Main.d \
./src/Main.d 


# Each subdirectory must supply rules for building sources it contributes
src/Copy\ of\ Main.o: ../src/Copy\ of\ Main.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/ramon/MEGAsync/prodgraph/code/RandomWalk/lib/eigen -O3 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fopenmp -MMD -MP -MF"src/Copy of Main.d" -MT"src/Copy\ of\ Main.d" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/ramon/MEGAsync/prodgraph/code/RandomWalk/lib/eigen -O3 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


