################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Sources/log.c \
../Sources/main.c \
../Sources/mcg.c \
../Sources/toa.c \
../Sources/uart.c 

OBJS += \
./Sources/log.o \
./Sources/main.o \
./Sources/mcg.o \
./Sources/toa.o \
./Sources/uart.o 

C_DEPS += \
./Sources/log.d \
./Sources/main.d \
./Sources/mcg.d \
./Sources/toa.d \
./Sources/uart.d 


# Each subdirectory must supply rules for building sources it contributes
Sources/%.o: ../Sources/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross ARM C Compiler'
	arm-none-eabi-gcc -mcpu=cortex-m0plus -mthumb -O0 -fmessage-length=0 -fsigned-char -ffunction-sections -fdata-sections  -g3 -I"../Sources" -I"../Includes" -std=c99 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


