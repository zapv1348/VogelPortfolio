/*
 * log.c
 *
 *  Created on: Apr 9, 2016
 *      Author: jimmy
 */
#ifndef DISABLE_LOGS


#include "log.h"
#include "stdio.h"
#include "toa.h"

void LOG_0(uint8_t * str, size_t length){
	int array_loop;

	for(array_loop = 0; array_loop < length; array_loop++){
		while(!(UART0->S1 & UART_S1_TDRE_MASK) && !(UART0->S1 & UART_S1_TC_MASK));
		UART0->D = str[array_loop];
	}

}

void LOG_1(uint8_t * str, size_t length, size_t param, size_t data_type_size, uint8_t int_len, float floaty){
	int array_loop;
	// Print the text part of the loop
	for(array_loop = 0; array_loop < length; array_loop++){
		while(!(UART0->S1 & UART_S1_TDRE_MASK) && !(UART0->S1 & UART_S1_TC_MASK));
		UART0->D = str[array_loop];
	}

	// Generate the integer part of the string for an 8-bit int
	if(data_type_size == 8){
		uint8_t * new_ascii[int_len + 1];

		// Convert integer to ascii
		UlToStr(&new_ascii[0], (uint32_t) param, int_len, data_type_size);

		// Transmit new ascii value over UART
		for(array_loop = 0; array_loop < int_len; array_loop++){
			while(!(UART0->S1 & UART_S1_TDRE_MASK) && !(UART0->S1 & UART_S1_TC_MASK));
			UART0->D = (uint32_t)new_ascii[0] >> 8*array_loop;
		}
	}

	// Generate the integer part of the string for an 16-bit int
	else if(data_type_size == 16){
		uint16_t * new_ascii[int_len + 1];

		// Convert integer to ascii
		UlToStr(&new_ascii[0], (uint32_t) param, int_len, data_type_size);

		// Transmit new ascii value over UART
		for(array_loop = 0; array_loop < int_len; array_loop++){
			while(!(UART0->S1 & UART_S1_TDRE_MASK) && !(UART0->S1 & UART_S1_TC_MASK));
			if(array_loop < 2)
				UART0->D = (uint32_t)new_ascii[0] >> 16*array_loop;
			else
				UART0->D = (uint32_t)new_ascii[1] >> 16*(array_loop - 2);

		}
	}

	// Generate the integer part of the string for an 32-bit int
	else if(data_type_size == 32){
		uint32_t * new_ascii[int_len + 1];

		// Convert integer to ascii
		UlToStr(&new_ascii[0], (uint32_t) param, int_len, data_type_size);

		// Transmit new ascii value over UART
		for(array_loop = 0; array_loop < int_len; array_loop++){
			while(!(UART0->S1 & UART_S1_TDRE_MASK) && !(UART0->S1 & UART_S1_TC_MASK));
			UART0->D = (uint32_t)new_ascii[array_loop];

		}
	}

	// Generate the floating point part of the string for an float
	else {
		char new_ascii[36];
		float newbie = floaty;

		// Convert float to ascii
		ftoa(newbie, new_ascii, 3);

		// Transmit new ascii value over UART
		for(array_loop = 0; array_loop < int_len; array_loop++){
			while(!(UART0->S1 & UART_S1_TDRE_MASK) && !(UART0->S1 & UART_S1_TC_MASK));
			UART0->D = (uint32_t)new_ascii[array_loop];
		}
	}

}


// Function that converts an integer to an ascii value
void UlToStr(void *s, uint32_t bin, uint8_t n, uint32_t data_size)
{
	if(data_size == 8){
		uint8_t * ptr = s;
		ptr += n;
		*ptr = '\0';

		while (n--)
		{
			*--ptr = (bin % 10) + '0';
			bin /= 10;
		}
		return ptr;
	}
	else if(data_size == 16){
		uint16_t * ptr = s;
		ptr += n;
		*ptr = '\0';

		while (n--)
		{
			*--ptr = (bin % 10) + '0';
			bin /= 10;
		}
		return ptr;
	}
	else if(data_size == 32){
		uint32_t * ptr = s;
		ptr += n;
		*ptr = '\0';

		while (n--)
		{
			*--ptr = (bin % 10) + '0';
			bin /= 10;
		}
		return ptr;
	}
}

#endif
