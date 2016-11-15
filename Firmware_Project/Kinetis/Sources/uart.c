/*
 * uart.c
 *
 *  Created on: Apr 9, 2016
 *      Author: jimmy
 */

#include "uart.h"

#define CORE_CLOCK 24000000
#define OVER_SAMPLE 7

uint8_t rx_data;
int j;

extern void UART0_IRQHandler(){
	rx_data = (uint8_t)UART0->D;
	j = j+1;
}

void uart_init(uint32_t baud_rate){
	SIM_SCGC5 |= SIM_SCGC5_PORTA_MASK;

	// Turn on uart0 clock and select 48MHz (PLL) clock
	SIM_SCGC4 |= SIM_SCGC4_UART0_MASK;
	SIM_SOPT2 |= SIM_SOPT2_UART0SRC(1) | SIM_SOPT2_PLLFLLSEL(1);

	// Select "Alt 2" UART0 functionality for pins A1 and A2
	PORTA_PCR1 = PORT_PCR_MUX(2) | PORT_PCR_ISF_MASK;
	PORTA_PCR2 = PORT_PCR_MUX(2) | PORT_PCR_ISF_MASK;

	// Set up correct baud rate = 9600
	uint16_t baud_divisor = CORE_CLOCK/(baud_rate * (OVER_SAMPLE + 1));
	UART0_BDH = (baud_divisor >> 8) & UART0_BDH_SBR_MASK;
	UART0_BDL = (baud_divisor & UART0_BDL_SBR_MASK);
	UART0_C4 = UART0_C4_OSR(OVER_SAMPLE);

	// Enable transmit, receive, and receive interrupts
	UART0_C2 = UART0_C2_TE_MASK | UART0_C2_RE_MASK | UART0_C2_RIE_MASK;
	NVIC->ICPR[0] |= 1 << 12;
	NVIC->ISER[0] |= (1 << 12);
}


