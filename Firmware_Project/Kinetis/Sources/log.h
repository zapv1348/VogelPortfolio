/*
 * log.h
 *
 *  Created on: Apr 9, 2016
 *      Author: jimmy
 */

#ifndef SOURCES_LOG_H_
#define SOURCES_LOG_H_

#include "MKL25Z4.h"
#include "stddef.h"

void LOG_0(uint8_t * str, size_t length);
void LOG_1(uint8_t * str, size_t length, size_t param, size_t data_type_size, uint8_t int_len, float floaty);
void UlToStr(void *s, uint32_t bin, uint8_t n, uint32_t data_size);


#endif /* SOURCES_LOG_H_ */
