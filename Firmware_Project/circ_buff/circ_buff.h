#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

typedef struct Circbuf {
    void * buffer; 
    void * head;
    void * tail;
    size_t size;
    size_t item_size;
    size_t num_items;
}Circbuf_t;

int Circbuff_add (Circbuf_t * buff, void * data);

int Circbuff_rem (Circbuf_t * buff, void * data);

int Circbuff_init (Circbuf_t * buff, size_t a_size, size_t a_item_size);

int Circbuff_delete(Circbuf_t * buff);
