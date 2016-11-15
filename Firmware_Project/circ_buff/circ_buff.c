#include "circ_buff.h"


//only going to assume 3 types of data.

//kind of have to assume you are going to give me proper data
//no way to check the data pointed to by the void pointer is the right type
int Circbuff_add (Circbuf_t * buff, void * data)
{
    if( !(buff)|!(data))
    {
        return -1; //bad input
    }
    if(buff->size==0|!(buff->item_size==1|buff->item_size==2|buff->item_size==4))
    {
        return -2;//buffer isn't properly initialized
    }
    if(buff->num_items+1>buff->size)
    {
        return -3;
    }
   //move data in then increment to next spot


    if(buff->item_size==1)
    {
        *(uint8_t *)buff->head=*(uint8_t *)data;
        if(((uint8_t *)buff->head)+1-(uint8_t *)buff->buffer>buff->size)
        {
            buff->head=buff->buffer;
            //(uint8_t *)buff->head=(uint8_t *)buff->buffer;
        }
        else
        {
            (uint8_t *)buff->head++;
        }
    }
    else if (buff->item_size==2)
    {
        *(uint16_t *)buff->head=*(uint16_t *)data;
        if(((uint16_t *)buff->head)+1-(uint16_t *)buff->buffer > buff->size)
        {
            buff->head=buff->buffer;
            //(uint16_t *)buff->head=(uint16_t *)buff->buffer;
        }
        else
        {
            (uint16_t *)buff->head++;
        }
    }
    else if (buff->item_size==4)
    {
        *(uint32_t *)buff->head=*(uint32_t *)data;
        if(((uint32_t *)buff->head)+1-(uint32_t *)buff->buffer >buff->size)
        {
            buff->head=buff->buffer;
            //(uint32_t *)buff->head=(uint32_t *)buff->buffer;
        }
        else
        {
            (uint32_t *)buff->head++;
        }
    }
    buff->num_items++;

    return 0;
}

int Circbuff_rem (Circbuf_t * buff, void * data)
{
    if( !(buff)|!(data))
    {
        return -1;// bad input
    }
    if(!(buff->tail)|!(buff->head)|!(buff->buffer)|!(buff->item_size==1|buff->item_size==2|buff->item_size==4)|buff->size==0|buff->num_items==0)
    {
        return -2;//buffer is either empty or not inited right
    }
    if(buff->item_size==1)
    {
        *(uint8_t *)data=*(uint8_t *)buff->tail;
        if(((uint8_t *)buff->tail)+1-(uint8_t *)buff->buffer > buff->size)
        {
            buff->tail=buff->buffer;
            //(uint8_t *)buff->tail=(uint8_t *)buff->buffer;
        }
        else
        {
            (uint8_t *)buff->tail++;
        }
    }
    
    else if (buff->item_size==2)
    {
        *(uint16_t *)data=*(uint16_t *)buff->tail;
        if(((uint16_t *)buff->tail)+1-(uint16_t *)buff->buffer > buff->size)
        {
            buff->tail=buff->buffer;
            //(uint16_t *)buff->tail=(uint16_t *)buff->buffer;
        }
        else
        {
            (uint16_t *)buff->tail++;
        }
    }

    else if (buff->item_size==4)
    {
        *(uint32_t *)data=*(uint32_t *)buff->tail;
        if(((uint32_t *)buff->tail)+1-(uint32_t *)buff->buffer > buff->size)
        {
            buff->tail=buff->buffer;
            //(uint32_t *)buff->tail=(uint32_t *)buff->buffer;
        }
        else
        {
            (uint32_t *)buff->tail++;
        }
    }

    buff->num_items--;
    return 0;
}


int Circbuff_delete(Circbuf_t * buff)
{
    if(buff!=NULL)
    {
        free(buff->buffer);
        return 0;
    }
    return -1;
}

//init all the pointers to NULL and 
int Circbuff_init (Circbuf_t * buff, size_t a_size, size_t a_item_size)
{
    if(!(buff)|a_size==0|!(a_item_size==1|a_item_size==2|a_item_size==4))
    {
        return -1;
    }
    buff->size=a_size;
    buff->item_size=a_item_size;
    buff->num_items=0;
    buff->buffer=malloc(a_item_size*a_size);
    buff->head=buff->buffer;
    buff->tail=buff->buffer;
    
    return 0;
}


#ifdef UNITTESTS
#define ITER 10
#define ITEM_S 1
void UNIT_TEST1()
{
    Circbuf_t a;
    Circbuff_init(&a,(size_t)ITER,(size_t)ITEM_S);
    int i=0,b, d=0;
    int e=5;
    for(i=0;i<ITER;i++)
    {
        b=Circbuff_add(&a, &e);
        if(b!=0)
        {
            d++;
            printf("returned from Circbuf_add with error %d\n",b);
        }
        if(a.num_items!=i+1)
        {
            d++;
        }
    }
    b=Circbuff_add(&a, &e);
    if(b==0)
    {
        d++;
        printf("added more elements than allowed should return -3, actual return=%d\n",b);
    }
    int * c=&e;
    
    for(i=0;i<ITER;i++)
    {
        b=Circbuff_rem(&a,c);
        if(b!=0)
        {
            d++;
            printf("returned error from Circbuf_rem with error %d, return val was %d\n",b,*c);
        }
        if(a.num_items!=(ITER-1-i))
        {
            d++;
        }
    }
    b=Circbuff_rem(&a,c);
    if(b==0)
    {
        printf("removed an element when there weren't any, should return -2, did return %d\n",b);
        d++;
    }
    if(d==0)
    {
        printf("Unit test 1, fill than empty worked perfectly\n");
    }
    else
    {
        printf("Unit test 1, fill than empty had %d failures\n", d);
    }
    Circbuff_delete(&a);
    return;
}

void UNIT_TEST2()
{
    Circbuf_t a;
    Circbuff_init(&a,ITER,ITEM_S);
    int i=0, b, d=0;
    int * c;
    int e=5;
    c=&e;
    for(i=0;i<ITER/2;i++)
    {
        b=Circbuff_add(&a,&e);
        b=Circbuff_rem(&a,c);
    }
    for(i=0;i<ITER;i++)
    {
        b=Circbuff_add(&a,&e);
        if(a.item_size==1)
        {
            if(((uint8_t*)a.head-(uint8_t*)a.buffer) !=(i+1+ITER/2)%(ITER+1))
            {
                d++;
            }
                 //printf("pointer difference is %d, wanted difference is %d\n",((uint8_t*)a.head-(uint8_t*)a.buffer),(i+1+ITER/2)%(ITER+1));
        }
        if(a.item_size==2)
        {
            if(((uint16_t*)a.head-(uint16_t*)a.buffer)!=(i+1+ITER/2)%(ITER+1))
            {
                d++;
            }
        }
        if(a.item_size==4)
        {
            if(((uint32_t*)a.head-(uint32_t*)a.buffer)!=(i+1+ITER/2)%(ITER+1))
            {
                d++;

            }
        }
    }
    for(i=0;i<ITER;i++)
    {
        b=Circbuff_rem(&a,c);
        if(b)
        {
            printf("b is %d\n",b);
        }
        if(a.item_size==1)
        {
            if(((uint8_t*)a.tail-(uint8_t*)a.buffer)!=(i+1+ITER/2)%(ITER+1))
            {
                d++;
            }
        }
        if(a.item_size==2)
        {
            if(((uint16_t*)a.tail-(uint16_t*)a.buffer)!=(i+1+ITER/2)%(ITER+1))
            {
                d++;
            }
        }
        if(a.item_size==4)
        {
            if(((uint32_t*)a.tail-(uint32_t*)a.buffer)!=(i+1+ITER/2)%(ITER+1))
            {
                d++;
            }
        }
    }
    if(d>0)
    {
        printf("failed UNITTEST2 with %d errors\n",d);
        return;
    }
    printf("passed UNITEST2\n");
    Circbuff_delete(&a);
    return;
}


int main()
{
   UNIT_TEST1();
   UNIT_TEST2();
}

#endif
