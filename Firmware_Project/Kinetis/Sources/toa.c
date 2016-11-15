#include "toa.h"
#include <math.h>

//converts integers to ascii, returns length of integer string if it worked, -1 for bad input
int itoa(int sup, char * buff, int radix)
{
    if( !buff | radix<=0)
    {
        return -1;
    }
    if(sup==0)
    {
        *buff='0';
        return 1;
    }
    unsigned int val,v;
    int sign=0;
    char temp[32];
    char *tmp=temp;
    if(sup<0&&radix==10)
    {
        sign=1;
        val=-sup;
    }
    else
    {
        val=sup;
    }
    //the string will be in backwards order.
    while( val || tmp==temp)
    {
        v=val%radix;
        val=val/radix;
        
        if(v<=10)
        {
            //should place it in the first location then increment
            *tmp++=v+'0';
        }
        else
        {
            //use capital chars after 9
            *tmp++=v+'A'-10;
        }
    }
    int i=tmp-temp;
    //putting in the sign
    if(sign)
    {
        *buff++ = '-';
        i++;
    }
    //string reverse to the rest of buffer
    while(tmp>temp)
    {
        *buff++ = *--tmp;//predecrement tmp because we have one extra from the while loop above
    }
    return i;
}

//returns length of string if it worked, otherwise -1 for bad input
int ftoa(float sup, char * buff, int after)
{
    if(!buff | after<0)
    {
        return -1;
    }
    int a=(int)sup;
    float b=sup-(float)a;
    
    int i=itoa(a,buff,10);
    if(after>0)
    {
        buff[i]='.';

        b=b*pow((float)10,(float)after); //need math.h
        i+=itoa(b,buff+i+1,10);
    }

    return i;
}


#ifdef UNITTESTS
int main()
{
    int r;
    int a=-820;
    char b[32];
    r=itoa(a,b,10);
    b[r]='\0';
    printf("-820 in base 10 is %s\n",b);
    char d[32];
    r=itoa(a,d,16);
    d[r]='\0';
    printf("-820 in base 16 is %s\n",d);
    char e[32];
    r=itoa(a,e,2);
    e[r]='\0';
    printf("-820 in base 2 is %s\n",e);

    float c=734.524;
    char f[32];
    r=ftoa(c,f,3);
    f[r]='\0';
    printf("734.524 in string form is %s\n",f);
    return 0;
}

#endif
