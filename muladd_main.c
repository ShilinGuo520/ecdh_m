#include "stdio.h"
#include "stdlib.h"

int muladd(unsigned int a,unsigned int b)
{
    unsigned int al,ah,bl,bh,cma,cmb,cm;
    unsigned long c ,re;
    al = a & 0x0000ffff;
    ah = (a >> 16) & 0xffff;
    bl = b & 0x0000ffff;
    bh = (b >> 16) & 0xffff;
    cma = ah * bl;
    cmb = al * bh;
    cm = cma + cmb;    
    if((cm < cma) && (cm < cmb))
    	c = (((unsigned long)(ah * bh)) << 32) + (unsigned long)(al * bl) + ((unsigned long)cm << 16) + 0x1000000000000;
    else 
    	c = (((unsigned long)(ah * bh)) << 32) + (unsigned long)(al * bl) + ((unsigned long)cm << 16);
    re = (unsigned long) a * b;
    if(re != c) {
	printf("re:\n");
        printf("%x\n",re);
	printf("%x\n",re>>32);

  	printf("c:\n");
	printf("%x\n",c);
	printf("%x\n",c>>32);
	return 1;
    }
    else
	return 0;
}


int main()
{
    int i = 0xffffffff;
    while(i--) {
	if(muladd(0xffffffff - i, 0xffffeeee)) {
	    printf("error:%x\n",i);
	    return 1;
  	}
    }
    printf("muladd is ok:\n");
    return 0;
}

