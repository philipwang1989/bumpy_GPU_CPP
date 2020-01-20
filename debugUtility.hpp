#ifndef _BUMPY_DBU_H_
#define _BUMPY_DBU_H_

#include <iostream>

void printArray(int Nc, double * x)
{
    // printf("m, I\n");
    for (int j=0;j<Nc;j++)
    {
        printf("%1.6f ", x[j]);
    }
    printf("\n");
    // return;
}

#endif
