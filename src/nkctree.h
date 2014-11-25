#ifndef NKCTREE_H
#define NKCTREE_H
#include "cvector.h"

typedef struct FR_nkctree{
    unsigned int row_count;
    unsigned int set_size;
    unsigned int k;
    unsigned int width;
    FR_cvector ** rows;
    void (*print)( struct FR_nkctree * );
    void (*free)( struct FR_nkctree * );
    void (*traverse)( struct FR_nkctree * , void * args );
}FR_nkctree;

FR_nkctree * FR_nkctree_new( unsigned int n , unsigned int k );
void FR_nkctree_print( FR_nkctree * self );

#endif

