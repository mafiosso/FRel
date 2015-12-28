#ifndef CVECTOR_H
#define CVECTOR_H

typedef struct FR_cvector{
    void * chunk;
    unsigned char elem_size;
    unsigned int allocd;
    unsigned int size;
}FR_cvector;

FR_cvector * FR_cvector_new( unsigned elem_size );
void FR_cvector_push( FR_cvector * self , void * item );
void FR_cvector_free( FR_cvector * self );
void FR_cvector_nth( FR_cvector * self, unsigned nth, void * out );
#endif
