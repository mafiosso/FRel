#ifndef BINARY_ARR
#define BINARY_ARR

typedef struct FR_barray{
  unsigned int len;
  unsigned char * chunk;
}FR_barray;

FR_barray * FR_barray_new( unsigned int bits_count );
FR_barray * FR_barray_new_arr( unsigned int count , unsigned char * arr );
void FR_barray_set( FR_barray * self , unsigned int index , unsigned char value );
unsigned char FR_barray_get( FR_barray * self , unsigned char index );
void FR_barray_copy_tuplet( FR_barray * dest , FR_barray * src , unsigned int dest_idx , unsigned int src_idx , unsigned int len );
void FR_barray_destroy( FR_barray * self );

#endif
