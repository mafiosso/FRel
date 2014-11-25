#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "binary_array.h"

FR_barray * FR_barray_new( unsigned int bits_count )
{
    unsigned sz = (bits_count/8) + ((bits_count%8) > 0);
    assert( sz > 0 );
    
    unsigned char * chunk = calloc( sz ,
                                    sizeof( unsigned char ) );
    FR_barray * out = malloc( sizeof( FR_barray ) );
    out->chunk = chunk;
    out->len = bits_count;
    return out;
}

FR_barray * FR_barray_new_arr( unsigned int count , unsigned char * arr ){
  FR_barray * out = FR_barray_new( count );

  for( int i = 0 ; i < count ; i++ ){
    FR_barray_set( out , i , arr[i] );
  }

  return out;
}

void FR_barray_set( FR_barray * self , unsigned int index , unsigned char value ){

    assert( self->len > index );
    
    unsigned int c_ptr = index/8;
    unsigned int r_ptr = index%8;

    unsigned char mask = 1 << r_ptr;

    if( value ){
        /* set bit */
        self->chunk[ c_ptr ] |= mask;
    }
    else{
        /* reset bit */
        self->chunk[ c_ptr ] &= (~mask);
    }
}

/* returns 1 if set 0 otherwise */
unsigned char FR_barray_get( FR_barray * self , unsigned char index ){
    assert( self->len > index );
    
    unsigned int c_ptr = index/8;
    unsigned int r_ptr = index%8;
    
    unsigned char mask = 1 << r_ptr;
    return ((self->chunk[ c_ptr ] & mask) > 0);
}

void FR_barray_copy_tuplet( FR_barray * dest , FR_barray * src , unsigned int dest_idx , unsigned int src_idx , unsigned int len ){
  for( int i = 0 ; i < len ; i++ ){
    FR_barray_set( dest , dest_idx + i , FR_barray_get( src , src_idx+i ) );
  }
}

void FR_barray_destroy( FR_barray * self ){
    free( self->chunk );
    free( self );
}
