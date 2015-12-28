#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cvector.h"

FR_cvector * FR_cvector_new( unsigned elem_size ){
    FR_cvector * out = calloc( 1 , sizeof( FR_cvector ) );
    out->elem_size = elem_size;
    return out;
}

void FR_cvector_free( FR_cvector * self ){
    free( self->chunk );
    free( self );
}

void FR_cvector_nth( FR_cvector * self, unsigned nth, void * out ){
    memcpy(out , self+(nth*self->elem_size) , self->elem_size);
}

void FR_cvector_push( FR_cvector * self , void * item ){
    if( self->allocd <= self->size ){
        if( !self->allocd ){
            self->chunk = malloc( 1 * self->elem_size );
            self->allocd = 1;
        }
        else{
            self->allocd *= 2;
            self->chunk = realloc( self->chunk , 
                                   self->allocd * self->elem_size );
        }
    }

    memcpy( self->chunk + (self->size * self->elem_size),
            item , self->elem_size );

    self->size++;
}
