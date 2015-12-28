#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "base.h"
#include "cvector.h"

FR_matrix * FR_matrix_new( unsigned int cols , unsigned int rows , unsigned int bytes_per_col ){
    FR_matrix * out = malloc( sizeof(FR_matrix) );
    out->cols = cols;
    out->rows = rows;
    out->bpc  = bytes_per_col;
    out->matrix_values = malloc(sizeof(void*)*rows);

    for( int i = 0 ; i < rows ; i++ ){
        out->matrix_values[i] = calloc( cols , bytes_per_col );
    }

    return out;
}

void * FR_matrix_get( FR_matrix * m , unsigned int x , unsigned int y ){
    unsigned char * row = m->matrix_values[y];
    return row + (m->bpc*x);
}

void FR_matrix_set( FR_matrix * m , unsigned int x , 
                    unsigned int y , void * val ){
    memcpy(((char*)m->matrix_values[y]) + (x * m->bpc) , val , m->bpc );
}

void FR_matrix_free( FR_matrix * self ){
    for( int y = 0 ; y < self->rows ; y++ ){
        free( self->matrix_values[y] );
    }

    free( self->matrix_values );
    free( self );
}

FR_matrix * FR_matrix_clone( FR_matrix * self ){
    FR_matrix * out = FR_matrix_new( self->cols , self->rows , self->bpc );

    for( int y = 0 ; y < self->rows ; y++ ){
        for( int x = 0 ; x < self->cols ; x++ ){
            memcpy( out->matrix_values[y]+(x*self->bpc) , self->matrix_values[y]+(x*self->bpc) , self->bpc );
        }
    }

    return out;
}

void FR_matrix_iprint( FR_matrix * self ){
    for( int y = 0 ; y < self->rows ; y++ ){
        printf("( ");
        for( int x = 0 ; x < self->cols ; x++ ){
            printf("%d " , ((int*)(self->matrix_values[y]))[x] );
        }
        printf(")\n");                
    }
}

void FR_matrix_fprint( FR_matrix * self ){
    for( int y = 0 ; y < self->rows ; y++ ){
        printf("%% ( ");
        for( int x = 0 ; x < self->cols ; x++ ){
            printf("%f " , ((float*)(self->matrix_values[y]))[x] );
        }
        printf(")\n");                
    }
}

FR_matrix * FR_matrix_transpose( FR_matrix * self ){
    FR_matrix * out = FR_matrix_new( self->rows , self->cols , self->bpc );
    
    for( int y = 0 ; y < self->rows ; y++ ){
        for( int x = 0 ; x < self->cols ; x++ ){
            memcpy( out->matrix_values[x]+(y*self->bpc) , self->matrix_values[y]+(x*self->bpc) , self->bpc );
        }
    }
    
    return out;    
}

void FR_matrix_imul( FR_matrix * m1  , FR_matrix * m2 ,
                     void * out_row , unsigned int row , unsigned int col){
    int * orow = (int*)out_row;
    int ** mv1 = (int**)m1->matrix_values;
    int ** mv2 = (int**)m2->matrix_values;
    
    int sum = 0;

    for( int rc  = 0 ; rc < m1->cols ; rc++ ){
        sum += (mv1[row][rc] * mv2[rc][col]);
    }

    orow[col] = sum;
}

float FR_godel_residuum( float l , float r ){
    if( l <= r ){
        return 1.0f;
    }
    return r;
}

float FR_godel_joint( float l , float r ){
    return (l < r ? l : r );
}

float FR_lukasiewicz_joint( float l , float r ){
    return MAX( l+r-1.0 , 0.0 );
}

float FR_lukasiewicz_residuum( float l , float r ){
    return MIN( 1.0-l+r , 1.0 );
}

void FR_matrix_gjoint( FR_matrix * m1  , FR_matrix * m2 ,
                       void * out_row , unsigned int row , 
                       unsigned int col)
{
    float * orow = (float*)out_row;
    float ** mv1 = (float**)m1->matrix_values;
    float ** mv2 = (float**)m2->matrix_values;

    float max = 0.0f;

    for( int rc  = 0 ; rc < m1->cols ; rc++ ){
        float lmax = FR_lukasiewicz_joint( mv1[row][rc] , mv2[rc][col] );
        max = ( max < lmax ? lmax : max );
    }

    orow[col] = max;
}

void FR_matrix_gresiduum( FR_matrix * m1  , FR_matrix * m2 ,
                          void * out_row , unsigned int row , 
                          unsigned int col)
{
    float * orow = (float*)out_row;
    float ** mv1 = (float**)m1->matrix_values;
    float ** mv2 = (float**)m2->matrix_values;

    float min = 1.0f;

    for( int rc  = 0 ; rc < m1->cols ; rc++ ){
        float lmin = FR_lukasiewicz_residuum( mv1[row][rc] , mv2[rc][col] );
        min = ( min > lmin ? lmin : min );
    }

    orow[col] = min;
}

/* compute composition of crisp relations */
void FR_matrix_irelcomp( FR_matrix * m1  , FR_matrix * m2 ,
                         void * out_row , unsigned int row , unsigned int col){
    FR_matrix_imul( m1  , m2 ,  out_row , row , col);
    
    int * orow = (int*)out_row;
    orow[col] = (orow[col] > 0 ? 1 : 0);
}

void FR_matrix_iimply( FR_matrix * m1  , FR_matrix * m2 ,
                       void * out_row , unsigned int row , unsigned int col){
    int * orow = (int*)out_row;
    int ** mv1 = (int**)m1->matrix_values;
    int ** mv2 = (int**)m2->matrix_values;
    
    for( int rc  = 0 ; rc < m1->cols ; rc++ ){
        if( mv1[row][rc] > mv2[rc][col] ){
            orow[col] = 0;
            return;
        }
    }

    orow[col] = 1;
}

FR_matrix * FR_matrix_o( FR_matrix * m1 , FR_matrix * m2 , 
                         void (*oper)( FR_matrix * , FR_matrix * ,
                                       void * , unsigned int, unsigned int) )
{
    FR_matrix * out = FR_matrix_new( m2->cols , m1->rows , m1->bpc );

    for( int r1 = 0 ; r1 < m1->rows ; r1++ ){
        for( int c2 = 0 ; c2 < m2->cols ; c2++ ){
            oper( m1 , m2 , out->matrix_values[r1] , r1 , c2 );
        }
    }

    return out;
}

int FR_matrix_eq( FR_matrix * m1, FR_matrix * m2 ){
    for( unsigned r = 0 ; r < m1->rows ; r++ ){
        for( unsigned c = 0 ; c < m->cols ; c++ ){
            void * m1_cell_address = ((char*)m1->matrix_values[r])+(c*m1->bpc);
            void * m2_cell_address = ((char*)m2->matrix_values[r])+(c*m1->bpc);
            if( memcmp( m1_cell_address, m2_cell_address, m1->bpc ) )
                return 0;
        }
    }
    return 1;
}

/* R o Q = T */
int FR_eq_solution_p( FR_matrix * R , FR_matrix * Q , FR_matrix * T ){
    FR_matrix * RoQ = FR_matrix_o( R , Q , FR_matrix_irelcomp );
    
    for( int r = 0 ; r < RoQ->rows ; r++ ){
        for( int c = 0 ; c < RoQ->cols ; c++ ){
            if( memcmp( RoQ->matrix_values[r]+(c*RoQ->bpc) ,
                        T->matrix_values[r]+(c*RoQ->bpc) , RoQ->bpc ) )
                return 0;
        }
    }

    return 1;
}

int FR_solution_fuzz_p( FR_matrix * R , FR_matrix * Q , FR_matrix * T )
{
    FR_matrix * RoQ = FR_matrix_o( R , Q , FR_matrix_gjoint );

    for( int r = 0 ; r < RoQ->rows ; r++ ){
        for( int c = 0 ; c < RoQ->cols ; c++ ){
            if( memcmp( RoQ->matrix_values[r]+(c*RoQ->bpc) ,
                        T->matrix_values[r]+(c*RoQ->bpc) , RoQ->bpc ) )
            {
                FR_matrix_free( RoQ );
                return 0;
            }
        }
    }

    FR_matrix_free( RoQ );
    return 1;
}

int FR_minsolution_fuzz_p( FR_matrix * R , FR_matrix * Q ,
                           FR_matrix * T , int kchain)
{
    if( !FR_solution_fuzz_p( R , Q , T ) )
        return 0;

    float step = 1.0f / (kchain-1);

    for( int i = 0 ; i < R->cols ; i++ ){
        float * cur = FR_matrix_get( R , i , 0 );
        float prev = *cur;
        if( (*cur - step) < 0.0f ){
            continue;
        }

        *cur = *cur - step;
     

        if( FR_solution_fuzz_p( R , Q , T ) )
        {
            *cur = prev;
            return 0;
        }
        
        *cur = prev;
    }

    return 1;   
}

FR_ntree * FR_eq_generate_sol( FR_matrix * biggest_sol , FR_matrix * Q , FR_matrix * T ){
    /* TODO: make it general not only for 1xN matrices ! */
    /* first check if biggest_sol is solution*/
    if( !FR_eq_solution_p( biggest_sol , Q , T ) ){
        return NULL;
    }

    FR_ntree * out = malloc( sizeof(FR_ntree) );
    out->childs_count = 0;
    out->childs = NULL;
    out->value  = biggest_sol;

    for( int x = 0 ; x < biggest_sol->cols ; x++ ){
        if( ((int*)(biggest_sol->matrix_values[0]))[x] == 0 ){
            continue;            
        }

        FR_matrix * subsol = FR_matrix_clone( biggest_sol );
        ((int*)(subsol->matrix_values[0]))[x] = 0;
        FR_ntree * child = FR_eq_generate_sol( subsol , Q , T );
        
        if( child != NULL ){
            out->childs = realloc( out->childs , 
                                   sizeof( FR_ntree**) * 
                                   (out->childs_count+1) );
            out->childs[out->childs_count] = child;
            out->childs_count++;
        }
        else{
            FR_matrix_free( subsol );
        }
    }

    return out;
}

FR_ntree * FR_eq_generate_sol_fuzz( FR_matrix * biggest_sol , FR_matrix * Q , 
                                    FR_matrix * T, unsigned k_chain ){
    /* TODO: make it general not only for 1xN matrices ! */
    /* first check if biggest_sol is solution*/
    if( !FR_eq_solution_fuzz_p( biggest_sol , Q , T ) ){
        return NULL;
    }

    float step = 1.0f / (float)k_chain;

    FR_ntree * out = malloc( sizeof(FR_ntree) );
    out->childs_count = 0;
    out->childs = NULL;
    out->value  = biggest_sol;

    for( int x = 0 ; x < biggest_sol->cols ; x++ ){
        if( ((int*)(biggest_sol->matrix_values[0]))[x] == 0 ){
            continue;            
        }

        FR_matrix * subsol = FR_matrix_clone( biggest_sol );
        ((float*)(subsol->matrix_values[0]))[x] -= step;
        FR_ntree * child = FR_eq_generate_sol_fuzz( subsol , Q , T, k_chain );
        
        if( child != NULL ){
            out->childs = realloc( out->childs , 
                                   sizeof( FR_ntree**) * 
                                   (out->childs_count+1) );
            out->childs[out->childs_count] = child;
            out->childs_count++;
        }
        else{
            FR_matrix_free( subsol );
        }
    }

    return out;
}

void FR_ntree2set( FR_ntree * n , FR_cvector * v){
    unsigned childs = n->childs_count;
    unsigned dups = 0;

    for( unsigned c = 0 ; c < childs; c++ ){
        /* filter out duplicates */
        FR_matrix * v = (FR_matrix*)n->value;
        for( unsigned i = 0; i < v->size; i++ ){
            FR_matrix * curval = NULL;
            FR_cvector_nth( v, i, &curval );
            if( FR_matrix_eq( v, curval )){
                dups++;
                goto duplicate;
            }            
        }
        
        FR_cvector_push( v, n->value );
        dups += FR_nkctree2set( n->childs[c], v );
        duplicate:
    }

    return dups;
}

FR_matrix * FR_eq_gtst_fuzz( FR_matrix * Q_m , FR_matrix * T_m ){
    FR_matrix * invT = FR_matrix_transpose( T_m );
    FR_matrix * bs   = FR_matrix_o( Q_m , invT , FR_matrix_gresiduum );
    FR_matrix * invbs = FR_matrix_transpose( bs );

    return invbs;
}

int FR_fuzz_subset_p( FR_matrix * m1 , FR_matrix * m2 ){

    for( unsigned y = 0 ; y < m1->rows ; y++ ){
        for( unsigned x = 0 ; x < m1->cols ; x++ ){
            float v1 = *((float*)FR_matrix_get( m1 , x , y ));
            float v2 = *((float*)FR_matrix_get( m2 , x , y ));
            if( v1 > v2 ){ return 0; }
        }
    }

    return 1;
}

int FR_fuzz_extra_subset_p( FR_matrix * m1 , FR_matrix * m2 ){
    for( unsigned y = 0 ; y < m1->rows ; y++ ){
        for( unsigned x = 0 ; x < m1->cols ; x++ ){
            float v1 = *((float*)FR_matrix_get( m1 , x , y ));
            float v2 = *((float*)FR_matrix_get( m2 , x , y ));
            if( !((v1 == v2) ||
                  (v1 == 0.0)) ){
                return 0;
            }
        }
    }

    return 1;
}


/* solves Q o R = T */
FR_ntree * FR_eq_solutions( FR_matrix * Q_matrix , FR_matrix * T_matrix ){
    /* find the greatest solution first */
    /* inv(Q o inv(T)) */
    FR_matrix * invT = FR_matrix_transpose( T_matrix );

    FR_matrix * bs = FR_matrix_o( Q_matrix , invT , FR_matrix_iimply );
    FR_matrix * invbs = FR_matrix_transpose( bs );

    FR_matrix_iprint( invbs );

    return FR_eq_generate_sol( invbs , Q_matrix , T_matrix );
}

void FR_eq_inspect( FR_ntree * lattice ){
    printf("< ");
    FR_matrix_iprint( (FR_matrix*)lattice->value );

    for( int c = 0 ; c < lattice->childs_count ; c++ ){
        FR_eq_inspect( lattice->childs[c] );
    }

    printf(">\n");
}

void FR_matrix_init( FR_matrix * self , const int m[self->rows][self->cols]  ){
    for( int y = 0 ; y < self->rows ; y++ ){
        for( int x = 0 ; x < self->cols ; x++ ){
            memcpy( self->matrix_values[y]+(x*self->bpc) , &(m[y][x]) , self->bpc );
        }
    }
}

void FR_matrix_init_f( FR_matrix * self , const float m[self->rows][self->cols]  ){
    for( int y = 0 ; y < self->rows ; y++ ){
        for( int x = 0 ; x < self->cols ; x++ ){
            memcpy( self->matrix_values[y]+(x*self->bpc) , &(m[y][x]) , self->bpc );
        }
    }
}

