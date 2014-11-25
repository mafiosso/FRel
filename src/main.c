#include <stdio.h>
#include <stdlib.h>

#include "base.h"
#include "binary_array.h"
#include "xmas_tree.h"
#include "nkctree.h"
#include "time.h"

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"


void xmas_print_sol( FR_xmas_tree * self , FR_matrix * Q , FR_matrix * T )
{
    printf("set size: %d\n" , self->set_size );

    unsigned int cols = self->set_size + 1;
    unsigned int space_width = (cols*self->set_size) +
        (cols-1) /* spaces between cols */;
    
    for( int r = 0 ; r < self->row_count ; r++ ){
        unsigned int s = (self->rows[r]->len / self->set_size);
        unsigned int needed_chars = self->rows[r]->len + (s-1);
        unsigned prelude_spaces = (space_width-needed_chars)/2;
        for( int s = 0 ; s < prelude_spaces ; s++ ){
            printf(" ");
        }
        
        for( int si = 0 ; si < s ; si++ ){
            int R[ self->set_size ];
/* COPY R */
            for( int ti = 0 ; ti < self->set_size ; ti++ ){
                R[ti] = FR_barray_get( self->rows[r] , 
                                       ti + (si * self->set_size ) );
            }
            
            FR_matrix * R_m = FR_matrix_new( self->set_size ,
                                             1 , sizeof(int ) );
            FR_matrix_init( R_m , R );
            
            if( FR_eq_solution_p( R_m , Q , T ) ){
                printf("%s" , KGRN);
            }
            
            for( int bi = 0 ; bi < self->set_size ; bi++ ){
                unsigned char bit = 
                    FR_barray_get( self->rows[r] ,
                                   bi + (si * self->set_size) );
                
                printf("%d" , (unsigned int)bit );
            }
            printf(" %s" , KNRM);
        }
        printf("\n");        
    }
}

int main( int argc , char ** argv ){

    const int Q_mat [5][5] = {{1,1,0,1,1},{1,1,0,1,0},{1,0,1,0,1},
                              {1,0,1,0,1},{0,1,1,1,1}};
    const int T_mat [1][5] = {{1, 0, 1 , 0 , 1}};
    
    FR_matrix * Q_m = FR_matrix_new( 5 , 5 , sizeof(int));
    FR_matrix * T_m = FR_matrix_new( 5 , 1 , sizeof(int));
    
    FR_matrix_init( Q_m , Q_mat );
    FR_matrix_init( T_m , T_mat );    

/*
    FR_matrix_iprint( Q_m );
    printf("\n");
    FR_matrix_iprint( T_m );
    printf("\n");

    FR_ntree * nt = FR_eq_solutions( Q_m , T_m );
    if( nt ){
        FR_eq_inspect( nt );
    }
    else{
        printf("no solution\n");
    }

    system( "clear" );
*/
    int x_depth = 5;
    if( argc >= 2 ){
        sscanf( argv[1] , "%d" , &x_depth );
        FR_matrix_free( Q_m );
        FR_matrix_free( T_m );
        Q_m = FR_matrix_new( x_depth , x_depth , sizeof(int) );
        T_m = FR_matrix_new( x_depth , 1 , sizeof( int ) );
        
        clock_t seed = clock();
        if( argc >= 3 ){
            sscanf( argv[2] , "%d" , &seed );
        }
        printf("%u\n" , seed);
        srand( seed );

        int q_arr [x_depth][x_depth];
        int t_arr [1][x_depth];
        for( int y = 0 ; y < x_depth ; y++ ){
            for( int x = 0 ; x < x_depth ; x++ ){
                q_arr[y][x] = rand() % 2;
            }
        }

        for( int x = 0 ; x < x_depth ; x++ ){
            t_arr[0][x] = rand()%2;
        }
        FR_matrix_init( Q_m , q_arr );
        FR_matrix_init( T_m , t_arr );
    }

    printf("\n*** Christmas tree ***\n");

    FR_xmas_tree * xt = FR_xmas_tree_new( x_depth );
//    xt->print( xt );

    xmas_print_sol( xt , Q_m , T_m );

    /* nk ctree */
    FR_nkctree * ct = FR_nkctree_new( 3 , 5 );
    FR_nkctree_print( ct );

/*    for( int r = 0 ; r < ct->row_count ; r++ ){
        for( int l = 0 ; l < ct->rows[r]->size ; l++ ){
            printf("%.2f " , ((float*)ct->rows[r]->chunk)[l] );
        }
        printf("\n");
    }
*/

    return 0;

}
