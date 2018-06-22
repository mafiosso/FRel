#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "base.h"
#include "binary_array.h"
#include "xmas_tree.h"
#include "nkctree.h"
#include "time.h"

unsigned FR_ntree2set( FR_ntree * n , FR_cvector * v);

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


FR_matrix * fmat_rand( int cols , int rows , int k ){
    float step = 1.0f / (k-1);

    FR_matrix * out = FR_matrix_new( cols , rows , sizeof(float) );

    for( int y = 0 ; y < rows ; y++ ){
        for( int x = 0 ; x < cols ; x++ ){
            ((float**)out->matrix_values)[y][x] = (rand()%(k-1)) * step;
        }
    }

    return out;
}

#define N 5
#define K 3
#define MIN_SOLS_COUNT 80
#define MIN_MIN_SC 4

FR_matrix ** find_example(){
    FR_matrix ** out = malloc(sizeof(FR_matrix*) * 2 );
    FR_nkctree * nk = FR_nkctree_new( N , K );

    unsigned count_sols = 0;
    unsigned min_count = 0;
    int tries = 0;

    do{
        tries++;
        if( (tries % 50000) == 0 ){
            // printf("tries: %d\n" , tries );
        }


        FR_matrix * Q = fmat_rand( 4 , 5 , K );
        FR_matrix * T = fmat_rand( 4 , 1 , K );

        float val = 1.0f;
        FR_matrix_set(Q , 2 , 1, &val );
        FR_matrix_set(Q , 1 , 0, &val );
        val = 0.0f;
        FR_matrix_set(T , 3 , 0 , &val );

        count_sols = FR_nkctree_sols_count( nk , Q , T , &min_count );
        if( min_count >= MIN_MIN_SC  && count_sols >= MIN_SOLS_COUNT ){

/*            for( int y = 0 ; y < Q->rows ; y++ ){
                for( int x = 0 ; x < Q->cols ; x++ ){
                    float val = *(float*)(FR_matrix_get( Q , x , y ));
                    float val2 = *(float*)(FR_matrix_get( T , x , 0 ));
                    if( val == 1.0f || val2 == 1.0f ){
                        goto succ;
                    }
                }
            }
            FR_matrix_free( Q );
            FR_matrix_free( T );
            continue;
            
            
        succ:
*/
            out[0] = Q;
            out[1] = T;
            return out;
        }


        FR_matrix_free( Q );
        FR_matrix_free( T );
    }while( 1 );

    FR_nkctree_free( nk );

    free( out );
    return NULL;
}

FR_matrix ** find_example2(){
    FR_matrix ** out = malloc(sizeof(FR_matrix*) * 2 );
    FR_nkctree * nk = FR_nkctree_new( 4 , 3 );

    unsigned count_sols = 0;
    unsigned min_count = 0;
    int tries = 0;

    do{
        tries++;

        FR_matrix * Q = fmat_rand( 4 , 5 , K );
        FR_matrix * T = fmat_rand( 4 , 1 , K );

        unsigned mcount = 6;

        if( (tries % 50000) == 0 ){
            printf("tries: %d, mcount: %d\n" , tries, nk->row_count );
        }

        count_sols = FR_nkctree_sols_count( nk , Q , T , &min_count );
        if( min_count >= mcount ){
            printf("mcount: %d\n", mcount);
            printf("count sols: %d, min_sols: %d\n", count_sols, min_count);

            out[0] = Q;
            out[1] = T;
            printf("out===>\n");
            FR_matrix_fprint( Q );
            printf("\n");
            FR_matrix_fprint( T );
            printf("<===end\n");

            return out;
        }

        FR_matrix_free( Q );
        FR_matrix_free( T );
    }while( 1 );

    FR_nkctree_free( nk );

    free( out );
    return NULL;
}

void binary_pattern(){
    const float Q_mat [5][4] = {{1, 1, 0, 1},
                                {1, 0, 1, 1},
                                {1, 1, 1, 1},
                                {1, 1, 1, 1},
                                {1, 1, 1, 1}};
    const float T_mat [1][4] = {{1,1,1,1}};

    FR_matrix * Q_fm = FR_matrix_new( 4 , 5 , sizeof(float) );
    FR_matrix * T_fm = FR_matrix_new( 4 , 1 , sizeof(float) );

    FR_matrix_init_f( Q_fm , Q_mat );
    FR_matrix_init_f( T_fm , T_mat );

    FR_nkctree * ct = FR_nkctree_new( 5 , 2 );
    
    FR_matrix_fprint( Q_fm );
    printf("\n" );
    FR_matrix_fprint( T_fm );

    FR_nkctree_2tikz( ct , 
                      Q_fm ,
                      T_fm );
    
}

void test_stack(){
    const float Q_fmat [5][4] = {{0.5,1.0,0.5,0.5},
                                 {0.5,0.5,1.0,0.0},
                                 {0,0.5,0.5,0},
                                 {0,0.5,0.5,0},
                                 {0,0.5,0.5,0}};

    const float T_fmat [1][4] = {{0.0, 0.5, 0.5 , 0.0}};
    
    FR_matrix * Q_fm = FR_matrix_new( 4 , 5 , sizeof(float) );
    FR_matrix * T_fm = FR_matrix_new( 4 , 1 , sizeof(float) );

    FR_matrix_init_f( Q_fm , Q_fmat );
    FR_matrix_init_f( T_fm , T_fmat );

    FR_cvector * cv = FR_cvector_new( sizeof(FR_matrix*) );

    FR_matrix * Q_tst = NULL;
    FR_matrix * T_tst = NULL;

    FR_cvector_push( cv, &Q_fm );
    FR_cvector_push( cv, &T_fm );

    FR_cvector_nth( cv, 0, &Q_tst );
    FR_cvector_nth( cv, 1, &T_tst );
    

    FR_matrix_fprint( Q_tst );
    FR_matrix_fprint( T_tst );

    printf("DONE test\n");
}

void FR_ntree_free( FR_ntree * nt ){
    for( unsigned c = 0 ; c < nt->childs_count ; c++ ){
        FR_ntree_free( nt->childs[c] );
    }

    free( nt->childs );
    free( nt );
}

/*
Memory intensive!

% ( 0.000000 0.000000 0.500000 0.000000 0.000000 0.500000 0.500000 0.000000 )
% ( 0.500000 0.000000 0.500000 0.500000 0.000000 0.000000 0.000000 0.000000 )
% ( 0.500000 0.500000 0.500000 0.000000 0.000000 0.000000 0.500000 0.500000 )
% ( 0.000000 0.500000 0.000000 0.500000 0.500000 0.000000 0.500000 0.500000 )
% ( 0.500000 0.000000 0.500000 0.500000 0.500000 0.500000 0.000000 0.500000 )
% ( 0.500000 0.500000 0.000000 0.000000 0.500000 0.000000 0.000000 0.000000 )
% ( 0.500000 0.000000 0.000000 0.500000 0.000000 0.500000 0.000000 0.000000 )
% ( 0.500000 0.500000 0.500000 0.000000 0.500000 0.000000 0.500000 0.000000 )
% ( 0.000000 0.500000 0.000000 0.000000 0.000000 0.000000 0.500000 0.500000 )
*/

unsigned duplicates_count( FR_matrix * Q, FR_matrix * T, unsigned k_chain ){
    FR_matrix * g_sol = FR_eq_gtst_fuzz( Q , T );
/*  FR_matrix_fprint( Q );
    FR_matrix_fprint( T );
    FR_matrix_fprint( g_sol );*/


    FR_ntree * t = FR_eq_generate_sol_fuzz( g_sol, Q , T, k_chain );
//    printf("ntree generated\n");

    FR_cvector * cv = FR_cvector_new( sizeof(FR_matrix*) );
    unsigned dups = FR_ntree2set( t, cv );

    for( unsigned d = 0; d < cv->size; d++ ){
        FR_matrix * sol = NULL;
        FR_cvector_nth( cv, d, &sol );
        if( sol ){
            FR_matrix_free( sol );
        }
    }

    FR_cvector_free( cv );
    //FR_matrix_free( g_sol );
    FR_ntree_free( t );
    return dups;
}

int * make_statistics( unsigned min_sols_count, unsigned n, unsigned k )
{
    const unsigned tries_cnt = 10000;
    unsigned d = n;
    FR_nkctree * nk = FR_nkctree_new( n , k );
    int * out = calloc(tries_cnt, sizeof(unsigned));

    unsigned count_sols = 0;
    unsigned min_count = 0;
    int tries = 0;

    FR_matrix * g_Q = NULL;
    FR_matrix * g_T = NULL;
    unsigned dup_max = 0;

    printf("[");

    do{
        tries++;
        if( (tries % 50000) == 0 ){
            printf("tries: %d\n" , tries );
        }


        FR_matrix * Q = fmat_rand( d , d , k );
        FR_matrix * T = fmat_rand( d , 1 , k );

        float val = 1.0f;

        count_sols = FR_nkctree_sols_count( nk , Q , T , &min_count );

        out[tries] = -1;

        if( count_sols > 0 ){
            unsigned dups = duplicates_count( Q, T, k );
            assert( dups <= count_sols );
//            printf("\{dups: %d, solutions: %d, min_sols: %d\},\n", dups, count_sols, min_count);
            out[tries] = dups;
            if( dups > dup_max ){
                dup_max = dups;
                printf("current max dup: %d\n", dup_max);
                printf("current sols: %d\n", count_sols);

                FR_matrix_fprint( Q );
                FR_matrix_fprint( T );
                
                g_Q = Q;
                g_T = T;
            }
            else{
                FR_matrix_free( Q );
                FR_matrix_free( T );
            }
        }

        if( !(tries % 1000) ){
            printf("try: %d\n", tries);
        }

        /*if( tries > 4000 ){
            printf("\n");
            FR_matrix_fprint(Q);
            FR_matrix_fprint(T);
            printf("\n");
            }*/


    }while( tries < tries_cnt );

    printf("*** result ***\n");
    FR_matrix_fprint( g_Q );
    printf("\n");
    FR_matrix_fprint( g_T );

    printf("max duplicates: %d\n", dup_max );

    printf("%c]", 8); /* backspacing, */

    FR_nkctree_free( nk );
    return out;
}

FR_matrix * fmat_zero(unsigned cols, unsigned rows){
    FR_matrix * out = FR_matrix_new( cols , rows , sizeof(float));
    float zero = 0.0f;

    for(unsigned r = 0 ; r < rows ; r++){
        for(unsigned c = 0; c < cols ; c++){
            FR_matrix_set( out , c , r, &zero );
        }
    }
    return out;
}

int main(){
    FR_matrix ** exmpl = find_example2();

    /* duplicates statistics */
//    int * ret = make_statistics( 4 , 4, 2 );
    for( unsigned k = 3; k < 7 ; k++ ){
        if( k == 4 ){ continue; }
        
        printf(">>>>\n");

        FR_nkctree * nk = FR_nkctree_new( 4 , k );
        FR_matrix * Q = exmpl[0];
        FR_matrix * T = exmpl[1];
        unsigned min_count = 0;
        unsigned count_sols = FR_nkctree_sols_count( nk , Q , T , &min_count );
        FR_matrix * g_sol = FR_eq_gtst_fuzz( Q , T );
        FR_ntree * t = FR_eq_generate_sol_fuzz( g_sol, Q , T, k );
        unsigned dup_sols = FR_ntree_size(t );
        unsigned dups = dup_sols - count_sols;
        
        printf("<<<< dim: %d, k: %d, sols: %d, min_sols: %d dups: %d\n", 4, k, count_sols, min_count, dups);            
    }
    
    return 0;
}

int dupmain(){
    const float Q_fmat [5][4] = {{0.5,1.0,0.5,0.5},
                                 {0.5,0.5,1.0,0.0},
                                 {0,0.5,0.5,0},
                                 {0,0.5,0.5,0},
                                 {0,0.5,0.5,0}};

    const float T_fmat [1][4] = {{0.0, 0.5, 0.5 , 0.0}};
    
    FR_matrix * Q_fm = FR_matrix_new( 4 , 5 , sizeof(float) );
    FR_matrix * T_fm = FR_matrix_new( 4 , 1 , sizeof(float) );

    FR_matrix_init_f( Q_fm , Q_fmat );
    FR_matrix_init_f( T_fm , T_fmat );

    printf("cmp: %d\n", FR_matrix_eq(Q_fm, Q_fm) );
    printf("cmp: %d\n", FR_matrix_eq(T_fm, T_fm) );

    printf("dups count: %d\n", duplicates_count( Q_fm, T_fm, 3 ));
    return 0;
}


int mymain( int argc , char ** argv ){

    const int Q_mat [5][5] = {{1,1,0,1,1},{1,1,0,1,0},{1,0,1,0,1},
                              {1,0,1,0,1},{0,1,1,1,1}};
    const int T_mat [1][5] = {{1, 0, 1 , 0 , 1}};

    const float Q_fmat [5][4] = {{0.5,1.0,0.5,0.5},
                                 {0.5,0.5,1.0,0.0},
                                 {0,0.5,0.5,0},
                                 {0,0.5,0.5,0},
                                 {0,0.5,0.5,0}};

    const float T_fmat [1][4] = {{0.0, 0.5, 0.5 , 0.0}};
    
    FR_matrix * Q_m = FR_matrix_new( 5 , 5 , sizeof(int));
    FR_matrix * T_m = FR_matrix_new( 5 , 1 , sizeof(int));

    FR_matrix * Q_fm = FR_matrix_new( 4 , 3 , sizeof(float) );
    FR_matrix * T_fm = FR_matrix_new( 4 , 1 , sizeof(float) );
    
    FR_matrix_init( Q_m , Q_mat );
    FR_matrix_init( T_m , T_mat );

    FR_matrix_init_f( Q_fm , Q_fmat );
    FR_matrix_init_f( T_fm , T_fmat );

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
        printf("%u\n" , (unsigned)seed);
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



//    printf("\n*** Christmas tree ***\n");

//    FR_xmas_tree * xt = FR_xmas_tree_new( x_depth );
//    xt->print( xt );

    //  xmas_print_sol( xt , Q_m , T_m );

    /* nk ctree */

//    FR_nkctree_2latex( ct );

/*    for( int r = 0 ; r < ct->row_count ; r++ ){
        for( int l = 0 ; l < ct->rows[r]->size ; l++ ){
            printf("%.2f " , ((float*)ct->rows[r]->chunk)[l] );
        }
        printf("\n");
    }
*/


//    FR_nkctree * ct = FR_nkctree_new( N , K );

/*    FR_matrix * gt = FR_eq_gtst( Q_fm , T_fm );
    printf("greatest solution: \n");
    FR_matrix_fprint( gt );

    FR_matrix * nil = FR_matrix_new( 3 , 1 , sizeof(float) );

    if( FR_solution_fuzz_p( nil , Q_fm , T_fm ) ){
        printf("IS SOLUTION!");

        if( !FR_minsolution_fuzz_p( gt , Q_fm , T_fm , 5 ) ){
            printf("\nISNT MINIMAL!\n");
        }
    }
    else{
        printf("ISNT SOLUTION!\n" );
        }*/

/*
    
    FR_matrix ** good = find_example();
    Q_fm = good[0];
    T_fm = good[1];

    srand( 12348463 );

    FR_matrix_fprint( Q_fm );
    printf("\n" );

    FR_matrix_fprint( T_fm );

    //FR_nkctree_sols( ct , Q_fm, T_fm );

    FR_nkctree_2tikz( ct , 
                      Q_fm ,
                      T_fm );
*/

//    binary_pattern();

    FR_xmas_tree * fx = FR_xmas_tree_new( 6 );
    FR_xmas_tree_print( fx );

    return 0;

}
