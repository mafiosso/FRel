#include <stdio.h>
#include <stdlib.h>
#include "cvector.h"
#include "nkctree.h"

static void _nkctree_alloc_rows( FR_nkctree * self , 
                                 unsigned int row_count ){
  self->row_count = row_count;
  self->rows = malloc( sizeof( float* ) * row_count );
}

void FR_nkctree_print( FR_nkctree * self ){
    unsigned int w = self->width;
    unsigned int hw = w/2;
    const int chars = 4;

    for( unsigned int r = 0 ; r < self->row_count ; r++ ){
        unsigned int aw = self->rows[r]->size;
        unsigned int skip = (w - aw)/2;
        skip *= (chars+1);

        for( int s = 0 ; s < skip ; s++ ){
            printf(" ");
        }

        for( unsigned int c = 0 ; c < aw ; c++ ){
            if( (c % self->set_size) != (self->set_size-1) ){
                printf("%.2f " , ((float*)self->rows[r]->chunk)[c] );
            }
            else{
                printf("%.2f|" , ((float*)self->rows[r]->chunk)[c] );
            }
        }


        printf("\n");
    }
}

void FR_nkctree_free( FR_nkctree * self ){
    for( int r = 0 ; r < self->row_count ; r++ ){
        FR_cvector_free( self->rows[r]  );
    }
    free( self->rows );
    free( self );
}

static float _a( unsigned int k , unsigned int csz ){

    return (float)k*(1.0/(float)(csz-1));
}

static FR_nkctree * _gen_1case( unsigned int k_chain ){
    FR_nkctree * out = calloc( 1 , sizeof( FR_nkctree ) );
    out->k = k_chain;

    _nkctree_alloc_rows( out , 1 );
    out->rows[0] = FR_cvector_new( sizeof(float) );
    out->set_size = 1;

    for( unsigned int kc = 0 ; kc < k_chain ; kc++ ){
        float val = _a( kc , k_chain );
        FR_cvector_push( out->rows[0] , (void*)&val );
    }

    out->width = k_chain;
    out->print = FR_nkctree_print;
    out->free = FR_nkctree_free;
    return out;
}

static FR_nkctree * _gen_2case( unsigned int k_chain ){
    FR_nkctree * out = calloc( 1 , sizeof( FR_nkctree ) );
    _nkctree_alloc_rows( out , k_chain );
    out->set_size = 2;
    out->k = k_chain;

    unsigned int lt_idx = k_chain -1;
    unsigned int chain_len = 1;

    for( int ik = 0 ; ik < k_chain ; ik ++ ){
        unsigned int r_idx = 0;
        unsigned int l_idx = lt_idx;

        out->rows[ik] = FR_cvector_new( sizeof(float) );

        for( unsigned int icl = 0 ; icl < chain_len ; icl++ ){
            float val = _a( l_idx , k_chain  );
            FR_cvector_push( out->rows[ik] , (void*)&val );
            val = _a( r_idx , k_chain );
            FR_cvector_push( out->rows[ik] , (void*)&val );
            
            if( icl < chain_len/2 ){
                r_idx++;
            }
            else
            {
                l_idx++;
            }
        }

        chain_len += 2;
        lt_idx--;
    }

    out->width = out->rows[ out->row_count-1 ]->size;
    out->print = FR_nkctree_print;
    out->free = FR_nkctree_free;
    return out;
}

static void copy_sigma( int sigma , FR_cvector * orig , 
                        FR_cvector * dest , unsigned int n )
{
    int src_idx = sigma*n;
    for( int c = 0 ; c < n ; c++ ){
        FR_cvector_push( dest , &((float*)orig->chunk)[c+src_idx] );
    }
}

/**/
FR_nkctree * FR_nkctree_new( unsigned int n , unsigned int k ){
    if( n == 1 ){
        return _gen_1case( k );
    }

    if( n == 2 ){
        return  _gen_2case( k );
    }

    FR_nkctree * prev = FR_nkctree_new( n-1 , k );

    FR_cvector ** tree = NULL;
    unsigned int rows_c = 0;
    unsigned int width = 0;    

    for( int r = 0 ; r < prev->row_count ; r++ ){
        int s = prev->rows[r]->size / (n-1);
        printf("row: %d\n" , r );
        printf("s: %d\n" , s );

        if( s < k ){
            tree = realloc( tree , (s + rows_c) * sizeof(FR_cvector*) );

            for( int sig = s-1 ; sig >= 0 ; sig-- ){
                int ks = 0;

                tree[rows_c] = FR_cvector_new( sizeof(float) );

                for( ks = 0 ; ks < (k-sig) ; ks++ ){
                    copy_sigma( sig , prev->rows[r] , tree[rows_c] , n-1 );
                    float aks = _a( ks , k );
                    FR_cvector_push( tree[rows_c] , (void*)&aks );
                }

                for( int sit = sig+1 ; sit < s ; sit++ ){
                    copy_sigma( sit , prev->rows[r],
                                tree[rows_c] , n-1 );
                    float aks = _a( ks , k );
                    FR_cvector_push( tree[rows_c] , (void*)&aks );
                }

                if( tree[rows_c]->size > width ){ width = tree[rows_c]->size; }
                    
                rows_c++;        
            }        
        }
        else{ /* s >= k */
            printf("s >= k\n");

            int a_target_i = 0;
            tree = realloc( tree , (k + rows_c) * sizeof(FR_cvector*) );

            for( int sig_k = k-1 ; sig_k >= 0 ; sig_k-- ){
                tree[rows_c] = FR_cvector_new( sizeof(float) );

                /* left chunk generation until break */
                int a_i = 0;
                for( a_i = 0 ; a_i < a_target_i ; a_i++ ){
                    copy_sigma( sig_k , prev->rows[r] , tree[rows_c] , n-1 );
                    float aks = _a( a_i , k );
                    FR_cvector_push( tree[rows_c] , (void*)&aks );
                }

                //a_i--;
                /* right chunk generation until reach sigma reaches s */
                for( int w_sig_k = sig_k ; w_sig_k < s ; w_sig_k++ ){
                    copy_sigma( w_sig_k , prev->rows[r] , tree[rows_c] , n-1 );
                    float aks = _a( a_i , k );
                    FR_cvector_push( tree[rows_c] , (void*)&aks );
                }                

                a_target_i++;
                if( tree[rows_c]->size > width ){ width = tree[rows_c]->size; }
                rows_c++;
            }
        }
    }

    FR_nkctree_free( prev );
    FR_nkctree * out = malloc( sizeof( FR_nkctree ) );
    out->free = FR_nkctree_free;
    out->print = FR_nkctree_print;
    out->traverse = NULL;

    out->row_count = rows_c;
    out->set_size = n;
    out->k = k;
    /* TODO: establish width in two cases s < k or s >= k 
       it is easy
     */
    out->width = width;    
    out->rows = tree;

    return out;
}
