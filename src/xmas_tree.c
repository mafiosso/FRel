#include <stdio.h>
#include <stdlib.h>
#include "binary_array.h"
#include "xmas_tree.h"

static void _xmas_alloc_rows( FR_xmas_tree * self , unsigned int row_count ){
  self->row_count = row_count;
  self->rows = malloc( sizeof( FR_barray* ) * row_count );
}

void FR_xmas_tree_destroy( FR_xmas_tree * self )
{
    for( unsigned int i = 0 ; i < self->row_count ; i++ ){
        FR_barray_destroy( self->rows[i] );
    }

    free( self->rows );
    free( self );
}

void FR_xmas_tree_print( FR_xmas_tree * self ){
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
            for( int bi = 0 ; bi < self->set_size ; bi++ ){
                unsigned char bit = 
                    FR_barray_get( self->rows[r] ,
                                   bi + (si * self->set_size) );

                printf("%d" , (unsigned int)bit );
            }
            printf(" ");
        }
        printf("\n");        
    }
}

void FR_xmas_tree_2latex( FR_xmas_tree * self ){
    unsigned int cols = self->set_size + 1;
    unsigned int space_width = (cols*self->set_size) +
        (cols-1) /* spaces between cols */;

    printf("\% Non human readeable latex output\n");

    for( int r = 0 ; r < self->row_count ; r++ ){
        unsigned int s = (self->rows[r]->len / self->set_size);
        int blank_cols = cols - s;
        int left_blank_cols = blank_cols/2;
        int right_blank_cols = left_blank_cols;

        for( int b = 0 ; b < left_blank_cols ; b++ ){
            printf(" & ");
        }

        for( int si = 0 ; si < s ; si++ ){
            for( int bi = 0 ; bi < self->set_size ; bi++ ){
                unsigned char bit = 
                    FR_barray_get( self->rows[r] ,
                                   bi + (si * self->set_size) );

                printf("%d" , (unsigned int)bit );
            }

            if( si != (s-1) )
                printf(" & ");
        }

        for( int b = 0 ; b < right_blank_cols ; b++ ){
            printf(" & ");
        }

        printf(" \\\\ \n");
    }
}

FR_xmas_tree * FR_xmas_tree_new( unsigned int set_size ){
  FR_xmas_tree * out = calloc( 1 , sizeof( FR_xmas_tree ) );
  out->free = FR_xmas_tree_destroy;
  out->print = FR_xmas_tree_print;
  out->to_latex = FR_xmas_tree_2latex;

  if( set_size <= 1 ){
      unsigned char order1[2] = {0 , 1};
      _xmas_alloc_rows( out , 1 );
      out->rows[0] = FR_barray_new_arr( 2 , order1 );
      out->set_size = 1;
      return out;
  }
  else if( set_size == 2 ){
      unsigned char row1[2] = {1 , 0};
      unsigned char row2[]  = { 0 , 0 ,    0 , 1 ,   1 , 1 };
      _xmas_alloc_rows( out , 2 );
      out->rows[0] = FR_barray_new_arr( 2 , row1 );
      out->rows[1] = FR_barray_new_arr( 6 , row2 );
      out->set_size = 2;
      return out;
  }
  
  /* potential to do it faster */
  FR_xmas_tree * prev = FR_xmas_tree_new( set_size -1 );
  out->set_size = set_size;

  _xmas_alloc_rows( out , prev->row_count*2 );

  /* new rows index */
  unsigned int nr_it = 0;
  
  for( int r = 0 ; r < prev->row_count ; r++ ){
      /* obtaining constant s according to article */
      unsigned int s = (prev->rows[r]->len / prev->set_size);
      /* 1st rule */
      if( s > 1 ){
          out->rows[nr_it] = FR_barray_new( (s-1) * set_size );

          for( int is = 1 ; is < s ; is++ ){
              /* copy tuplet and add 0 starting from 1 */
              FR_barray_copy_tuplet( out->rows[nr_it] , prev->rows[r] ,
                                     (is-1) * set_size ,
                                     prev->set_size * (is) ,
                                     prev->set_size );
              /* a_1 a_2 ... a_s-1 -> a_2.0 ... a_s.0 */
              /* set trailing 0 */
              FR_barray_set( out->rows[nr_it] ,
                             (is*set_size)-1 ,
                             0 );
          }
          nr_it++;
      } /* else omit */

      /* 2nd rule */
      /* s_1 . 0 */
      out->rows[nr_it] = FR_barray_new( (s+1) * set_size );
      FR_barray_copy_tuplet( out->rows[nr_it] , prev->rows[r] , 0 , 0 , prev->set_size );
      FR_barray_set( out->rows[nr_it] , set_size-1 , 0 );
      
      /* s_1 . 1 | s_2 . 1 | ... | s_s-1 . 1 */
      for( int is = 1 ; is <= s ; is++ ){
          FR_barray_copy_tuplet( out->rows[nr_it] , prev->rows[r] , 
                                 (is) * set_size ,
                                 prev->set_size * (is-1) , 
                                 prev->set_size );

          FR_barray_set( out->rows[nr_it] ,
                             (is * set_size) + set_size -1 ,
                             1 );
      }

      nr_it++;
  }
  
  if( nr_it < out->row_count ){
      out->row_count = nr_it;
      /* shrink array - can be optimized by exact computation of array first */
      out->rows = realloc( out->rows , sizeof( FR_barray* ) * out->row_count );
  }
  
  prev->free( prev );
  return out;
}

