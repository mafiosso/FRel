#ifndef XMAS
#define XMAS

typedef struct FR_xmas_tree{
  unsigned int set_size;
  unsigned int row_count;
  FR_barray ** rows;
  void (*print)( struct FR_xmas_tree * );
  void (*to_latex)( struct FR_xmas_tree * );
  void (*free) ( struct FR_xmas_tree * );
  void (*traverse_hensl) ( struct FR_xmas_tree * , unsigned char * (*hensl_func)(void *) );
}FR_xmas_tree;


FR_xmas_tree * FR_xmas_tree_new( unsigned int set_size );
#endif
