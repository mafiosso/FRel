#ifndef BASE
#define BASE

typedef struct FR_relation{
    /* TODO: unsigned int dim; i.e.: A1 x ... x AN have dimension N */
    void * data_ptr;

}FR_relation;

typedef struct FR_ntree{
    unsigned int childs_count;
    struct FR_ntree ** childs;
    void * value;
}FR_ntree;

typedef struct FR_matrix{
    unsigned int cols;
    unsigned int rows;
    unsigned int bpc;
    void ** matrix_values;
}FR_matrix;


FR_matrix * FR_matrix_new( unsigned int cols , unsigned int rows , unsigned int bytes_per_col );
void FR_matrix_free( FR_matrix * self );
FR_matrix * FR_matrix_clone( FR_matrix * self );
void FR_matrix_iprint( FR_matrix * self );
FR_matrix * FR_matrix_transpose( FR_matrix * self );
void FR_matrix_imul( FR_matrix * m1  , FR_matrix * m2 ,
                     void * out_row , unsigned int row , unsigned int col);
void FR_matrix_irelcomp( FR_matrix * m1  , FR_matrix * m2 ,
                         void * out_row , unsigned int row , unsigned int col);
void FR_matrix_iimply( FR_matrix * m1  , FR_matrix * m2 ,
                       void * out_row , unsigned int row , unsigned int col);
FR_matrix * FR_matrix_o( FR_matrix * m1 , FR_matrix * m2 , 
                         void (*oper)( FR_matrix * , FR_matrix * , void * , unsigned int, unsigned int) );
int FR_eq_solution_p( FR_matrix * R , FR_matrix * Q , FR_matrix * T );
FR_ntree * FR_eq_generate_sol( FR_matrix * biggest_sol , FR_matrix * Q , FR_matrix * T );
FR_ntree * FR_eq_solutions( FR_matrix * Q_matrix , FR_matrix * T_matrix );
void FR_eq_inspect( FR_ntree * lattice );
void FR_matrix_init( FR_matrix * self , const int m[self->rows][self->cols]  );

#endif
