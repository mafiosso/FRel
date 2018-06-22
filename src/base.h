#ifndef BASE
#define BASE

/* terminal color defs */
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define MAX( a , b ) ((a) < (b) ? (b) : (a))
#define MIN( a , b ) ((a) <= (b) ? (a) : (b))


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
int FR_matrix_eq( FR_matrix * m1, FR_matrix * m2 );

/* fuzzy */
int FR_solution_fuzz_p( FR_matrix * R , FR_matrix * Q , FR_matrix * T );
int FR_minsolution_fuzz_p( FR_matrix * R , FR_matrix * Q ,
                           FR_matrix * T , int kchain);
FR_matrix * FR_eq_gtst_fuzz( FR_matrix * Q_m , FR_matrix * T_m );

FR_matrix * FR_eq_gtst( FR_matrix * Q_m , FR_matrix * T_m );
void FR_matrix_fprint( FR_matrix * self );
float FR_godel_joint( float l , float r );
float FR_godel_residuum( float l , float r );
void FR_matrix_init_f( FR_matrix * self , const float m[self->rows][self->cols]  );

int FR_fuzz_subset_p( FR_matrix * m1 , FR_matrix * m2 );
int FR_fuzz_extra_subset_p( FR_matrix * m1 , FR_matrix * m2 );
FR_ntree * FR_eq_generate_sol_fuzz( FR_matrix * biggest_sol , FR_matrix * Q , 
                                    FR_matrix * T, unsigned k_chain );

unsigned FR_ntree_size( FR_ntree * n );


#endif
