#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "base.h"
#include "binary_array.h"
#include "xmas_tree.h"
#include "nkctree.h"
#include "time.h"

void help()
{
    printf("syntax: xmas_print [cardinality]\n");
}

int mymain2(int argc, char ** argv){
    int card = 5;
    if( argc == 2 ){
        if(!sscanf(argv[1], "%d", &card)){
            help();
            return -1;
        }
    }

    FR_xmas_tree * xt = FR_xmas_tree_new( card );
    xt->to_latex( xt );

    return 0;
}
