#include "tests.h"
#include "masked.h"
#include "masked_utils.h"
#include "masked_representations.h"
#include "hal.h"
#include <stdio.h>
#include "poly.h"

static void report_test(char *msg,int err){
    char buf[128];
    if(err==0){
        sprintf(buf,"%s -> OK",msg);
    }else{
        sprintf(buf,"%s -> /!\\ Error /!\\",msg);
    }
    hal_send_str(buf);
}
unsigned int test_function(){

    int err = 1;
    report_test("test_function",err);
    err = 0;
    report_test("test_function",err);
    return err;
}
unsigned int test_convertions_APoly(){
    poly x;
    StrAPoly strmasked_x,strmasked_y;
    APoly masked_x;

    masked_poly(strmasked_x,&x);
    StrAPoly2APoly(masked_x,strmasked_x);
    APoly2StrAPoly(strmasked_y,masked_x);

    int err = 0;
    for(int i=0;i<KYBER_N;i++){
        for(int d=0;d<NSHARES;d++){
            err += (strmasked_x[d][i] != strmasked_y[d][i]);
        }
    }

    report_test("test_convertions_APoly",err);
    return err;
}

unsigned int test_convertions_bitslice(){
    poly x;
    StrAPoly strmasked_x;
    APoly masked_x,masked_y;
    uint32_t bitslice[COEF_NBITS][NSHARES];
    uint32_t *bitslice_ptr[COEF_NBITS];
    int16_t *dense_ptr[BSSIZE];

    for(int i=0;i<COEF_NBITS;i++){
        bitslice_ptr[i] = &bitslice[i][0];
    }

    for(int i=0;i<BSSIZE;i++){
        dense_ptr[i] = &masked_x[i][0];
    }

    // create dense masked polynomial
    masked_poly(strmasked_x,&x);
    StrAPoly2APoly(masked_x,strmasked_x);

    // map to bitslice
    masked_dense2bitslice(bitslice_ptr,
            dense_ptr,
            COEF_NBITS,
            BSSIZE,
            NSHARES);



    for(int i=0;i<BSSIZE;i++){
        dense_ptr[i] = &masked_y[i][0];
    }
    masked_bitslice2dense(dense_ptr,
            bitslice_ptr,
            COEF_NBITS,
            BSSIZE,
            NSHARES);

    int err = 0;
    for(int i=0;i<COEF_NBITS;i++){
        for(int d=0;d<NSHARES;d++){
            err += (masked_y[i][d] != masked_x[i][d]);
        }
    }

    report_test("test_convertions_bitslice",err);
    return err;
}
