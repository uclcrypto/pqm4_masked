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
