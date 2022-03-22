#include <stdint.h>
#include <stdlib.h>
#include "bs_gadgets.h"
#include "masked_utils.h"

#define KYBER_MASKING_ORDER (D-1)
#define D NSHARES
#define Q 3329
#define GET_BIT_SHARING(i,j) ((i)*(j)) 

//#define SNI
void RefreshXOR(uint32_t *z,
        uint32_t d){

    uint32_t r;
    for(uint32_t i=0;i<(d-1);i++){
        for(uint32_t j=i+1; j<d;j++){
            r = get_random();
            z[i] ^= (r);
            z[j] ^= (r); 
        }
    }
}

void RefreshIOS(uint32_t *z,
        uint32_t d){
    
    uint32_t r;
    for(uint32_t i=0;i<(d-1);i++){
        for(uint32_t j=i+1; j<d;j++){
            r = get_random();
            z[i] ^= (r);
            z[j] ^= (r); 
        }
    }
}


void RefreshXOR16(uint32_t *z,
        uint32_t d){

    uint32_t r;
    for(uint32_t i=0;i<(d-1);i++){
        for(uint32_t j=i+1; j<d;j++){
            r = rand16();
            z[i] ^= (r);
            z[j] ^= (r); 
        }
    }
}


void FullRefreshXOR(uint32_t *x,
        uint32_t d){

    uint32_t r;
    
    uint32_t x0 = x[0];
    for(uint32_t i=0;i<(d);i++){
        for(uint32_t j=1; j<d;j++){
            r = get_random();
            x0 ^= (r);
            x[j] ^= (r); 
        }
    }
    x[0] = x0;
}

static inline uint32_t pini_and_core(uint32_t a, uint32_t b, uint32_t r) {
    uint32_t temp;
    uint32_t s;
    asm(
    "eor %[temp], %[b], %[r]\n\t"
    "and %[temp], %[a], %[temp]\n\t"
    "bic %[s], %[r], %[a]\n\t"
    "eor %[s], %[s], %[temp]"
    :[s]"=r"(s), [temp]"=&r"(temp) /* outputs, use temp as an arbitrary-location clobber */
    :[a]"r"(a),  [b]"r"(b), [r]"r"(r)  /* inputs */
    );
    return s;
}

static inline void SecAndfix16(const uint32_t *x,
        const uint32_t *y,
        uint32_t *z,
        uint32_t d){
    uint32_t r;
    uint32_t ztmp[d];
    uint32_t xi,yi,ztmpi;
    for(uint32_t i=0;i<d;i++){
        ztmp[i] = x[i] & y[i];
    }
    for(uint32_t i=0;i<(d-1);i++){
        xi = x[i]; yi = y[i]; ztmpi = ztmp[i];
        for(uint32_t j=i+1; j<d;j++){
            r = rand16();
#ifdef SNI
            ztmpi ^= (r) ^ (xi & y[j]);
            ztmp[j] ^= (r) ^ (x[j] & yi);
#else
            ztmpi ^= pini_and_core(xi,y[j],r);
            ztmp[j] ^= pini_and_core(x[j],yi,r);
#endif
        }
        z[i] = ztmpi;
    }
    z[d-1] = ztmp[d-1];
}

void SecAnd16(uint32_t *x,
        uint32_t *y,
        uint32_t *z,
        uint32_t d){
    if(d==1){
        z[0] = y[0] & x[0];
    }else if(d==2){
        uint32_t r = rand16();
        uint32_t tmp;
        tmp = (x[0] & y[0]) ^ (x[0] & y[1]) ^ r;
        z[1] = (x[1] & y[1]) ^ (x[1] & y[0]) ^ r;
        z[0] = tmp;
    }else if(d==3){
        SecAndfix16(x,y,z,3);
    }else if(d==4){
        SecAndfix16(x,y,z,4);
    }else if(d==5){
        SecAndfix16(x,y,z,5);
    }else if(d==6){
        SecAndfix16(x,y,z,6);
    }else if(d==7){
        SecAndfix16(x,y,z,7);
    }else if(d==8){
        SecAndfix16(x,y,z,8);
    }else if(d==9){
        SecAndfix16(x,y,z,9);
    }else if(d==10){
        SecAndfix16(x,y,z,10);
    }else if(d==11){
        SecAndfix16(x,y,z,11);
    }else if(d==12){
        SecAndfix16(x,y,z,12);
    }else if(d==13){
        SecAndfix16(x,y,z,13);
    }else if(d==14){
        SecAndfix16(x,y,z,14);
    }else if(d==15){
        SecAndfix16(x,y,z,15);
    }else if(d==16){
        SecAndfix16(x,y,z,16);
    }else{
        SecAndfix16(x,y,z,d);
    }
}

static inline void SecAndfix(const uint32_t *x,
        const uint32_t *y,
        uint32_t *z,
        uint32_t d){
    uint32_t r;
    uint32_t ztmp[d];
    uint32_t xi,yi,ztmpi;
    for(uint32_t i=0;i<d;i++){
        ztmp[i] = x[i] & y[i];
    }
    for(uint32_t i=0;i<(d-1);i++){
        xi = x[i]; yi = y[i]; ztmpi = ztmp[i];
        for(uint32_t j=i+1; j<d;j++){
            r = get_random();
#ifdef SNI
            ztmpi ^= (r) ^ (xi & y[j]);
            ztmp[j] ^= (r) ^ (x[j] & yi);
#else
            ztmpi ^= pini_and_core(xi,y[j],r);
            ztmp[j] ^= pini_and_core(x[j],yi,r);
#endif
        }
        z[i] = ztmpi;
    }
    z[d-1] = ztmp[d-1];
}

void SecAnd(uint32_t *x,
        uint32_t *y,
        uint32_t *z,
        uint32_t d){
    if(d==1){
        z[0] = y[0] & x[0];
    }else if(d==2){
        uint32_t r = get_random();
        uint32_t tmp;
        tmp = (x[0] & y[0]) ^ (x[0] & y[1]) ^ r;
        z[1] = (x[1] & y[1]) ^ (x[1] & y[0]) ^ r;
        z[0] = tmp;
    }else if(d==3){
        SecAndfix(x,y,z,3);
    }else if(d==4){
        SecAndfix(x,y,z,4);
    }else if(d==5){
        SecAndfix(x,y,z,5);
    }else if(d==6){
        SecAndfix(x,y,z,6);
    }else if(d==7){
        SecAndfix(x,y,z,7);
    }else if(d==8){
        SecAndfix(x,y,z,8);
    }else if(d==9){
        SecAndfix(x,y,z,9);
    }else if(d==10){
        SecAndfix(x,y,z,10);
    }else if(d==11){
        SecAndfix(x,y,z,11);
    }else if(d==12){
        SecAndfix(x,y,z,12);
    }else if(d==13){
        SecAndfix(x,y,z,13);
    }else if(d==14){
        SecAndfix(x,y,z,14);
    }else if(d==15){
        SecAndfix(x,y,z,15);
    }else if(d==16){
        SecAndfix(x,y,z,16);
    }else{
        SecAndfix(x,y,z,d);
    }
}
void SecANDbs(uint32_t *z,
        const uint32_t *a,
        const uint32_t *b){
    uint32_t ztmp[D];
    uint32_t r;
    uint32_t i,j;
    for(i=0;i<D;i++){
        ztmp[i] = a[i] & b[i];
    }
    for(i=0;i<(D-1);i++){
        for(j=i+1; j<D;j++){
            r = get_random();

#ifdef SNI
            ztmp[i] ^= (r) ^ (a[i] & b[j]);
            ztmp[j] ^= (r) ^ (a[j] & b[i]);
#else
            // PINI
            ztmp[i] ^= pini_and_core(a[i],b[j],r); 
            ztmp[j] ^= pini_and_core(a[j],b[i],r); 
#endif
        }
    }
    for(i=0;i<D;i++){
        z[i] = ztmp[i];
    }
}

void SecXORbs(uint32_t *z,
        const uint32_t *a,
        const uint32_t *b){
    for(uint32_t i=0;i<D;i++){
        z[i] = a[i] ^ b[i];
    }
}
