#include "masked.h"
#include "gadgets_legacy.h"
#include "masked_utils.h"
static uint32_t deriveW(uint32_t k);
static uint32_t deriveW(uint32_t k){
    k -= 1;
    uint32_t r = 0;
    uint32_t hw = k & 0x1;
    k >>=1;
    while(k>0){
        hw += (k & 0x1);
        k >>= 1;
        r++;
    }
    return ((r + (hw != 1))-1); 
}
void my_SecAdd(uint32_t *x,
        uint32_t *y,
        uint32_t *z,
        uint32_t d,
        uint32_t k){
    
    uint32_t mask = (1<<k) -1;
    if(k==32){
        mask = 0xffffffff;
    }
    uint32_t W = deriveW(k);
    uint32_t p[d];
    uint32_t g[d];
    uint32_t tmp[d];
    uint32_t a[d];

    uint32_t pow;
    for(uint32_t i = 0;i<d;i++){
        p[i] = x[i] ^ y[i];
    }

    if(k > 16){SecAnd(x,y,g,d);}else{SecAnd16(x,y,g,d);}
    for(uint32_t j=0; j<W ;j++){
        pow = 1<<j;
        
        for(uint32_t i =0;i<d;i++){
            a[i] = (g[i] << pow)&mask;
        }

        if(k>16){SecAnd(a,p,tmp,d);}else{SecAnd16(a,p,tmp,d);};

        for(uint32_t i=0;i<d;i++){
            a[i] = tmp[i];
            g[i] ^= a[i];
        }

        for(uint32_t i =0;i<d;i++){
            a[i] = (p[i]<< pow)&mask;
        }
        if(k>16){RefreshXOR(a,d);}else{RefreshXOR16(a,d);};
        
        if(k>16){SecAnd(p,a,tmp,d);}else{SecAnd16(p,a,tmp,d);};

        for(uint32_t i=0;i<d;i++){p[i] = tmp[i];}
    }
    
    pow = 1<<W;
    for(uint32_t i=0;i<d;i++){
        a[i] = (g[i]<<pow)&mask;
    }

    if(k>16){SecAnd(a,p,tmp,d);}else{SecAnd16(a,p,tmp,d);};
    for(uint32_t i=0;i<d;i++){
        g[i] = tmp[i] ^ g[i];
        z[i] = x[i] ^ y[i] ^ ((g[i]<<1)&mask);
    }
}

void SecA2BModpow2(uint32_t *a,
        uint32_t *out,
        uint32_t d,
        uint32_t k){
    
    if (d == 1){
        out[0] = a[0];
        return;
    }
    uint32_t lim = d/2;
    uint32_t y[d];
    uint32_t z[d];
    for(uint32_t i=0;i<d;i++){
        y[i] = 0; 
        z[i] = 0;
    }
    
    SecA2BModpow2(a,y,lim,k);
    if(k>16){RefreshXOR(y,d);}else{RefreshXOR16(y,d);;}
    
    SecA2BModpow2(&a[lim],&z[lim],d-lim,k);
    if(k>16){RefreshXOR(z,d);}else{RefreshXOR16(z,d);;}
    my_SecAdd(y,z,out,d,k);
}
