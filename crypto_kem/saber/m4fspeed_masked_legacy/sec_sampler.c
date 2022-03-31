#include <stdint.h>
#include <stddef.h>
#include "bs_gadgets.h"
#include "gadgets_legacy.h"
#include "bench.h"
#include "sec_sampler.h"
#include "masked_representations.h"
#include "masked_utils.h"
#define KYBER_MASKING_ORDER (D-1)
#define D NSHARES
#define GET_BIT_SHARING(i,j) ((i)*(j)) 

void map_bitslice_masked(const uint32_t *x, // [w*D]
        uint32_t *bitslice, // [k*D]
        uint32_t k, // 13 -> 20
        uint32_t w, // 32
        uint32_t d){

    for(uint32_t j=0;j<k*d;j++){
        bitslice[j] = 0;
    }
    
    if((w & 0x3) == 0){
        for(uint32_t s=0;s<d;s++){
            for(uint32_t j=0; j<w;j+=4){
                uint32_t xd = x[s+j*d];
                uint32_t xd2 = x[s+(j+1)*d];
                uint32_t xd3 = x[s+(j+2)*d];
                uint32_t xd4 = x[s+(j+3)*d];
                for(uint32_t i=0; i<k;i++){
                    bitslice[s + i*d] |= (xd &0x1)<<j;
                    bitslice[s + i*d] |= (xd2 &0x1)<<(j+1);
                    bitslice[s + i*d] |= (xd3 &0x1)<<(j+2);
                    bitslice[s + i*d] |= (xd4 &0x1)<<(j+3);
                    xd = xd >> 1;
                    xd2 = xd2 >> 1;
                    xd3 = xd3 >> 1;
                    xd4 = xd4 >> 1;
                }
            }
        }
    }else{
        for(uint32_t s=0;s<d;s++){
            for(uint32_t j=0; j<w;j+=1){
                uint32_t xd = x[s+j*d];
                for(uint32_t i=0; i<k;i++){
                    bitslice[s + i*d] |= (xd &0x1)<<j;
                    xd = xd >> 1;
                }
            }
        }

    }
}

void unmap_bitslice_masked(const uint32_t *bitslice,
        uint32_t *x,
        uint32_t k,
        uint32_t w,
        uint32_t d){
    uint32_t tmp[k*d];
    uint32_t x0,x1,y;
    for(uint32_t j=0; j<k*d;j++){
        tmp[j] = bitslice[j];
    }
    for(uint32_t j=0; j<w;j+=2){
        for(uint32_t s=0;s<d;s++){
            x0=0;x1=0;
            for(uint32_t i=0; i<k;i++){
                y = tmp[s+i*d];
                x0 |= ((y>>j) & 0x1)<<i;
                x1 |= ((y>>(j+1)) & 0x1)<<i;
            }
            x[s+j*d] = x0;
            x[s+(j+1)*d] = x1;
        }
    }

}

void map_bitslice(uint32_t *x,
        uint32_t *bitslice,
        uint32_t k,
        uint32_t w){
    uint32_t tmp[w];
    for(uint32_t j=0;j<w;j++){
        tmp[j] = x[j];
    }
    for(uint32_t i=0; i<k;i++){
        bitslice[i] = 0;
        for(uint32_t j=0; j<w;j++){
            bitslice[i] |= ((tmp[j]) &0x1)<<j;
            tmp[j] = tmp[j]>>1;
        }
    }
}

void unmap_bitslice(uint32_t *bitslice,
        uint32_t *x,
        uint32_t k,
        uint32_t w){

    for(uint32_t j=0; j<w;j++){
        x[j] = 0;
        for(uint32_t i=0; i<k;i++){
            x[j] |= ((bitslice[i]>>j) & 0x1)<<i;
        }
    }
}

/*
 * output: z[D*kappa] sum(z) = HW(x)
 * input: x[D*lambda]
 * 
 * x is shared on kappa bits
 * z is shared on lambda bits
 * See algo 11
 */

void secbitadd(uint32_t *z, const uint32_t *x, const size_t lambda, const size_t kappa){
    uint32_t t[D*lambda],w[D];
    size_t i,j,l;

    for(i=0;i<D*lambda;i++){
        t[i] = 0;
        z[i] = 0;
    }

    for(j=0;j<kappa;j++){
        SecXORbs(&t[GET_BIT_SHARING(0,D)],
                    &z[GET_BIT_SHARING(0,D)],
                    &x[GET_BIT_SHARING(j,D)]);

        for(i=0;i<D;i++){
            w[i] = (&(x[GET_BIT_SHARING(j,D)]))[i]; 
        }
        for(l=1;l<lambda;l++){
            SecANDbs(w,
                    w,
                    &z[GET_BIT_SHARING(l-1,D)]);
            SecXORbs(&t[GET_BIT_SHARING(l,D)],
                    &z[GET_BIT_SHARING(l,D)],
                    w);
        }
        for(i=0;i<(D*lambda);i++){
            z[i] = t[i];
        }
    }
}


/// algo 12
// z = z - HW(x)
void secbitsub(uint32_t *z, 
        const uint32_t *x, 
        const size_t lambda, const size_t kappa){
    uint32_t t[D*lambda],w[D],u[D];
    size_t i,j,l;
    for(i=0;i<D*lambda;i++){
        t[i] = 0;
    }

    for(j=0;j<kappa;j++){
        SecXORbs(&t[GET_BIT_SHARING(0,D)],
                    &z[GET_BIT_SHARING(0,D)],
                    &x[GET_BIT_SHARING(j,D)]);

        for(i=0;i<D;i++){
            w[i] = (&(x[GET_BIT_SHARING(j,D)]))[i]; 
        }
        for(l=1;l<lambda;l++){
            for(i=0;i<D;i++){
                u[i] = z[GET_BIT_SHARING(l-1,D)+i];
            }
            u[0] ^= 0xffffffff;
            SecANDbs(w,
                    w,
                    u);
            SecXORbs(&t[GET_BIT_SHARING(l,D)],
                    &z[GET_BIT_SHARING(l,D)],
                    w);
        }
        for(i=0;i<(D*lambda);i++){
            z[i] = t[i];
        }
    }
}

void secconstadd2(uint32_t *z){
    SecXORbs(&z[GET_BIT_SHARING(2,D)],
                &z[GET_BIT_SHARING(2,D)],
                &z[GET_BIT_SHARING(1,D)]);
    z[GET_BIT_SHARING(1,D)] ^= 0xffffffff;
}

void secconstadd4(uint32_t *z){
    SecXORbs(&z[GET_BIT_SHARING(3,D)],
                &z[GET_BIT_SHARING(3,D)],
                &z[GET_BIT_SHARING(2,D)]);
    z[GET_BIT_SHARING(2,D)] ^= 0xffffffff;
}
// algo 9
#define SPOG_GADGETS
void sec_b2a_q(uint32_t *a, uint32_t *x,size_t k, uint32_t kb){
    uint32_t b[D],tmp[D];
    size_t i,j;
    
    uint32_t mask = (1<<kb)-1;
    uint32_t Q = 1<<kb;

    for(i=0;i<D;i++){
        tmp[i] = (x[i] >> (k-1))&0x1;
    }
#ifdef SPOG_GADGETS
    sec_b2a_qbit(a,tmp,kb);
#else
    convert_B2A_ui32(tmp,a,1,kb);
#endif
    for(j=1;j<k;j++){
        for(i=0;i<D;i++){
            tmp[i] = (x[i] >> (k-j-1))&0x1;
        }
#ifdef SPOG_GADGETS
        sec_b2a_qbit(b,tmp,kb);
#else
        convert_B2A_ui32(tmp,b,1,kb);
#endif
        for(i=0;i<D;i++){
            a[i] = (2 * a[i] + b[i])&mask;
        }
    }
}
  
void sec_sampler2(uint32_t *a, const uint32_t *x, 
                        const uint32_t *y,
                        const uint32_t kappa,const uint32_t kmodulus){

    uint32_t ta,tb,tc,td,te,tf,tg;
    uint32_t mask = (1<<kmodulus)-1;

    int lambda = 4;
    uint32_t tmp[lambda*D];
    uint32_t tmpb[kmodulus*D];
    uint32_t tmpa[kmodulus*D];

    for(int i=0;i<D;i++){
        tmp[GET_BIT_SHARING(lambda-1,D)+i] = 0;
    }

    start_bench(cbd_bool);
    secbitadd(tmp,x,lambda-1,kappa);
    secbitsub(tmp,y,lambda,kappa);
    secconstadd4(tmp);
    secconstadd4(tmp);

    uint32_t tmp_nbs[32*D];

    // deserialize bitslice
    unmap_bitslice_masked(tmp,tmp_nbs,lambda,32,D);

    stop_bench(cbd_bool);

    // convert bool to 
    start_bench(cbd_b2a);
    for(size_t i=0;i<32;i++){
        sec_b2a_q(&a[GET_BIT_SHARING(i,D)],
                &tmp_nbs[GET_BIT_SHARING(i,D)],
                lambda,kmodulus);

        uint32_t lol = a[GET_BIT_SHARING(i,D)];
        a[GET_BIT_SHARING(i,D)] = (lol + (1<<kmodulus) - 8)&mask;
    }

    stop_bench(cbd_b2a);
}

// algo 7 in SOPG 
void sec_b2a_qbit_n(uint32_t *c,const uint32_t *a,const uint32_t xn,const uint32_t n, const uint32_t kb){
    uint32_t b[n];
    uint32_t r;
    uint32_t Q = 1<<kb;
    uint32_t mask =  (1<<kb)-1;
    size_t i,j;

    b[n-1] = rand16()&mask;
    b[0] = (a[0] - b[n-1] + Q)&mask;

    for(j=1;j<n-1;j++){
        r = rand16()&mask;
        b[j] = (a[j] - r + Q)&mask;
        b[n-1] = (b[n-1] + r)&mask;
    }

    for(i=0;i<n;i++){
        c[i] = b[i] + 2*Q;
        c[i] -= (2*b[i]*xn);
        c[i] = c[i]&mask;
    }
    c[0] = (c[0] + xn)&mask;
}



// algo 6
void b2a_qbit(uint32_t *a,uint32_t *x, const uint32_t kb){
    a[0] = x[0];

    int id =0;
    for(int i=1;i<D;i++){
        sec_b2a_qbit_n(a,a,x[i],i+1,kb);
    }
}

// algo 5
void sec_b2a_qbit(uint32_t *a, uint32_t *x,uint32_t kb){
    b2a_qbit(a,x,kb);
    refresh_add(a,kb);
}

// algo 8
void refresh_add(uint32_t *a,uint32_t kb){
    size_t i,j;
    uint32_t rnd;
    uint32_t mask = (1<<kb)-1;
    uint32_t Q = 1<<kb;
    for(i=0;i<D-1;i++){
        for(j=i+1;j<D;j++){
            rnd = rand16() & mask;
            a[i] = (a[i] + rnd);
            a[j] = (a[j] + Q - rnd) & mask;
        }
    }
    for(i=0;i<D;i++){
        a[i] = a[i] & mask;
    }
}


