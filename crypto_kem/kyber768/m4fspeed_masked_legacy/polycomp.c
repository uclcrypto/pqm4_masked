#include "polycomp.h"
#include "gadgets_legacy.h"
#include "masked_utils.h"
#include "params.h"
#include <stdint.h>

static unsigned umodulus_switch(unsigned x, unsigned q_start, unsigned q_end){
  return (2*q_end*x+q_start)/(2*q_start);
}

static unsigned compress(unsigned x, unsigned q, unsigned d){
  return umodulus_switch(x, q, 1<<d)%(1<<d);
}

static unsigned decompress(unsigned x, unsigned q, unsigned d){
  return umodulus_switch(x, 1<<d, q);
}



static void sec_and(Masked* x, Masked* y, Masked* res, int k){

#if KYBER_MASKING_ORDER == 1
    uint16_t u = rand16()&((1<<k)-1);
    uint16_t z;
    z = u ^ (x->shares[0] & y->shares[0]);
    z = z ^ (x->shares[0] & y->shares[1]);
    z = z ^ (x->shares[1] & y->shares[0]);
    z = z ^ (x->shares[1] & y->shares[1]);
    res->shares[0] = z;
    res->shares[1] = u;

#else
    Masked r;
    uint16_t i, j, z_ij, z_ji;
    for(i=0; i < KYBER_MASKING_ORDER + 1; ++i) r.shares[i] = x->shares[i] & y->shares[i];
    for(i=0; i < KYBER_MASKING_ORDER + 1; ++i)
        for(j=i+1; j < KYBER_MASKING_ORDER + 1; ++j){
            z_ij  = rand16()&((1<<k)-1);
            z_ji  = (x->shares[i] & y->shares[j]) ^ z_ij;
            z_ji ^= (x->shares[j] & y->shares[i]);
            r.shares[i] ^= z_ij;
            r.shares[j] ^= z_ji;            
        }
    for(i=0; i < KYBER_MASKING_ORDER + 1; ++i) res->shares[i] = r.shares[i];
#endif

}

static void sec_mult(Masked* a, Masked* b, Masked* c, unsigned q){
  uint16_t r[2], t;
  int pool_size = (KYBER_MASKING_ORDER*(KYBER_MASKING_ORDER+1))/2;
  uint16_t pool[pool_size];
  for(int i=0; i < pool_size; i += 2){
    rand_q(r);
    pool[i  ] = r[0];
    pool[i+1] = r[1];
  }
  if ((pool_size%2) == 1){
    rand_q(r);
    pool[pool_size-1] = r[0];
  }

  int cpt = 0;

  for(int i=0; i < KYBER_MASKING_ORDER+1; ++i) c->shares[i] = (a->shares[i]*b->shares[i])%q;

  for(int i=0; i < KYBER_MASKING_ORDER+1; ++i){
    for(int j=i+1; j < KYBER_MASKING_ORDER+1; ++j){
      t = ((pool[cpt]+a->shares[i]*b->shares[j])+a->shares[j]*b->shares[i])%q;
      c->shares[i] = (c->shares[i] + q - pool[cpt])%q;
      c->shares[j] = (c->shares[j] + t)%q;
      cpt++;
    }
  }
}



static uint32_t genrand(int l)
{
  if (l==32) return rand32();
  return rand32() & ((1 << l)-1);
}



static void Expand(uint32_t *x,uint32_t *xp,int k,int n2,int n)
{
  for(int i=0;i<n/2;i++)
  {
    uint32_t r=genrand(k);
    xp[2*i]=x[i] ^ r;
    xp[2*i+1]=r;
  }
  if ((n & 1)==1) 
  {
    if (n2==n/2)
      xp[n-1]=0;
    else
      xp[n-1]=x[n2-1];
  }
}

static void SecAnd(uint32_t *a,uint32_t *b,uint32_t *c,int k,int n)
{
  for(int i=0;i<n;i++)
    c[i]=a[i] & b[i];

  for(int i=0;i<n;i++)
  {
    for(int j=i+1;j<n;j++)
    {
      uint32_t tmp=rand32(); //rand();
      uint32_t tmp2=(tmp ^ (a[i] & b[j])) ^ (a[j] & b[i]);
      c[i]^=tmp;
      c[j]^=tmp2;
    }
  }
  for(int i=0;i<n;i++) c[i]=c[i] % (1 << k);
}

static void SecAdd(uint32_t *x,uint32_t *y,uint32_t *z,int k,int n)
{
  uint32_t u[n];
  for(int i=0;i<n;i++) u[i]=0;
  uint32_t w[n];
  SecAnd(x,y,w,k,n);
  uint32_t a[n];
  for(int i=0;i<n;i++) a[i]=x[i] ^ y[i];
  for(int j=0;j<k-1;j++)
  {
    uint32_t ua[n];
    SecAnd(u,a,ua,k,n);
    for(int i=0;i<n;i++) u[i]=(2*(ua[i] ^ w[i])) % (1 << k);
  }
  for(int i=0;i<n;i++) z[i]=x[i] ^ y[i] ^ u[i];
}


uint32_t GoubinAB(uint32_t A,uint32_t r,int k)
{
  uint32_t G=rand32();
  uint32_t T=G << 1;
  uint32_t x=G ^ r;
  uint32_t O=G & x;
  x=T ^ A;
  G=G ^ x;
  G=G & r;
  O=O ^ G;
  G=T & A;
  O=O ^ G;
  for(int i=1;i<k;i++)
  {
    G=T & r;
    G=G ^ O;
    T=T & A;
    G=G ^ T;
    T=G << 1;
  }
  x=x ^ T;
  return x;
}


void ConvertAB(uint32_t *A,uint32_t *z,int k,int n)
{
  if(n==1)
  {
    z[0]=A[0];
    return;
  }

  if(n==2)
  {
    z[0]=GoubinAB(A[0],A[1],k);
    z[1]=A[1];
    return;
  }

  uint32_t x[n/2];
  ConvertAB(A,x,k,n/2);
  uint32_t xp[n];
  Expand(x,xp,k,n/2,n);
  
  uint32_t y[(n+1)/2];
  ConvertAB(A+n/2,y,k,(n+1)/2);
  uint32_t yp[n];
  Expand(y,yp,k,(n+1)/2,n);

  SecAdd(xp,yp,z,k,n);
}


static void compress_inv(unsigned c, unsigned range[2], unsigned q, unsigned d){
  unsigned B, y, a, b;

  B = ((q+(1<<d))>>(d+1));
  y = decompress(c, q, d);
  a = (q + y - B) % q;
  b = (y + B)     % q;
  if (compress(a, q, d) != c) a = (a+1)%q;
  if (compress(b, q, d) != c) b = (b+q-1)%q;
  range[0] = a; range[1] = b;
}


void range_compare(Masked* x, Masked* z, unsigned c, unsigned k, unsigned q){
  Masked temp;
  unsigned range[2], delta;
  compress_inv(c, range, q, k);
  
  z->shares[0] = 1;
  for(int i=1; i < KYBER_MASKING_ORDER+1; ++i) z->shares[i] = 0;
  x->shares[0] = (x->shares[0] + q - range[0])%q;

  if (range[1] < range[0]) range[1] += q;
  delta = range[1] - range[0] + 1;
  for(unsigned i=0; i < delta; ++i){
    sec_mult(x, z, &temp, q);
    for(int j=0; j < KYBER_MASKING_ORDER+1; ++j) {
      z->shares[j] = temp.shares[j];
    }
    x->shares[0] = (x->shares[0] + q - 1)%q;
  }
}

void convert_A2B_CGV14_32bits(u32Masked* x, u32Masked* y, unsigned k){
  ConvertAB(x->shares, y->shares, k, KYBER_MASKING_ORDER+1);
}

int zero_test_mod_mult(Masked* a, int q){
  uint16_t u_j, B;
  for(int j = 0; j < KYBER_MASKING_ORDER+1; ++j){
      // FIXME use uniform randomness %Q
    u_j = 1 + rand16()%(q-1);
    for(int i=0; i < KYBER_MASKING_ORDER+1; ++i) a->shares[i] = (u_j*a->shares[i])%q;
    linear_arithmetic_refresh(a, q);
  }
  B = 0;
  for(int i=0; i < KYBER_MASKING_ORDER+1; ++i) B = (B + a->shares[i])%q;
  if (B == 0) return 1;
  else        return 0;
}

int zero_testing_prime_multi(Masked* ppoly, int q, const int SIZE){
  Masked a, t, C;
  for(int i=0; i < KYBER_MASKING_ORDER+1; ++i) a.shares[i] = rand16()%q;
  sec_mult(&a, &ppoly[0], &C, q);

  for(int i=1; i < SIZE; ++i){
    for(int j=0; j < KYBER_MASKING_ORDER+1; ++j) a.shares[j] = rand16()%q;
    sec_mult(&a, &ppoly[i], &t, q);
    for(int j=0; j < KYBER_MASKING_ORDER+1; ++j) C.shares[j] = (C.shares[j] + t.shares[j])%q;
  }
  return zero_test_mod_mult(&C, q);
}

int zero_test_poly_mul(Masked* ppoly, int q, int lambda, const int SIZE){
  int b=1;
  for(int i=0; (i < lambda); ++i){
    b &= zero_testing_prime_multi(ppoly, q, SIZE);
  }
  return b;
}

int zero_test_poly_mul_with_reduction(Masked* ppoly, int q, int kappa, const int SIZE){
  uint16_t a;
  Masked y[kappa];
  
  for(int k=0; k < kappa; ++k){
    for(int i=0; i < KYBER_MASKING_ORDER+1; ++i) y[k].shares[i] = 0;
    for(int j=0; j < SIZE; j+=1){
      a = sampleq();
      for(int i=0; i < KYBER_MASKING_ORDER+1; ++i) y[k].shares[i] = (y[k].shares[i] + a*ppoly[j].shares[i])%q; 
    }
  }

  return zero_test_poly_mul(y, q, kappa, kappa);
}


void boolean_zero_test(Masked* x, Masked* y, int k, int logk){
  Masked z;
  y->shares[0] = (~x->shares[0]) | ((1<<(1<<logk))-(1<<k));
  for(int i=1; i < KYBER_MASKING_ORDER+1; ++i) y->shares[i] = x->shares[i];

  for(int i=0; i < logk; ++i){
    for(int j=0; j < KYBER_MASKING_ORDER+1; ++j) z.shares[j] = y->shares[j] >> (1 << i);
    boolean_refresh(&z, k);
    sec_and(y, &z, y, k); 
  }
  for(int i=0; i < KYBER_MASKING_ORDER+1; ++i) y->shares[i] = (y->shares[i])&1;

  boolean_refresh(y, k);

}


void bool_poly_zero_test(Masked* ppoly, Masked* b, int k, int logk, const int SIZE){
  Masked temp1, temp2;
  Masked *p1=&temp1, *p2=&temp2, *swap;
  temp1.shares[0] = 0xFFFF; 
  for(int i=1; i < KYBER_MASKING_ORDER+1; ++i) temp1.shares[i] = 0;

  for(int i=0; i < SIZE; ++i){
    ppoly[i].shares[0] = ~ppoly[i].shares[0]; 
    sec_and(&(ppoly[i]), p1, p2, k);
    swap = p2;
    p2 = p1;
    p1 = swap;
  }
  temp1.shares[0] = ~temp1.shares[0];

  boolean_zero_test(&temp1, b, k, logk);

}


void high_order_compress(Masked* x, Masked* y, unsigned q, unsigned k, unsigned ell){
  u32Masked z;
  uint64_t temp;
  
  temp = (uint64_t)(x->shares[0]);
  z.shares[0] = (temp<<(k+ell+1))/(2*q);
  z.shares[0] = (z.shares[0] + (1<<(ell-1)))%(1<<(k+ell));

  for(int i=1; i < KYBER_MASKING_ORDER+1; ++i) {
    temp = (uint64_t)(x->shares[i]);
    z.shares[i] = ((temp<<(k+ell+1))/(2*q))%(1<<(k+ell));
  }

  u32Masked t;
  convert_A2B_CGV14_32bits(&z, &t, k+ell);
  for(int i=0; i < KYBER_MASKING_ORDER+1; ++i) y->shares[i] = (t.shares[i] >> ell)%(1<<(k));
} 

int kyber_poly_comp_hybrid(Masked* mmasked_poly, uint16_t* ppoly){
  int l1=KYBER_K*KYBER_N, l2=KYBER_N;
  Masked z[l1+1];
  Masked b2;
  unsigned d1=KYBER_DU, d2=KYBER_DV, alpha=16, q=KYBER_Q;


    alpha=0;
    uint32_t prod = KYBER_Q*D;
    while(prod > 0){
        prod = prod >> 1;
        alpha++;
    }
  for(int i=0; i < l1; ++i) range_compare(&(mmasked_poly[i]), &(z[i]), ppoly[i], d1, q);
  for(int i=0; i < l2; ++i) high_order_compress(&(mmasked_poly[l1 + i]), &(mmasked_poly[l1 + i]), q, d2, alpha);

  for(int i=0; i < l2; ++i) mmasked_poly[l1 + i].shares[0] = (mmasked_poly[l1 + i].shares[0] ^ ppoly[l1 + i])&((1<<d2)-1); 


  bool_poly_zero_test(mmasked_poly+l1, &b2, d2, 2, l2);
  b2.shares[0] = (~b2.shares[0])&1;
  for(int i=1; i < KYBER_MASKING_ORDER+1; ++i) b2.shares[i] &=1;
  convert_B2A(&b2, z+l1, 1, q);
  return zero_test_poly_mul_with_reduction(z, q, 11, l1+1);

}
