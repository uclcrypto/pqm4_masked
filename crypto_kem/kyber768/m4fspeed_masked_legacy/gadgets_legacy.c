#include "masked.h"
#include "gadgets_legacy.h"
#include "masked_utils.h"



void linear_arithmetic_refresh(Masked* x, unsigned q){
  int16_t r[2];
  for(int i=0; i < (KYBER_MASKING_ORDER-1); i+=2){
    rand_q(r);
    x->shares[i] = (x->shares[i] + r[0])%q;
    x->shares[KYBER_MASKING_ORDER] = (x->shares[KYBER_MASKING_ORDER] - r[0] + q)%q;

    x->shares[i+1] = (x->shares[i+1] + r[1])%q;
    x->shares[KYBER_MASKING_ORDER] = (x->shares[KYBER_MASKING_ORDER] - r[1] + q)%q;
  }
  #if KYBER_MASKING_ORDER%2 == 1
  rand_q(r);
  x->shares[KYBER_MASKING_ORDER-1] = (x->shares[KYBER_MASKING_ORDER-1] + r[0])%q;
  x->shares[KYBER_MASKING_ORDER] = (x->shares[KYBER_MASKING_ORDER] - r[0] + q)%q;
  #endif
}
void convert_B2A(Masked* x, Masked* y, unsigned k, unsigned q){
  Masked T[(1<<k)];
  Masked T_p[(1<<k)];
  
  for(int u=0; u < (1<<k); ++u){
    T[u].shares[0] = u%q;
    for(int i=1; i < KYBER_MASKING_ORDER+1; ++i) T[u].shares[i] = 0;
  }
  for(int i=0; i < KYBER_MASKING_ORDER; ++i){
    for(int u=0; u < (1<<k); ++u){
      for(int j=0; j < KYBER_MASKING_ORDER+1; ++j){
        T_p[u].shares[j] = T[u^(x->shares[i])].shares[j]; 
      }
    }
    for(int u=0; u < (1<<k); ++u){
      linear_arithmetic_refresh(&(T_p[u]), q); 
      for(int j=0; j < KYBER_MASKING_ORDER+1; ++j) T[u].shares[j] = T_p[u].shares[j];
    }
  }
  for(int i=0; i < KYBER_MASKING_ORDER+1; ++i) y->shares[i] = T[x->shares[KYBER_MASKING_ORDER]].shares[i]; 
  linear_arithmetic_refresh(y, q);
}


void convert_2_l_to_1bit_bool(Masked* x, Masked* b, unsigned l){
  if (l == 4){
    uint16_t T[KYBER_MASKING_ORDER+1];
    uint16_t r;
    T[0] = 0x0FF0;

    for(int i=1; i < KYBER_MASKING_ORDER+1; ++i) T[i] = 0;
    
    for(int i=0; i < KYBER_MASKING_ORDER; ++i){
      for(int j=0; j < KYBER_MASKING_ORDER+1; ++j) T[j] = (T[j] << (16-x->shares[i])) + (T[j]>>(x->shares[i]));
      for(int j=0; j < KYBER_MASKING_ORDER; ++j){
        r = rand16();
        T[j] ^= r;
        T[KYBER_MASKING_ORDER] ^= r;
      }
    }
    for(int i=0; i < KYBER_MASKING_ORDER+1; ++i) b->shares[i] = (T[i]>>(x->shares[KYBER_MASKING_ORDER]))&1; 
    for(int j=0; j < KYBER_MASKING_ORDER; ++j){
      r = rand16();
      b->shares[j] ^= r;
      b->shares[KYBER_MASKING_ORDER] ^= r;
    }
  }
  else if (l == 5){
    uint32_t T[KYBER_MASKING_ORDER+1];
    uint32_t r;


    T[0] = 0x00FFFF00;

    for(int i=1; i < KYBER_MASKING_ORDER+1; ++i) T[i] = 0;
    
    for(int i=0; i < KYBER_MASKING_ORDER; ++i){
      for(int j=0; j < KYBER_MASKING_ORDER+1; ++j) T[j] = (T[j] << (32-x->shares[i])) + (T[j]>>(x->shares[i]));
      for(int j=0; j < KYBER_MASKING_ORDER; ++j){
        r = rand32();
        T[j] ^= r;
        T[KYBER_MASKING_ORDER] ^= r;
      }
    }
    for(int i=0; i < KYBER_MASKING_ORDER+1; ++i) b->shares[i] = (T[i]>>(x->shares[KYBER_MASKING_ORDER]))&1; 
    for(int j=0; j < KYBER_MASKING_ORDER; ++j){
      r = rand32();
      b->shares[j] ^= r;
      b->shares[KYBER_MASKING_ORDER] ^= r;
    }
  }

  else if (l == 6){
    uint64_t T[KYBER_MASKING_ORDER+1];
    uint64_t r;

    T[0] = 0x0000FFFFFFFF0000LLU;

    for(int i=1; i < KYBER_MASKING_ORDER+1; ++i) T[i] = 0;
    
    for(int i=0; i < KYBER_MASKING_ORDER; ++i){
      for(int j=0; j < KYBER_MASKING_ORDER+1; ++j) T[j] = (T[j] << (64-x->shares[i])) + (T[j]>>(x->shares[i]));
      for(int j=0; j < KYBER_MASKING_ORDER; ++j){
        r = rand64();
        T[j] ^= r;
        T[KYBER_MASKING_ORDER] ^= r;
      }
    }
    for(int i=0; i < KYBER_MASKING_ORDER+1; ++i) b->shares[i] = (T[i]>>(x->shares[KYBER_MASKING_ORDER]))&1; 
    for(int j=0; j < KYBER_MASKING_ORDER; ++j){
      r = rand64();
      b->shares[j] ^= r;
      b->shares[KYBER_MASKING_ORDER] ^= r;
    }
  }

  else if (l == 7){
    uint64_t T1[KYBER_MASKING_ORDER+1];
    uint64_t T2[KYBER_MASKING_ORDER+1];
    uint64_t r;
    unsigned shift;

  

    T1[0] = 0xFFFFFFFF00000000LLU;
    T2[0] = 0x00000000FFFFFFFFLLU;
 

    for(int i=1; i < KYBER_MASKING_ORDER+1; ++i) {
      T1[i] = 0;
      T2[i] = 0;
    }
    
    for(int i=0; i < KYBER_MASKING_ORDER; ++i){
      for(int j=0; j < KYBER_MASKING_ORDER+1; ++j){
        shift = x->shares[i];
        if (shift%64 != 0){
          r = T1[j];
          T1[j] = (T1[j] >> (shift%64)) + (T2[j] << (64-(shift%64)));
          T2[j] = (T2[j] >> (shift%64)) + (r << (64-(shift%64)));
        }
        if (shift >= 64){
          r = T2[j];
          T2[j] = T1[j];
          T1[j] = r;
         }
      }
      for(int j=0; j < KYBER_MASKING_ORDER; ++j){
        r = rand64();
        T1[j] ^= r;
        T1[KYBER_MASKING_ORDER] ^= r;
        r = rand64();
        T2[j] ^= r;
        T2[KYBER_MASKING_ORDER] ^= r;
      }
    }
    

    
    for(int i=0; i < KYBER_MASKING_ORDER+1; ++i){
      shift = x->shares[KYBER_MASKING_ORDER];
      if (shift < 64) b->shares[i] = (T1[i]>>(shift))&1;
      else            b->shares[i] = (T2[i]>>(shift-64))&1;
    }
    for(int j=0; j < KYBER_MASKING_ORDER; ++j){
      r = rand32();
      b->shares[j] ^= r;
      b->shares[KYBER_MASKING_ORDER] ^= r;
    }

  }

  else if (l == 8){

    uint64_t T[4][KYBER_MASKING_ORDER+1];
    uint64_t T_p[4][KYBER_MASKING_ORDER+1];
    uint64_t r;
    unsigned shift, jump, small_shift;

    T[0][0] = 0x0000000000000000LLU;
    T[1][0] = 0xFFFFFFFFFFFFFFFFLLU;
    T[2][0] = 0xFFFFFFFFFFFFFFFFLLU;
    T[3][0] = 0x0000000000000000LLU;

    for(int i=1; i < KYBER_MASKING_ORDER+1; ++i) {
      T[0][i] = 0;
      T[1][i] = 0;
      T[2][i] = 0;
      T[3][i] = 0;
    }
    
    for(int i=0; i < KYBER_MASKING_ORDER; ++i){
      for(int j=0; j < KYBER_MASKING_ORDER+1; ++j){
        shift = x->shares[i];
        jump = shift/64;
        small_shift = shift%64;

        for(int k=0; k < 4; ++k)
          T_p[(k-jump+4)%4][j] = T[k][j];
        for(int k=0; k < 4; ++k){
          if(small_shift != 0) 
            T[k][j] = (T_p[k][j] >> small_shift) | (T_p[(k+1)%4][j] << (64-small_shift));
          else
            T[k][j] = T_p[k][j];
        }
        
      }
      for(int j=0; j < KYBER_MASKING_ORDER; ++j){
        r = rand64();
        T[0][j] ^= r;
        T[0][KYBER_MASKING_ORDER] ^= r;
        r = rand64();
        T[1][j] ^= r;
        T[1][KYBER_MASKING_ORDER] ^= r;
        r = rand64();
        T[2][j] ^= r;
        T[2][KYBER_MASKING_ORDER] ^= r;
        r = rand64();
        T[3][j] ^= r;
        T[3][KYBER_MASKING_ORDER] ^= r;

      }
    }
    
    for(int i=0; i < KYBER_MASKING_ORDER+1; ++i){
      shift = x->shares[KYBER_MASKING_ORDER];
      b->shares[i] = (T[shift/64][i] >> (shift%64))&1;
    }

    for(int j=0; j < KYBER_MASKING_ORDER; ++j){
      r = rand32();
      b->shares[j] ^= r;
      b->shares[KYBER_MASKING_ORDER] ^= r;
    }
  }
} 


static unsigned switch_table[9] =  {6, 7, 7, 7, 8, 8, 8, 8, 8}; // Value of \ell in the paper

void modulus_switch(Masked* x, unsigned q, unsigned shift){
  /* 
   * Modulus switch between Z_q and Z_{2^shift} 
   * round((x<<shift)/q) = ((x<<(shift+1) + q1)//(2*q)
   * No overflow should appear for the values we use in the paper
   */
  int64_t temp;
  for(int i =0; i < KYBER_MASKING_ORDER+1; ++i) {
    temp = (int64_t)(x->shares[i]) << (shift+1);
    temp = (temp+q)/(2*q);
    x->shares[i] = (int)temp&((1<<shift)-1);
  }
}

void kyber_decryption(Masked* x, Masked* b){
  unsigned l = switch_table[KYBER_MASKING_ORDER-1];
  modulus_switch(x, KYBER_Q, l);
  convert_2_l_to_1bit_bool(x, b, l);
}
void boolean_refresh(Masked* x, unsigned k){
    int r;
    for(int i=0; i< KYBER_MASKING_ORDER+1; ++i){
        for(int j=i+1; j < KYBER_MASKING_ORDER+1; ++j){
            r = (int) rand32() & ((1<<k)-1);
            x->shares[i] = (x->shares[i] ^ r);
            x->shares[j] = (x->shares[j] ^ r);
        }
    }
}
