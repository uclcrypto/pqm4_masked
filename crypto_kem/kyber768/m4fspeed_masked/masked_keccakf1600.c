/* Based on the implementation "libkeccak-tiny" by David Leon Gil.
 * available at https://github.com/coruus/keccak-tiny under CC0 License.
 * */

#include <stdint.h>
#include <assert.h>
#include "masked.h"
#include "masked_keccakf1600.h"
#include "gadgets.h"

#define NROUNDS 24


void MaskedKeccakF1600_StateExtractBytes(
        const MaskedKeccakState state,
        unsigned char *data,
        unsigned int offset,
        unsigned int length,
        size_t data_msk_stride,
        size_t data_data_stride
        )
{
    for(size_t i=0;i<length;i++) {
        for(size_t j=0; j<NSHARES; j++) {
            data[i*data_data_stride+j*data_msk_stride] =
                state[j][(offset + i) >> 3] >> (8*((offset + i) & 0x07));
        }
    }
}

void MaskedKeccakF1600_StateXORBytes(
        MaskedKeccakState state,
        const unsigned char *data,
        unsigned int offset,
        unsigned int length,
        size_t data_msk_stride,
        size_t data_data_stride
        )
{
    for(size_t i = 0; i < length; i++) {
        for(size_t j=0; j<NSHARES; j++) {
            state[j][(offset + i) >> 3] ^=
                (uint64_t)data[i*data_data_stride+j*data_msk_stride] << (8 * ((offset + i) & 0x07));
        }
    }
}

void MaskedKeccakF1600_StateXORPublicBytes(
        MaskedKeccakState state,
        const unsigned char *data,
        unsigned int offset,
        unsigned int length
        ) {
    for(size_t i = 0; i < length; i++) {
        state[0][(offset + i) >> 3] ^= (uint64_t)data[i] << (8 * ((offset + i) & 0x07));
    }
}

/******** The Keccak-f[1600] permutation ********/

/*** Constants. ***/
static const uint8_t rho[24] = \
  { 1,  3,   6, 10, 15, 21,
    28, 36, 45, 55,  2, 14,
    27, 41, 56,  8, 25, 43,
    62, 18, 39, 61, 20, 44};
static const uint8_t pi[24] = \
  {10,  7, 11, 17, 18, 3,
    5, 16,  8, 21, 24, 4,
   15, 23, 19, 13, 12, 2,
   20, 14, 22,  9, 6,  1};
static const uint64_t RC[24] = \
  {1ULL, 0x8082ULL, 0x800000000000808aULL, 0x8000000080008000ULL,
   0x808bULL, 0x80000001ULL, 0x8000000080008081ULL, 0x8000000000008009ULL,
   0x8aULL, 0x88ULL, 0x80008009ULL, 0x8000000aULL,
   0x8000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL, 0x8000000000008003ULL,
   0x8000000000008002ULL, 0x8000000000000080ULL, 0x800aULL, 0x800000008000000aULL,
   0x8000000080008081ULL, 0x8000000000008080ULL, 0x80000001ULL, 0x8000000080008008ULL};

/*** Helper macros to unroll the permutation. ***/
#define ROL(a, offset) ((a << offset) ^ (a >> (64-offset)))
#define REPEAT6(e) e e e e e e
#define REPEAT24(e) REPEAT6(e e e e)
#define REPEAT5(e) e e e e e
#define FOR5(v, s, e) \
  v = 0;            \
  REPEAT5(e; v += s;)

static void unmask_state(KeccakState st, const MaskedKeccakState msk_st) {
    for (size_t i=0; i<KECCAK_NWORDS; i++) {
        st[i] = 0;
        for (size_t j=0; j<NSHARES; j++) {
            st[i] ^= msk_st[j][i];
        }
    }
}
void disp_keccak_state(uint64_t *st) {
    (void)(st); // for warning
#if 0
    char buf[120];
    for (size_t j=0; j<KECCAK_NWORDS; j+=5) {
        sprintf(buf, "\t%llx %llx %llx %llx %llx", st[j+0], st[j+1], st[j+2], st[j+3], st[j+4]);
        hal_send_str(buf);
    }
#endif
}

void MaskedKeccakF1600_StatePermute(MaskedKeccakState state) {
    KeccakState umsk;
    unmask_state(umsk, state);
    disp_keccak_state(umsk);
  uint8_t x, y;

  for (int i = 0; i < NROUNDS; i++) {
      // Sharewise implementation for Theta, Rho and phi
      for (int j=0; j<NSHARES; j++) {
        uint64_t* a = &state[j][0];
        uint64_t b[5];
        uint64_t t = 0;
        // Theta
        FOR5(x, 1,
             b[x] = 0;
             FOR5(y, 5,
                  b[x] ^= a[x + y]; ))
        FOR5(x, 1,
             FOR5(y, 5,
                  a[y + x] ^= b[(x + 4) % 5] ^ ROL(b[(x + 1) % 5], 1); ))
        // Rho and pi
        t = a[1];
        x = 0;
        REPEAT24(b[0] = a[pi[x]];
                 a[pi[x]] = ROL(t, rho[x]);
                 t = b[0];
                 x++; )
      }
      // Chi: non-linear -> not sharewise.
      // Masked gadgets are implemented on 32-bit words and Chi does not contain rotations,
      // so we can work on 32-bit words
      for (y=0; y<25; y+=5)
      for (int off=0; off<2; off++)
      {
          uint32_t sb_state[5*NSHARES];
          size_t sb_state_msk_stride = 1; // in 32-bit words
          size_t sb_state_data_stride = NSHARES; // in 32-bit words
          uint32_t *sb_in = ((uint32_t *) state)+ 2*y+off;
          size_t sb_in_data_stride = 2; // in 32-bit words
          size_t sb_in_msk_stride = 2*25; // in 32-bit words


          for (x=0;x<5;x++) {
              copy_sharing(
                      NSHARES,
                      sb_state + x * sb_state_data_stride, sb_state_msk_stride,
                      sb_in + ((x+1)%5) * sb_in_data_stride, sb_in_msk_stride
                      );
              sb_state[x*sb_state_data_stride] = ~sb_state[x*sb_state_data_stride]; // NOT: on a single share
              masked_and(
                      NSHARES,
                      sb_state + x * sb_state_data_stride, sb_state_msk_stride,
                      sb_state + x * sb_state_data_stride, sb_state_msk_stride,
                      sb_in + ((x+2)%5) * sb_in_data_stride, sb_in_msk_stride
                      );
          }
          for (x=0;x<5;x++) {
              masked_xor(
                      NSHARES,
                      sb_in + x * sb_in_data_stride, sb_in_msk_stride,
                      sb_in + x * sb_in_data_stride, sb_in_msk_stride,
                      sb_state + x * sb_state_data_stride, sb_state_msk_stride
                      );
          }
      }
      // Iota
      // Add constant: on a single share
      state[0][0] ^= RC[i];
  }
    unmask_state(umsk, state);
    disp_keccak_state(umsk);
}
