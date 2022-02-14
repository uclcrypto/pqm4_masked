#include "masked_cbd.h"
#include "gadgets.h"

void masked_cbd_seed(
    size_t nshares,
    uint16_t *s, 
    size_t s_msk_stride,
    size_t s_data_stride,
    const uint8_t *buf,
    size_t buf_msk_stride,
    size_t buf_data_stride)
{

#if SABER_MU == 8 // Saber
  int offset;
  uint16_t cbd[2*BSSIZE*NSHARES];

  for(int j=0; j < SABER_N; j+=2*BSSIZE){
    uint32_t a[4*NSHARES*2];
    uint32_t b[4*NSHARES*2];

    for(int d = 0; d<(2*4*NSHARES);d++){
      a[d] = 0; b[d] = 0;
    }

    for(int n=0;n<32;n+=1){
      offset = j+n;
      for(int d=0;d<NSHARES;d++){
        uint8_t byte = buf[offset * buf_data_stride +  d*buf_msk_stride];
        a[0*NSHARES+d] |= ((byte>>0) & 0x1)<<n;
        a[1*NSHARES+d] |= ((byte>>1) & 0x1)<<n;
        a[2*NSHARES+d] |= ((byte>>2) & 0x1)<<n;
        a[3*NSHARES+d] |= ((byte>>3) & 0x1)<<n;

        b[0*NSHARES+d] |= ((byte>>4) & 0x1)<<n;
        b[1*NSHARES+d] |= ((byte>>5) & 0x1)<<n;
        b[2*NSHARES+d] |= ((byte>>6) & 0x1)<<n;
        b[3*NSHARES+d] |= ((byte>>7) & 0x1)<<n;
      }
      offset = j+n+32;
      for(int d=0;d<NSHARES;d++){
        uint8_t byte = buf[offset * buf_data_stride +  d*buf_msk_stride];
        a[0*NSHARES+d+NSHARES*4] |= ((byte>>0) & 0x1)<<n;
        a[1*NSHARES+d+NSHARES*4] |= ((byte>>1) & 0x1)<<n;
        a[2*NSHARES+d+NSHARES*4] |= ((byte>>2) & 0x1)<<n;
        a[3*NSHARES+d+NSHARES*4] |= ((byte>>3) & 0x1)<<n;

        b[0*NSHARES+d+NSHARES*4] |= ((byte>>4) & 0x1)<<n;
        b[1*NSHARES+d+NSHARES*4] |= ((byte>>5) & 0x1)<<n;
        b[2*NSHARES+d+NSHARES*4] |= ((byte>>6) & 0x1)<<n;
        b[3*NSHARES+d+NSHARES*4] |= ((byte>>7) & 0x1)<<n;
      }
    }

    masked_cbd(nshares,4,
              SABER_EQ,
              cbd,1,NSHARES,
              a,1,NSHARES,
              b,1,NSHARES);

    for(int n=0;n<32;n+=1){
      for(int d=0;d<NSHARES;d++){
        s[d*s_msk_stride + n*s_data_stride] = cbd[d + n*NSHARES];
        s[d*s_msk_stride + (n+32)*s_data_stride] = cbd[d + (n+32)*NSHARES];
      }
    }
  }
#else
#error "Unsupported SABER parameter."
#endif
}
