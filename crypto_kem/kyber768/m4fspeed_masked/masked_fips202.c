/* Based on the public domain implementation in
 * crypto_hash/keccakc512/simple/ from http://bench.cr.yp.to/supercop.html
 * by Ronny Van Keer
 * and the public domain "TweetFips202" implementation
 * from https://twitter.com/tweetfips202
 * by Gilles Van Assche, Daniel J. Bernstein, and Peter Schwabe */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "masked.h"
#include "masked_keccakf1600.h"
#include "masked_fips202.h"

#define SHAKE256_RATE 136
#define SHA3_512_RATE 72

typedef struct {
    // Keccack masked state
    MaskedKeccakState state;
    // Either the number of absorbed bytes that have not been permuted, or
    // number of not-yet-squeezed bytes.
    size_t xored_non_absorbed;
} keccak_inc_ctx;


/*************************************************
 * Name:        masked_keccak_inc_init
 *
 * Description: Initializes the incremental Keccak state to zero.
 *
 * Arguments:   - keccak_inc_ctx *s: pointer to input/output incremental state.
 **************************************************/
static void masked_keccak_inc_init(keccak_inc_ctx *s) {
    memset(s->state, 0, sizeof(s->state));
    s->xored_non_absorbed = 0;
}

/*************************************************
 * Name:        masked_keccak_inc_absorb
 *
 * Description: Incremental keccak absorb
 *              Preceded by keccak_inc_init, succeeded by keccak_inc_finalize
 *
 * Arguments:   - keccak_inc_ctx *s: pointer to input/output incremental state.
 *              - uint32_t r: rate in bytes (e.g., 168 for SHAKE128)
 *              - const uint8_t *m: pointer to input to be absorbed into s_inc
 *              - size_t mlen: length of input in bytes
 *              - size_t m_msk_stride: stride of m shares.
 *              - size_t m_data_stride: stride of m data bytes.
 **************************************************/
static void masked_keccak_inc_absorb(
        keccak_inc_ctx *s,
        uint32_t r,
        const uint8_t *m,
        size_t mlen,
        size_t m_msk_stride,
        size_t m_data_stride
        ) {
    while (mlen + s->xored_non_absorbed >= r) {
        MaskedKeccakF1600_StateXORBytes(
                s->state,
                m,
                s->xored_non_absorbed,
                r-s->xored_non_absorbed,
                m_msk_stride,
                m_data_stride
                );
        mlen -= (size_t)(r - s->xored_non_absorbed);
        m += r - s->xored_non_absorbed;
        s->xored_non_absorbed = 0;

        MaskedKeccakF1600_StatePermute(s->state);
    }
    MaskedKeccakF1600_StateXORBytes(
            s->state,
            m,
            s->xored_non_absorbed,
            mlen,
            m_msk_stride,
            m_data_stride
            );
    s->xored_non_absorbed += mlen;
}

/*************************************************
 * Name:        masked_keccak_inc_finalize
 *
 * Description: Finalizes Keccak absorb phase, prepares for squeezing
 *
 * Arguments:   - keccak_inc_ctx *s: pointer to input/output incremental state.
 *              - uint32_t r: rate in bytes (e.g., 168 for SHAKE128)
 *              - const uint8_t *m: pointer to input to be absorbed into s_inc
 *              - uint8_t p: domain-separation byte for different
 *                                 Keccak-derived functions
 **************************************************/
static void masked_keccak_inc_finalize(keccak_inc_ctx *s, uint32_t r, uint8_t p) {
    /* After keccak_inc_absorb, we are guaranteed that s->xored_non_absorbed < r,
       so we can always use one more byte for p in the current state. */
    if(s->xored_non_absorbed == r-1){
      p |= 128;
      MaskedKeccakF1600_StateXORPublicBytes(s->state, &p, s->xored_non_absorbed, 1);
    } else {
      MaskedKeccakF1600_StateXORPublicBytes(s->state, &p, s->xored_non_absorbed, 1);
      p = 128;
      MaskedKeccakF1600_StateXORPublicBytes(s->state, &p, r-1, 1);
    }
    s->xored_non_absorbed = 0;
}

/*************************************************
 * Name:        masked_keccak_inc_squeeze
 *
 * Description: Incremental Keccak squeeze; can be called on byte-level
 *
 * Arguments:   - keccak_inc_ctx *s: pointer to input/output incremental state.
 *              - uint32_t r: rate in bytes (e.g., 168 for SHAKE128)
 *              - uint8_t *h: pointer to output bytes
 *              - size_t outlen: number of bytes to be squeezed
 *              - size_t h_msk_stride: stride of h shares.
 *              - size_t h_data_stride: stride of h data bytes.
 **************************************************/
static void masked_keccak_inc_squeeze(
        keccak_inc_ctx *s,
        uint32_t r,
        uint8_t *h,
        size_t outlen,
        size_t h_msk_stride,
        size_t h_data_stride
        ) {
    size_t len;
    if(outlen < s->xored_non_absorbed) {
        len = outlen;
    } else {
        len = s->xored_non_absorbed;
    }
    MaskedKeccakF1600_StateExtractBytes(
            s->state,
            h,
            r-s->xored_non_absorbed,
            len,
            h_msk_stride,
            h_data_stride);
    h += len;
    outlen -= len;
    s->xored_non_absorbed -= len;

    /* Then squeeze the remaining necessary blocks */
    while (outlen > 0) {
        MaskedKeccakF1600_StatePermute(s->state);
        if(outlen < r) {
            len = outlen;
        } else {
            len = r;
        }
        MaskedKeccakF1600_StateExtractBytes(
                s->state,
                h,
                0,
                len,
                h_msk_stride,
                h_data_stride);
        h += len;
        outlen -= len;
        s->xored_non_absorbed = r - len;
    }
}


/*************************************************
 * Name:        masked_shake256
 *
 * Description: SHAKE256 XOF with non-incremental API
 *
 * Arguments:   - uint8_t *output:      pointer to output
 *              - size_t outlen:        requested output length in bytes
 *              - size_t out_msk_stride: stride of output shares.
 *              - size_t out_data_stride: stride of output data bytes.
 *              - const uint8_t *input: pointer to input
 *              - size_t inlen:         length of input in bytes
 *              - size_t in_msk_stride: stride of input shares.
 *              - size_t in_data_stride: stride of input data bytes.
 **************************************************/
void masked_shake256(
        uint8_t *output, size_t outlen, size_t out_msk_stride, size_t out_data_stride,
        const uint8_t *input, size_t inlen, size_t in_msk_stride, size_t in_data_stride
        ) {
  keccak_inc_ctx state;
  masked_keccak_inc_init(&state);
  /* Absorb input */
  masked_keccak_inc_absorb(&state, SHAKE256_RATE, input, inlen, in_msk_stride, in_data_stride);
  masked_keccak_inc_finalize(&state, SHAKE256_RATE, 0x1F);
  /* Squeeze output */
  masked_keccak_inc_squeeze(&state, SHAKE256_RATE, output, outlen, out_msk_stride, out_data_stride);
}

/*************************************************
 * Name:        masked_sha3_512
 *
 * Description: SHA3-512 with non-incremental API
 *
 * Arguments:   - uint8_t *output:      pointer to output
 *              - size_t out_msk_stride: stride of output shares.
 *              - size_t out_data_stride: stride of output data bytes.
 *              - const uint8_t *input: pointer to input
 *              - size_t inlen:         length of input in bytes
 *              - size_t in_msk_stride: stride of input shares.
 *              - size_t in_data_stride: stride of input data bytes.
 **************************************************/
void masked_sha3_512(
        uint8_t *output, size_t out_msk_stride, size_t out_data_stride,
        const uint8_t *input, size_t inlen, size_t in_msk_stride, size_t in_data_stride
        ) {
  keccak_inc_ctx state;
  masked_keccak_inc_init(&state);
  /* Absorb input */
  masked_keccak_inc_absorb(&state, SHA3_512_RATE, input, inlen, in_msk_stride, in_data_stride);
  masked_keccak_inc_finalize(&state, SHA3_512_RATE, 0x06);
  /* Squeeze output */
  masked_keccak_inc_squeeze(&state, SHA3_512_RATE, output, 64, out_msk_stride, out_data_stride);
}
