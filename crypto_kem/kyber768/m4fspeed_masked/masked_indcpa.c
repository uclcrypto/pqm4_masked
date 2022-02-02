#include "indcpa.h"
#include "ntt.h"
#include "poly.h"
#include "polyvec.h"
#include "randombytes.h"
#include "symmetric.h"
#include "masked.h"
#include "masked_utils.h"
#include "masked_poly.h"

#include <string.h>
#include <stdint.h>

extern void doublebasemul_asm_acc(int16_t *r, const int16_t *a, const int16_t *b, int16_t zeta);
/*************************************************
 * Name:        masked matacc
 *
 * Description: Multiplies a row of A or A^T, generated on-the-fly,
 *              with a vector of polynomials and accumulates into the result.
 *
 * Arguments:   - StrAPoly *r_masked:                    pointer to output polynomial to accumulate in
 *              - StrAPolyVec *b_masked:                 pointer to input vector of polynomials to multiply with
 *              - unsigned char i:            byte to indicate the index < KYBER_K of the row of A or A^T
 *              - const unsigned char *seed:  pointer to the public seed used to generate A
 *              - int transposed:             boolean indicatin whether A or A^T is generated
 **************************************************/
static void masked_matacc(StrAPoly r_masked, StrAPolyVec b_masked, unsigned char i, const unsigned char *seed, int transposed) {
    unsigned char buf[XOF_BLOCKBYTES+2];
    unsigned int buflen, off;
    xof_state state;
    unsigned int ctr, pos, k, l, d;
    uint16_t val0, val1;
    int16_t c[4];

    for(d=0;d<NSHARES;d++){
        for(l=0;l<KYBER_N;l++){
            r_masked[d][l] = 0;
        }
    }

    for(int j=0;j<KYBER_K;j++) {
        ctr = pos = 0;
        if (transposed)
            xof_absorb(&state, seed, i, j);
        else
            xof_absorb(&state, seed, j, i);

        xof_squeezeblocks(buf, 1, &state);
        buflen = XOF_BLOCKBYTES;

        k = 0;
        while (ctr < KYBER_N/4)
        {
            val0 = ((buf[pos+0] >> 0) | ((uint16_t)buf[pos+1] << 8)) & 0xFFF;
            val1 = ((buf[pos+1] >> 4) | ((uint16_t)buf[pos+2] << 4)) & 0xFFF;
            pos += 3;

            if (val0 < KYBER_Q) {
                c[k++] = (int16_t) val0;
                if (k == 4) {
                    for(d=0;d<NSHARES;d++){
                        doublebasemul_asm_acc(&r_masked[d][4*ctr], &b_masked[j][d][4*ctr], c, zetas[ctr]);
                    }
                    ctr++;
                    k = 0;
                }
            }

            if (val1 < KYBER_Q && ctr < KYBER_N/4) {
                c[k++] = (int16_t) val1;
                if (k == 4) {
                    for(d=0;d<NSHARES;d++){
                        doublebasemul_asm_acc(&r_masked[d][4*ctr], &b_masked[j][d][4*ctr], c, zetas[ctr]);
                    }
                    ctr++;
                    k = 0;
                }
            }

            if (pos + 3 > buflen && ctr < KYBER_N/4) {
                off = buflen % 3;
                for(l = 0; l < off; l++)
                    buf[l] = buf[buflen - off + l];
                xof_squeezeblocks(buf + off, 1, &state);
                buflen = off + XOF_BLOCKBYTES;
                pos = 0;
            }
        }
    }
}

/*************************************************
 * Name:        masked_indcpa_enc_cmp
 *
 * Description: Re-encryption function.
 *              Compares the re-encypted ciphertext with the original ciphertext byte per byte.
 *              The comparison is performed in a constant time manner.
 *
 *
 * Arguments:   - unsigned char *ct:         pointer to input ciphertext to compare the new ciphertext with (of length KYBER_INDCPA_BYTES bytes)
 *              - const unsigned char *m:    pointer to input message (of length KYBER_INDCPA_MSGBYTES bytes)
 *              - const unsigned char *pk:   pointer to input public key (of length KYBER_INDCPA_PUBLICKEYBYTES bytes)
 *              - const unsigned char *coin: pointer to input random coins used as seed (of length KYBER_SYMBYTES bytes)
 *                                           to deterministically generate all randomness
 * Returns:     - boolean byte indicating that re-encrypted ciphertext is NOT equal to the original ciphertext
 **************************************************/
unsigned char masked_indcpa_enc_cmp(const unsigned char *c,
        const unsigned char *m,
        const unsigned char *pk,
        const unsigned char *coins) {
    uint64_t rc = 0;
    polyvec sp;
    poly bp;
    poly *pkp = &bp;
    poly *k = &bp;
    poly *v = &sp.vec[0];
    const unsigned char *seed = pk+KYBER_POLYVECBYTES;
    int i,d ;
    unsigned char nonce = 0;
    
    StrAPolyVec masked_sp;
    StrAPoly masked_bp;
    StrAPoly masked_v;

    for (i = 0; i < KYBER_K; i++){

        // TODO protected this noise sampling
        poly_getnoise(sp.vec + i, coins, nonce++);
        masked_poly(masked_sp[i], sp.vec + i);
        
        masked_poly_ntt(masked_sp[i]);
    }

    for (i = 0; i < KYBER_K; i++) {

        masked_matacc(masked_bp, masked_sp, i, seed, 1);
        masked_poly_invntt(masked_bp);
        unmasked_poly(&bp,masked_bp);


        // TODO protect this noise sampling
        poly_addnoise(&bp, coins, nonce++);
        poly_reduce(&bp);

        // TODO protect polynomial comparison
        rc |= cmp_poly_packcompress(c, &bp, i);
    }

    poly_frombytes(pkp, pk);
    for(d=0;d<NSHARES;d++){
        poly_basemul_i16(masked_v[d], pkp->coeffs, masked_sp[0][d]);
    }

    for (i = 1; i < KYBER_K; i++) {
        poly_frombytes(pkp, pk + i*KYBER_POLYBYTES);
        for(d=0;d<NSHARES;d++){
            poly_basemul_acc_i16(masked_v[d], pkp->coeffs, masked_sp[i][d]);
        }
    }

    masked_poly_invntt(masked_v);
    unmasked_poly(v,masked_v);

    poly_addnoise(v, coins, nonce++);
    poly_frommsg(k, m);
    poly_add(v, v, k);
    poly_reduce(v);

    rc |= cmp_poly_compress(c + KYBER_POLYVECCOMPRESSEDBYTES, v);

    rc = ~rc + 1;
    rc >>= 63;
    return (unsigned char)rc;
}
