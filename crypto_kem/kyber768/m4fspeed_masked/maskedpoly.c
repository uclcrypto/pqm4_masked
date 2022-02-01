#include "poly.h"
#include "maskedpoly.h"
#include "masked.h"

#include "cbd.h"
#include "ntt.h"
#include "params.h"
#include "symmetric.h"

#include <stdint.h>

/*************************************************
* Name:        maskedpoly_ntt
*
* Description: Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in bitreversed order
*
* Arguments:   - uint16_t *r: pointer to in/output polynomial
**************************************************/
void maskedpoly_ntt(StrAPoly r) {
    for(int d=0;d<NSHARES;d++){
        ntt(r[d]);
    }
}

/*************************************************
* Name:        maskedpoly_invntt
*
* Description: Computes inverse of negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in bitreversed order, output in normal order
*
* Arguments:   - uint16_t *a: pointer to in/output polynomial
**************************************************/
void maskedpoly_invntt(StrAPoly r) {
    for(int d=0; d<NSHARES;d++){
        invntt(r[d]);
    }
}
