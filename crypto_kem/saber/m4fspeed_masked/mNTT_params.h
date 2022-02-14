#ifndef MNTT_PARAMS_H
#define MNTT_PARAMS_H

#define ARRAY_N 256

#define NTT_N 64
#define LOGNTT_N 6

#define Q1 44683393
#define Q1pr 7
// omegaQ1 = Q1pr^((Q1 - 1) / (NTT_N << 1)) mod Q1
#define omegaQ1 22875492
// invomegaQ1 = omegaQ1^{-1} mod Q1
#define invomegaQ1 7412207
// RmodQ1 = 2^32 mod^{+-} Q1
#define RmodQ1 5361568
// Q1prime = -Q1^{-1} mod^{+-} 2^32
#define Q1prime -604401537
// invNQ1 = NTT_N^{-1} mod Q1
#define invNQ1 43985215


#define Q2 769
#define Q2pr 11
// omegaQ2 = Q2pr^((Q2 - 1) / (NTT_N << 1)) mod Q2
#define omegaQ2 554
// invomegaQ2 = omegaQ2^{-1} mod Q2
#define invomegaQ2 676
// RmodQ2 = 2^16 mod^{+-} Q2
#define RmodQ2 171
// Q2prime = -Q2^{-1} mod^{+-} 2^16
#define Q2prime 767
// Q2Q2prime = Q2 || Q2prime
#define Q2Q2prime 50397951
// invNQ2 = NTT_N^{-1} mod Q2
#define invNQ2 757

#define Q1half 22341696
#define Q2invRmod 8548531
#define Q2bar 5585133

#endif



