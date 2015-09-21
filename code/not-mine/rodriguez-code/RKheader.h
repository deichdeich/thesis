#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define csc(x)
#define cot(x)
#define E
#define PI
#define PIover2
#define SQRT2
#define SOL
#define SOLARMASS
#define GRAV_CONST
(1/sin((double)(x)))
(1/tan((double)(x)))
2.71828182845904523536029
3.14159265358979323846264
1.57079632679489661923132
1.41421356237309504880168
299792458.
1476.71618
7.53462359e-21
/*The sigma values for 2 sigma*/
#define ENERGY_SIGMA
#define JZ_SIGMA
#define CARTER_SIGMA
#define THETA_SIGMA
#define R_SIGMA
1.41421337/*0.207107*/
1690.01
13313.1
0.0957049
4.9445e6

/*The Cash-Karp Parameters for Embedded Runga-Kutta*/
#define A_2 0.2
#define A_3 0.3
#define A_4 0.6
#define A_5 1.
#define A_6 0.875
#define B_21 0.2
#define B_31 0.075
#define B_32 0.225
#define B_41 0.3
#define B_42 (-0.9)
#define B_43 1.2
#define B_51 (-0.203703704)
#define B_52 2.5
#define B_53 2.592592593
#define B_54 (-1.296296296)
#define B_61 0.029495804
#define B_62 0.341796875
#define B_63 0.041594329
#define B_64 0.400345414
#define B_65 0.061767578
#define C_1 0.097883598
#define C_2 0
#define C_3 0.40257649
#define C_4 0.21043771
#define C_5 0
#define C_6 0.289102202
#define C_1STAR  0.102177373
#define C_2STAR 0
#define C_3STAR 0.383907903
#define C_4STAR 0.244592737
#define C_5STAR 0.019321987
#define C_6STAR 0.25
#define CDIFF_1 (-0.004293775)
#define CDIFF_2 0
#define CDIFF_3 0.018668587
#define CDIFF_4 (-0.034155027)
#define CDIFF_5 (-0.019321987)
#define CDIFF_6 0.039102202