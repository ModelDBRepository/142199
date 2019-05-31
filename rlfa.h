#ifndef min
#define min(x,y) ( ((x) < (y)) ? (x) : (y) )
#define max(x,y) ( ((x) > (y)) ? (x) : (y) )
#endif

#ifndef Pi
#define Pi 3.1415926535897931
#endif

#ifndef NR_END
#define NR_END 1
#endif

#ifndef FREE_ARG
#define FREE_ARG char*
#endif

#ifndef div_nz
#define div_nz(x,y) ( (fabs(y) < (runpar.epsilon)) ? (0) : (x/y) )
#endif

#ifndef Gammaf
#define Gammaf(VV, theta,sigma) ( 1.0/(1.0+exp(-(VV-(theta))/(sigma))) )
#endif


/* Structure Declaration */

typedef struct fl_st{
  FILE *in, *tmp, *out, *col, *avr;
} fl_st;

typedef struct avr_val{
  double MRend, MLend, MFend;
  double MLmin, MLmax;
  double MFmin, MFmax;
  double MRmin, MRmax;
  double TT, freq, Tup, dutyc;
  double t_ML_pos;
  char pattern;
} avr_val;

