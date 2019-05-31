#ifndef Mline
#define Mline 120
#endif

#ifndef tghfunc
#define tghfunc(x) ( ((x) > (0.0)) ? (tanh(x)) : (0.0) )
#endif

#ifndef linthr
#define linthr(x) ( ((x) > (0.0)) ? (x) : (0.0) )
#endif

typedef struct nrn_par{
  double beta, theta, Mx, Iapp, nu, IappF;
} nrn_par;

typedef struct syn_par{
  double gg, taus, kf, UV, taur, tauf;
} syn_par;

typedef struct net_par{
  nrn_par R, L, F;
  syn_par RR, LR, RL, FR, RF, FF, FL, LF;
  double sig;
  char tghlin, sfact, nuif;
} net_par;

typedef struct run_par{
  double epsilon;
  double deltat, nt;
  int twrite, tmcol, ttrans, tupdown;
  int sm;
  char method, incond, fpcal, smforce;
} run_par;

typedef struct Varb_ar{
  double *RR, *LR, *RL, *FR, *RF, *FF, *FL, *LF;
  double *ptr;
  int neq;
} Varb_ar;

typedef struct stat_run{
  double MRmin, MRmax;
  double MLmin, MLmax;
  double MFmin, MFmax;
  double MRthr;
  double *tud;
  double TT, freq, Tup, dutyc;
  int mud, nud, iud;
  char *sud;
} stat_run;


/* Function Declaration */

void read_input(net_par *netpar, run_par *runpar, int sm, fl_st fl);
void create_variables(Varb_ar *Varb, int sm, fl_st fl);
void free_variables(Varb_ar *Varb, int sm, fl_st fl);
void in_con(Varb_ar *Varbar, net_par *netpar, run_par *runpar,
     int sm, fl_st fl);
void n_run(Varb_ar *Varbar, net_par *netpar, run_par *runpar,
     int sm, avr_val *av, fl_st fl);
void one_integration_step(net_par *netpar, run_par *runpar, Varb_ar *Varbar,
     Varb_ar *kin, Varb_ar *kout, double delt, int it, double time, 
     Varb_ar *Varc, fl_st fl);
double MRcal(net_par *netpar, run_par *runpar, double sRR, double sRL, 
       double sRF, fl_st fl);
double MLcal(net_par *netpar, run_par *runpar, double sLR, double sLF,
       fl_st fl);
double MFcal(net_par *netpar, run_par *runpar, double sFR, double sFL,
       double sFF, fl_st fl);
double Mcal(double Iall, double beta, double Mx, double sig, char tghlin,
       fl_st fl);
void one_syn_cal(double *Varo_syn, double *kout_syn, double MM, syn_par synpar,
     double time, int it, run_par *runpar, fl_st fl);
void pr_fct(Varb_ar *Varbar, net_par *netpar, run_par *runpar, double time,
     int it, fl_st fl);
void stat_ud_cal(stat_run *statrun, net_par *netpar, run_par *runpar,
     fl_st fl);
double lininter(double x1, double x2, double xc, double y1, double y2);
