#define Mpar 20000
#define Mraw 200
#define Mline 2000
#define Mword 50

/* Structure Declaration */

typedef struct scan_val{
  double eps;
  char scan_type, bhv;
  int  seed;

  /* double parmina, parmaxa, para_ar[Mpar]; */
  /* int npara, ipara, npta, nrepeata; */
  /* char par1a[Mword], par2a[Mword]; */

  double parmin, parmax, par_ar[Mpar];
  int npar, ipar, npt, nrepeat, irepeat;
  char par1[Mword], par2[Mword];

  double dpar;
  char ptype; 
} scan_val;


/* Function Declaration */

void find_borders(char **ar_input, scan_val *svala, scan_val *sval, 
     int sm, avr_val *av, fl_st fl);
void find_one_border(double xl, double xr, 
     char **ar_input, scan_val *svala, scan_val *sval, int sm,
     avr_val *av, fl_st fl);
void update_and_run(char **ar_input, scan_val *svala, scan_val *sval, 
     int sm, avr_val *av, fl_st fl);
void read_first_input_line(scan_val *sval, fl_st fl);
void read_second_input_line(scan_val *sval, fl_st fl);
void read_file_old(char **ar_input, int skip_lines, int *nraw, fl_st fl);
void update_file_old(scan_val *sval, char **ar_input, int nraw, fl_st fl);
void write_file_old(char **ar_input, int nraw, fl_st fl);
int process_line_old(scan_val *sval, char line[], FILE *ftmp);
int process_line_no_colon_old(scan_val *sval, char line[]);
int process_line_yes_colon_old(scan_val *sval, char line[]);
void one_par(fl_st fl, int sm, avr_val *av);
void write_avr(scan_val svala, scan_val sval, avr_val av, int savr, fl_st fl);
