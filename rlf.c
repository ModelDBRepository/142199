/* This program simulates a rate model coupled RS and LTS neuronal    */
/* populations.                                                       */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "rlfa.h"
#include "rlf.h"

main(int argc, char*argv[])
{
  scan_val svala, sval;
  avr_val av;
  fl_st fl;
  int sm, savr;
  int nraw;
  char suffix[3]="a1", file_name[36], fnm[8];
  char **ar_input;

  ar_input = cmatrix(1, Mraw, 0, Mline-1);

  if (argc >= 2) strcpy(suffix,argv[1]);

  fl.in  = fopen(strcat(strcat(strcpy(file_name,"rlf."),"n."  ),suffix), "r");
  fl.tmp = fopen(strcat(strcat(strcpy(file_name,"rlf."),"tmp."),suffix), "w+");
  fl.out = fopen(strcat(strcat(strcpy(file_name,"rlf."),"out."),suffix), "w");
  fl.col = fopen(strcat(strcat(strcpy(file_name,"rlf."),"col."),suffix), "w+");

  fscanf(fl.in, "scan=%c\n", &sval.scan_type);

  if (sval.scan_type == 'n')
  {
    /* Simulating one parameter set */
    savr = sm = 0;
    fscanf(fl.in, "seed=%d\n", &sval.seed);
    fprintf(fl.out, "seed=%d\n", sval.seed);
    fl.avr = fl.out;

    read_file_old(ar_input, 2, &nraw, fl);
    write_file_old(ar_input, nraw, fl);

    one_par(fl, sm, &av);
    write_avr(svala, sval, av, savr, fl);

  }
  else if ((sval.scan_type == 'e') || (sval.scan_type == 'u'))
  {
    /* Simulating several parameter sets */
    savr = sm = 1;
    read_first_input_line(&sval, fl);
    fscanf(fl.in, "seed=%d\n", &sval.seed);
    fprintf(fl.out, "seed=%d\n", sval.seed);
    fl.avr = fopen(strcat(strcat(strcpy(file_name,"rlf."),"avr."),suffix), "w");

    for (sval.ipar=0; sval.ipar <= sval.npar; sval.ipar++)
    {
      for (sval.irepeat=1; sval.irepeat <= sval.nrepeat; sval.irepeat++)
      {
        printf("ipar=%d npar=%d irepeat=%d nrepeat=%d\n", sval.ipar, 
        sval.npar, sval.irepeat, sval.nrepeat);
        if (sval.ipar > 0) fprintf(fl.out, "\n");

        read_file_old(ar_input, 3, &nraw, fl);
        update_file_old(&sval, ar_input, nraw, fl);
        write_file_old(ar_input, nraw, fl);
 
        one_par(fl, sm, &av);
        write_avr(svala, sval, av, savr, fl);
        fflush(fl.avr);
      }
    }
    fclose(fl.avr);
  }
  else if ((sval.scan_type == 't'))
  {
    /* Scanning two parameters */
    savr = sm = 1;
    svala.scan_type = sval.scan_type;
    read_first_input_line(&svala, fl);
    read_first_input_line(&sval, fl);
    fscanf(fl.in, "seed=%d\n", &sval.seed);
    fprintf(fl.out, "seed=%d\n", sval.seed);
    fl.avr = fopen(strcat(strcat(strcpy(file_name,"rlf."),"avr."),suffix), "w");

    for (svala.ipar=0; svala.ipar <= svala.npar; svala.ipar++)
    {
      for (sval.ipar=0; sval.ipar <= sval.npar; sval.ipar++)
      {
        for (sval.irepeat=1; sval.irepeat <= sval.nrepeat; sval.irepeat++)
        {
          printf("svala: ipar=%d npar=%d\n", svala.ipar, svala.npar);
          printf("sval: ipar=%d npar=%d\n", sval.ipar, sval.npar);
          printf("irepeat=%d nrepeat=%d\n", sval.irepeat, sval.nrepeat);
          if (svala.ipar > 0 || sval.ipar) fprintf(fl.out, "\n");

          read_file_old(ar_input, 4, &nraw, fl);
          update_file_old(&svala, ar_input, nraw, fl);
          update_file_old(&sval, ar_input, nraw, fl);
          write_file_old(ar_input, nraw, fl);

	printf("a ipar=%d ipar=%d irepeat=%d\n", svala.ipar,
          sval.ipar, sval.irepeat);
          one_par(fl, sm, &av);
          write_avr(svala, sval, av, savr, fl);
          fflush(fl.avr);
        }
      }
    }
    fclose(fl.avr);
  }
  else if ((sval.scan_type == 's'))
  {
    /* Scanning one parameter and searching another one */
    savr = sm = 1;
    svala.scan_type = sval.scan_type;
    read_first_input_line(&svala, fl);
    read_second_input_line(&sval, fl);
    fscanf(fl.in, "seed=%d\n", &sval.seed);
    fprintf(fl.out, "seed=%d\n", sval.seed);
    fl.avr = fopen(strcat(strcat(strcpy(file_name,"rlf."),"avr."),suffix), "w");
    for (svala.ipar=0; svala.ipar <= svala.npar; svala.ipar++)
    {
      printf("svala: ipar=%d npar=%d\n", svala.ipar, svala.npar);

      find_borders(ar_input, &svala, &sval, sm, &av, fl);
     
      /* write_avr(svala, sval, av, savr, fl); */
      fflush(fl.avr);
    }
    fclose(fl.avr);
  }


  free_cmatrix(ar_input, 1, Mraw, 0, Mline-1);

  fclose(fl.in);
  fclose(fl.tmp);
  fclose(fl.out);
  fclose(fl.col);
}

void find_borders(char **ar_input, scan_val *svala, scan_val *sval,
     int sm, avr_val *av, fl_st fl)
{
  double par_ar[Mpar], parl, parr;
  int jpar;
  char spk_pat_l, spk_pat_r;
  

  sval->ipar = 0;

  if (sval->npar == 0)
  {
    printf("npar=%d does not fit scan=s\n", sval->npar);
    exit (0);
  }
  else
  {
    for (jpar=0; jpar<=sval->npar; jpar++)
    {
      par_ar[jpar] = sval->parmin + (sval->parmax - sval->parmin) *
      jpar / sval->npar;
    } 
  }

  parl = par_ar[0];
  sval->par_ar[0] = parl;
  printf("jpar=0 par=%lf\n", sval->par_ar[0]);
  update_and_run(ar_input, svala, sval, sm, av, fl);
  printf("pattern=%c\n", av->pattern);
  spk_pat_l = av->pattern;

  for (jpar=1; jpar<=sval->npar; jpar++)
  {
    parr = par_ar[jpar];
    sval->par_ar[0] = parr;
    printf("jpar=%d par=%lf\n", jpar, sval->par_ar[0]);
    update_and_run(ar_input, svala, sval, sm, av, fl);
    printf("pattern=%c\n", av->pattern);
    spk_pat_r = av->pattern;

    printf("pat=%c %c %c\n", spk_pat_l, sval->ptype,  spk_pat_r);
    if ((spk_pat_l == sval->ptype) && (spk_pat_r != sval->ptype))
    {
      printf("l!\n");
      find_one_border(par_ar[jpar-1], par_ar[jpar], 
      ar_input, svala, sval, sm, av, fl);
    }
    else if ((spk_pat_l != sval->ptype) && (spk_pat_r == sval->ptype))
    {
      printf("r!\n");
      find_one_border(par_ar[jpar], par_ar[jpar-1],
      ar_input, svala, sval, sm, av, fl);
    }
    else
    {
      printf("n!\n");
    }

    spk_pat_l = spk_pat_r;
  }
}

void find_one_border(double xl, double xr, 
     char **ar_input, scan_val *svala, scan_val *sval,
     int sm, avr_val *av, fl_st fl)
{
  double xm;
  char spk_pat_m;

  while (fabs(xr - xl) > sval->dpar/2.0)
  {
    xm = (xl + xr) / 2.0; 
    printf("xl=%lf xr = %lf xm=%lf dpar=%lf\n", xl, xr, xm, sval->dpar);
    sval->par_ar[0] = xm;
    update_and_run(ar_input, svala, sval, sm, av, fl);
    printf("pattern=%c\n", av->pattern);
    spk_pat_m = av->pattern;

    if (spk_pat_m == sval->ptype)
    {
      xl = xm;
    }
    else
    {
      xr = xm;
    }
  }

  printf("res: %lf %lf %c\n", svala->par_ar[svala->ipar], xm, sval->ptype);
  fprintf(fl.avr, "%lf %lf %c\n", svala->par_ar[svala->ipar], xm, sval->ptype);

  fflush(fl.avr);
}

void update_and_run(char **ar_input, scan_val *svala, scan_val *sval,
     int sm, avr_val *av, fl_st fl)
{
  int nraw;

  read_file_old(ar_input, 4, &nraw, fl);
  update_file_old(svala, ar_input, nraw, fl);
  update_file_old(sval, ar_input, nraw, fl);
  write_file_old(ar_input, nraw, fl);

  one_par(fl, sm, av);
}

/* This function reads and processes the first input line */
void read_first_input_line(scan_val *sval, fl_st fl)
{
  int ipar;
  char line[Mline], *p1, par2all[Mword], *p2, *p3;

  if (fgets(line, Mline, fl.in) == NULL)
  {
    printf("empty input file!!!\n");
    exit(0);
  }
  
  p1 = strtok(line, " ");
  sscanf(p1, " %s", sval->par1);
  p1 = strtok(NULL, " ");
  sscanf(p1, " %s", par2all);

  if (strstr(par2all, "#") == NULL)
  {
    strcpy(sval->par2, par2all);
    sval->npt = 1;
  }
  else
  {
    p3 = &sval->par2[0];
    for (p2 = &par2all[0]; *p2 != '#'; p2++) *p3++ = *p2;
    *p3 = '\0';
    p2++;
   sscanf(p2, "%d", &sval->npt);
  }

  fprintf(fl.out, " par1=%s par2=%s npt=%d", sval->par1, sval->par2,
  sval->npt);

  if ((sval->scan_type == 'e') || (sval->scan_type == 't') ||
      (sval->scan_type == 's'))
  {
    p1 = strtok(NULL, " ");
    sscanf(p1, "parmin=%lf", &sval->parmin);
    p1 = strtok(NULL, " ");
    sscanf(p1, "parmax=%lf", &sval->parmax);
    p1 = strtok(NULL, " ");
    sscanf(p1, "npar=%d", &sval->npar);
    p1 = strtok(NULL, " ");
    sscanf(p1, "nrepeat=%d", &sval->nrepeat);
    fprintf(fl.out, " parmin=%lf parmax=%lf\n npar=%d nrepeat=%d\n", 
    sval->parmin, sval->parmax, sval->npar, sval->nrepeat);

    if (sval->npar == 0)
    {
      sval->par_ar[0] = sval->parmin;
    }
    else
    {
      for (ipar=0; ipar<=sval->npar; ipar++)
      {
        sval->par_ar[ipar] = sval->parmin + (sval->parmax - sval->parmin) *
        ipar / sval->npar;
      } 
    }
  }
  else if (sval->scan_type == 'u')
  {
    ipar=-1;
    while((p1 = strtok(NULL, " ")) != NULL)
    {
      sscanf(p1, "%lf", &sval->par_ar[++ipar]);
    }
    sval->npar=ipar;
  }
  else
  {
    printf("wrong scan_type=%c!!!\n", sval->scan_type);
    exit(0);
  }

  for (ipar=0; ipar<=sval->npar; ipar++)
    fprintf(fl.out, "ipar=%d par_ar=%lf\n", ipar, sval->par_ar[ipar]);

  return;
}

/* This function reads and processes the second input line */
void read_second_input_line(scan_val *sval, fl_st fl)
{
  int ipar;
  char line[Mline], *p1, par2all[Mword], *p2, *p3;

  if (fgets(line, Mline, fl.in) == NULL)
  {
    printf("empty input file!!!\n");
    exit(0);
  }
  
  p1 = strtok(line, " ");
  sscanf(p1, " %s", sval->par1);
  p1 = strtok(NULL, " ");
  sscanf(p1, " %s", par2all);

  if (strstr(par2all, "#") == NULL)
  {
    strcpy(sval->par2, par2all);
    sval->npt = 1;
  }
  else
  {
    p3 = &sval->par2[0];
    for (p2 = &par2all[0]; *p2 != '#'; p2++) *p3++ = *p2;
    *p3 = '\0';
    p2++;
   sscanf(p2, "%d", &sval->npt);
  }

  fprintf(fl.out, "par1=%s par2=%s npt=%d", sval->par1, sval->par2,
  sval->npt);

  p1 = strtok(NULL, " ");
  sscanf(p1, "parmin=%lf", &sval->parmin);
  p1 = strtok(NULL, " ");
  sscanf(p1, "parmax=%lf", &sval->parmax);
  p1 = strtok(NULL, " ");
  sscanf(p1, "npar=%d", &sval->npar);
  p1 = strtok(NULL, " ");
  sscanf(p1, "dpar=%lf", &sval->dpar);
  p1 = strtok(NULL, " ");
  sscanf(p1, "ptype=%c", &sval->ptype);
  fprintf(fl.out, " parmin=%lf parmax=%lf\nnpar=%d dpar=%lf ptype=%c\n", 
  sval->parmin, sval->parmax, sval->npar, sval->dpar, sval->ptype);

  return;
}

/* This function reads the input file into the array ar_input         */
void read_file_old(char **ar_input, int skip_lines, int *nraw, fl_st fl)
{
  int iraw;
  char line[Mline];

  rewind(fl.in);
  rewind(fl.tmp);

  for (iraw=1; iraw<=skip_lines; iraw++) fgets(line, Mline, fl.in);

  iraw = 0;
  while (fgets(line, Mline, fl.in) != NULL) 
  {
    iraw++;
    strcpy(ar_input[iraw], line);
  }

  *nraw = iraw;
}

/* This function updates the array ar_input                           */
void update_file_old(scan_val *sval, char **ar_input, int nraw, fl_st fl)
{
  int iraw, jraw, nchange;


  /* scanning - multiplying all the occurrences of the specific parameter */
  /* value                                                                */

  if (strcmp(sval->par1, "ALL") == 0)
  {
    nchange=0;

    for (iraw=1; iraw<= nraw; iraw++)    
    {
      if (process_line_old(sval, ar_input[iraw], fl.tmp)) nchange++;
    }
    return;
  }

  /* scanning - changing the specific parameter value */

  iraw = 1;
  while (iraw <= nraw)
  {
    if (strncmp(ar_input[iraw], sval->par1, strlen(sval->par1)) != 0)
      iraw++;
    else 
      break;
  }

  if (iraw > nraw)
  {
    printf("par1 not found!!!\n");
    exit(0);
  }

  for (jraw=iraw; jraw<=nraw; jraw++)
  {
    if (process_line_old(sval, ar_input[jraw], fl.tmp))
    {
       break;
    }
  }

  /* checking for end of file */

  if (jraw > nraw)
  {
    printf("match not found!!!\n");
    exit(0);
  }
}

/* This function writes the array ar_input on the tmp file            */
void write_file_old(char **ar_input, int nraw, fl_st fl)
{
  int iraw;

  rewind(fl.tmp);

  for (iraw=1; iraw<=nraw; iraw++) fprintf(fl.tmp, "%s", ar_input[iraw]);

  rewind(fl.tmp);
}

/* This function processes one line and relplace a value by sval.par2 or */
/* multiplies it by sval->par2                                           */
int process_line_old(scan_val *sval, char line[], FILE *ftmp)
{
  int ret_val;

  if (strstr(sval->par2, ":") == NULL)
  {
    ret_val = process_line_no_colon_old(sval, line);
  }
  else
  {
    ret_val = process_line_yes_colon_old(sval, line);
  }

  return(ret_val);
}

/* This function processes one line and relplace a value by sval.par2 or */
/* multiplies it by sval->par2 .                                         */
/* There is no colon (:) in sval->par2 .                                 */
int process_line_no_colon_old(scan_val *sval, char line[])
{
  double par_ref;
  int  il, il_end, im, ipt, iw, inew;
  int cond1, cond2, cond3;
  char *pline, line_new[Mline];

  il_end = -1;
  while (line[il_end+1] != '\n' && il_end < Mline-2) il_end++;
  pline = line;
  il = -1;

  while (++il <= il_end)
  {
    /* Condition: matched pattern, '=' at the end, ' ' or beginning of  */
    /* line at the beginning                                            */
    cond1 = strncmp(pline+il, sval->par2, strlen(sval->par2)) == 0;
    cond2 = line[il + strlen(sval->par2)] == '=';
    if (il == 0)
      cond3 = 1;
    else
      cond3 = line[il - 1] == ' ';
    if (cond1 && cond2 && cond3) break;
  }

  if (il >= il_end-1)
  /* par2 does not appear in line */
  {
    return(0);
  }
  else
  /* par2 appears in line */
  {
    inew = -1;
    for (im=0; im<=il+strlen(sval->par2); im++) 
    {
      inew++;
      line_new[inew] = line[im];
    }

    while (line[il-1] != '=') il++;
    ipt=0;
    while (ipt < sval->npt-1)
    {
      il++;
      inew++;
      line_new[inew] = line[il];
      if (line[il-1] != ' ' && line[il] == ' ') ipt++;
    }
    
    while (line[il] == ' ') il++;
    iw=-1;
    while ((line[il] != ' ') && (il <= il_end))
    {
      il++;
    }

    if (strcmp(sval->par1, "ALL") != 0)
      sprintf(&line_new[inew+1], "%lf", sval->par_ar[sval->ipar]);
    else
    {
      sprintf(&line_new[inew+1], "%lf", sval->par_ar[sval->ipar] * par_ref);
    }

    inew = strlen(line_new);

    for (im=il; im<=il_end; im++)
    {
      line_new[inew-il+im] = line[im];
    }

    line_new[inew-il+il_end+1] = '\n';
    line_new[inew-il+il_end+2] = '\0';
    strcpy(line, line_new);

    return(1);
  }
}

/* This function processes one line and relplace a value by sval.par2 or */
/* multiplies it by sval->par2 .                                         */
/* There is a colon (:) in sval->par2 .                                  */
int process_line_yes_colon_old(scan_val *sval, char line[])
{
  double par_ref;
  int  il, il_end, im, ipt, iw, inew;
  int cond1, cond2, cond3;
  int len2a, len2an, num_par, count_par;
  char line_cp[Mword], par2a[Mword], par2an[Mword], par2b[Mword];
  char *pline, *p1, line_new[Mline];;

  il_end = -1;
  while (line[il_end+1] != '\n' && il_end < Mline-2) il_end++;
  pline = line;

  strcpy(line_cp, sval->par2);
  p1 = strtok(line_cp, ":");
  sscanf(p1, " %s", par2a);
  len2a = strlen(par2a);

  strcat(strcpy(par2an, par2a),":");
  p1 = strtok(NULL, " ");
  sscanf(p1, " %s", par2b);
  len2an = strlen(par2an);
  sscanf(par2b, "%d", &num_par);

  /* Checking for par2a */

  il = -1;

  while (++il <= il_end)
  {
    /* Condition: matched pattern, '=' at the end, ' ' or beginning of  */
    /* line at the beginning                                            */
    cond1 = strncmp(pline+il, par2a, len2a) == 0;
    cond2 =  *(pline + il + len2a) == '=';
    if (il == 0)
      cond3 = 1;
    else
      cond3 = line[il - 1] == ' ';
    if (cond1 && cond2 && cond3) break;
  }

  if (il >= il_end-1)
  /* par2 does not appear in line */
  {
    return(0);
  }

  if (num_par == 1)
  {
   while (++il <= il_end)
    {
      if (line[il] == '=') break;
    }
  }
  else
  {
    /* skipping the parameter values before num_par */
    count_par = 1;
    while (++il <= il_end)
    {
      if (line[il] == ' ') count_par++;
      if (count_par == num_par) break;
    }
  }

  /* printing the line before the parameter that is updated */
 
  inew = -1;
  for (im=0; im<=il; im++) 
  {
    inew++;
    line_new[inew] = line[im];
  }

 /* skipping the parameter that is updated */
  while (++il <= il_end)
  {
    if (line[il] == ' ') break;
  }

  if (strcmp(sval->par1, "ALL") != 0)
    sprintf(&line_new[inew+1], "%lf", sval->par_ar[sval->ipar]);
  else
  {
    sprintf(&line_new[inew+1], "%lf", sval->par_ar[sval->ipar] * par_ref);
  }

  inew = strlen(line_new);
  for (im=il; im<=il_end; im++)
  {
    line_new[inew-il+im] = line[im];
  }

  line_new[inew-il+il_end+1] = '\n';
  line_new[inew-il+il_end+2] = '\0';

  strcpy(line, line_new);

  return(1);
}

/* This function writes population-average results on fl.avr */
void write_avr(scan_val svala, scan_val sval, avr_val av, int savr, fl_st fl)
{
  fprintf(fl.out, "\nMRend=%lf MLend=%lf\n", av.MRend, av.MLend);

  if (savr >= 1) 
  {
     if (sval.scan_type == 't') 
       fprintf(fl.avr, "%lf ", svala.par_ar[svala.ipar]);

    fprintf(fl.avr, "%lf", sval.par_ar[sval.ipar]);
    /* fprintf(fl.avr, "%d", sval.irepeat); */
  }

  fprintf(fl.avr, " %lf %lf %lf", av.MRend, av.MLend, av.MFend);
  fprintf(fl.avr, " %lf", av.t_ML_pos);
  fprintf(fl.avr, " %lf %lf ", av.MRmin, av.MRmax);
  fprintf(fl.avr, " %lf %lf %lf %lf", av.TT, av.freq, av.Tup, av.dutyc);
  fprintf(fl.avr, " %c", av.pattern);
  fprintf(fl.avr, " %lf %lf ", av.MLmin, av.MLmax);
  fprintf(fl.avr, " %lf %lf ", av.MFmin, av.MFmax);

  fprintf(fl.avr, "\n");
  fflush(fl.avr);
}
