/* This part of the program simulates the model of one cell             */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nr.h"
#include "rlfa.h"
#include "rlfn.h"

/* This function simulates the system for one parameter set */
void one_par(fl_st fl, int sm, avr_val *av)
{
  net_par netpar;
  run_par runpar;
  Varb_ar Varbar;

  read_input(&netpar, &runpar, sm, fl);
  create_variables(&Varbar, sm, fl);
  in_con(&Varbar, &netpar, &runpar, sm, fl);
  n_run(&Varbar, &netpar, &runpar, sm, av, fl);
  free_variables(&Varbar, sm, fl);
}

/* This function reads the input parameter file */
void read_input(net_par *netpar, run_par *runpar, int sm, fl_st fl)
{
  runpar->epsilon = 1.0e-10;

  /* Reading the data */

  fscanf(fl.tmp, "NEURON\n");

  fscanf(fl.tmp, "R: beta=%lf theta=%lf Mx=%lf\n", &netpar->R.beta,
  &netpar->R.theta, &netpar->R.Mx);
  fscanf(fl.tmp, "Iapp=%lf\n", &netpar->R.Iapp);

  fscanf(fl.tmp, "L: beta=%lf theta=%lf Mx=%lf\n", &netpar->L.beta,
  &netpar->L.theta, &netpar->L.Mx);

  fscanf(fl.tmp, "F: beta=%lf theta=%lf Mx=%lf\n", &netpar->F.beta,
  &netpar->F.theta, &netpar->F.Mx);
  fscanf(fl.tmp, "IappF=%lf nu=%lf\n", &netpar->F.IappF, &netpar->F.nu);

  fscanf(fl.tmp, "SYNAPSES\n");

  fscanf(fl.tmp, "RR: gg=%lf taus=%lf kf=%lf UV=%lf taur=%lf tauf=%lf\n",
  &netpar->RR.gg, &netpar->RR.taus, &netpar->RR.kf, &netpar->RR.UV,
  &netpar->RR.taur, &netpar->RR.tauf); 

  fscanf(fl.tmp, "LR: gg=%lf taus=%lf kf=%lf UV=%lf taur=%lf tauf=%lf\n",
  &netpar->LR.gg, &netpar->LR.taus, &netpar->LR.kf, &netpar->LR.UV,
  &netpar->LR.taur, &netpar->LR.tauf); 

  fscanf(fl.tmp, "RL: gg=%lf taus=%lf kf=%lf UV=%lf taur=%lf tauf=%lf\n",
  &netpar->RL.gg, &netpar->RL.taus, &netpar->RL.kf, &netpar->RL.UV,
  &netpar->RL.taur, &netpar->RL.tauf); 

  fscanf(fl.tmp, "FR: gg=%lf taus=%lf kf=%lf UV=%lf taur=%lf tauf=%lf\n",
  &netpar->FR.gg, &netpar->FR.taus, &netpar->FR.kf, &netpar->FR.UV,
  &netpar->FR.taur, &netpar->FR.tauf); 

  fscanf(fl.tmp, "RF: gg=%lf taus=%lf kf=%lf UV=%lf taur=%lf tauf=%lf\n",
  &netpar->RF.gg, &netpar->RF.taus, &netpar->RF.kf, &netpar->RF.UV,
  &netpar->RF.taur, &netpar->RF.tauf); 

  fscanf(fl.tmp, "FF: gg=%lf taus=%lf kf=%lf UV=%lf taur=%lf tauf=%lf\n",
  &netpar->FF.gg, &netpar->FF.taus, &netpar->FF.kf, &netpar->FF.UV,
  &netpar->FF.taur, &netpar->FF.tauf); 

  fscanf(fl.tmp, "FL: gg=%lf taus=%lf kf=%lf UV=%lf taur=%lf tauf=%lf\n",
  &netpar->FL.gg, &netpar->FL.taus, &netpar->FL.kf, &netpar->FL.UV,
  &netpar->FL.taur, &netpar->FL.tauf); 

  fscanf(fl.tmp, "LF: gg=%lf taus=%lf kf=%lf UV=%lf taur=%lf tauf=%lf\n",
  &netpar->LF.gg, &netpar->LF.taus, &netpar->LF.kf, &netpar->LF.UV,
  &netpar->LF.taur, &netpar->LF.tauf); 

  fscanf(fl.tmp, "tghlin=%c sfact=%c nuif=%c sig=%lf\n", &netpar->tghlin,
  &netpar->sfact, &netpar->nuif, &netpar->sig);

  fscanf(fl.tmp, "GENERAL\n");

  fscanf(fl.tmp, "deltat=%lf nt=%lf\n", &runpar->deltat, &runpar->nt);
  fscanf(fl.tmp, "twrite=%d tmcol=%d ttrans=%d tupdown=%d\n", &runpar-> twrite,
  &runpar->tmcol, &runpar->ttrans, &runpar->tupdown);
  fscanf(fl.tmp, "method=%c incond=%c fpcal=%c smforce=%c\n", &runpar->method,
  &runpar->incond, &runpar->fpcal, &runpar->smforce);

  /* printing flag */

  runpar->sm = sm;    /* 0 - print as if there is only one parameter set */
                      /* 1 - no such printing                            */
  if (runpar->smforce == 'p')      /* always print    */
    runpar->sm = 0;
  else if (runpar->smforce == 'a') /* always print (regardless adj_str)  */
    runpar->sm = 0;
  else if (runpar->smforce == 'n') /* always no print */
    runpar->sm = 1;
  else if (runpar->smforce != 'l') /* leave as is     */
  {
    printf("smforce should be either p or n or l or a !!! smforce=%c\n",
    runpar->smforce);
    exit(0);
  }
  fprintf(fl.out, "sm=%d\n", sm);

  /* determining IappF */
  if (netpar->nuif == 'i')
  {
    netpar->F.Iapp = netpar->F.IappF;
  }
  if (netpar->nuif == 'n')
  {
    netpar->F.Iapp = netpar->F.nu * netpar->R.Iapp;
  }
  else if (netpar->nuif != 'i')
  {
    printf("nuif=%c should be i or n!\n", netpar->nuif);
    exit(0);
  }

  /* Writing the data */

  fprintf(fl.out, "NEURON\n");

  fprintf(fl.out, "R: beta=%lf theta=%lf Mx=%lf\n", netpar->R.beta,
  netpar->R.theta, netpar->R.Mx);
  fprintf(fl.out, "Iapp=%lf\n", netpar->R.Iapp);

  fprintf(fl.out, "L: beta=%lf theta=%lf Mx=%lf\n", netpar->L.beta,
  netpar->L.theta, netpar->L.Mx);

  fprintf(fl.out, "F: beta=%lf theta=%lf Mx=%lf\n", netpar->F.beta,
  netpar->F.theta, netpar->F.Mx);
  fprintf(fl.out, "Iapp=%lf nu=%lf\n", netpar->F.Iapp, netpar->F.nu);

  fprintf(fl.out, "SYNAPSES\n");

  fprintf(fl.out, "RR: gg=%lf taus=%lf kf=%lf\nUV=%lf taur=%lf tauf=%lf\n",
  netpar->RR.gg, netpar->RR.taus, netpar->RR.kf, netpar->RR.UV,
  netpar->RR.taur, netpar->RR.tauf); 

  fprintf(fl.out, "LR: gg=%lf taus=%lf kf=%lf\nUV=%lf taur=%lf tauf=%lf\n",
  netpar->LR.gg, netpar->LR.taus, netpar->LR.kf, netpar->LR.UV,
  netpar->LR.taur, netpar->LR.tauf); 

  fprintf(fl.out, "RL: gg=%lf taus=%lf kf=%lf\nUV=%lf taur=%lf tauf=%lf\n",
  netpar->RL.gg, netpar->RL.taus, netpar->RL.kf, netpar->RL.UV,
  netpar->RL.taur, netpar->RL.tauf); 

  fprintf(fl.out, "FR: gg=%lf taus=%lf kf=%lf\nUV=%lf taur=%lf tauf=%lf\n",
  netpar->FR.gg, netpar->FR.taus, netpar->FR.kf, netpar->FR.UV,
  netpar->FR.taur, netpar->FR.tauf); 

  fprintf(fl.out, "RF: gg=%lf taus=%lf kf=%lf\nUV=%lf taur=%lf tauf=%lf\n",
  netpar->RF.gg, netpar->RF.taus, netpar->RF.kf, netpar->RF.UV,
  netpar->RF.taur, netpar->RF.tauf); 

  fprintf(fl.out, "FF: gg=%lf taus=%lf kf=%lf\nUV=%lf taur=%lf tauf=%lf\n",
  netpar->FF.gg, netpar->FF.taus, netpar->FF.kf, netpar->FF.UV,
  netpar->FF.taur, netpar->FF.tauf); 

  fprintf(fl.out, "FL: gg=%lf taus=%lf kf=%lf\nUV=%lf taur=%lf tauf=%lf\n",
  netpar->FL.gg, netpar->FL.taus, netpar->FL.kf, netpar->FL.UV,
  netpar->FL.taur, netpar->FL.tauf); 

  fprintf(fl.out, "LF: gg=%lf taus=%lf kf=%lf\nUV=%lf taur=%lf tauf=%lf\n",
  netpar->LF.gg, netpar->LF.taus, netpar->LF.kf, netpar->LF.UV,
  netpar->LF.taur, netpar->LF.tauf); 

  fprintf(fl.out, "tghlin=%c sfact=%c nuif=%c sig=%lf\n", netpar->tghlin,
  netpar->sfact, netpar->nuif, netpar->sig);

  fprintf(fl.out, "GENERAL\n");

  fprintf(fl.out, "deltat=%lf nt=%lf\n", runpar->deltat, runpar->nt);
  fprintf(fl.out, "twrite=%d tmcol=%d ttrans=%d tupdown=%d\n", runpar-> twrite,
  runpar->tmcol, runpar->ttrans, runpar->tupdown);
  fprintf(fl.out, "method=%c incond=%c fpcal=%c smforce=%c\n", runpar->method,
  runpar->incond, runpar->fpcal, runpar->smforce);

  fprintf(fl.out, "\n");
  fflush(fl.out);
}

/* This function creates the structure Varbar */
void create_variables(Varb_ar *Varb, int sm, fl_st fl)
{
  int ieq;

  Varb->neq = 8 * 3;
  Varb->ptr = dvector(1, Varb->neq);

  Varb->RR = Varb->ptr;
  Varb->LR = Varb->ptr + 3;
  Varb->RL = Varb->ptr + 2 * 3;
  Varb->FR = Varb->ptr + 3 * 3;
  Varb->RF = Varb->ptr + 4 * 3;
  Varb->FF = Varb->ptr + 5 * 3;
  Varb->FL = Varb->ptr + 6 * 3;
  Varb->LF = Varb->ptr + 7 * 3;

  for (ieq=1; ieq<=Varb->neq; ieq++) Varb->ptr[ieq] = 0.0;
}

/* This function frees the structure Varbar */
void free_variables(Varb_ar *Varb, int sm, fl_st fl)
{
  free_dvector(Varb->ptr, 1, Varb->neq);
}

/* This function substitutes the initial conditions */
void in_con(Varb_ar *Varbar, net_par *netpar, run_par *runpar,
     int sm, fl_st fl)
{
  int ieq;
  char line[Mline];

  if (runpar->incond == 'r')
  {
    fscanf(fl.tmp, "\n");
    fscanf(fl.tmp, "INITIAL CONDITIONS\n");

    fgets(line, Mline, fl.tmp);
    /* "sRR xRR uRR  sLR xLR uLR  sRL xRL uRL\n" */
    fscanf(fl.tmp, "%lf", &Varbar->ptr[1]);
    for (ieq=2; ieq<=9; ieq++)
      fscanf(fl.tmp, " %lf", &Varbar->ptr[ieq]);
    fscanf(fl.tmp, "\n");

    fgets(line, Mline, fl.tmp);
    /* "sFR xFR uFR  sRF xRF uRF  sFF xFF uFF\n" */
    fscanf(fl.tmp, "%lf", &Varbar->ptr[10]);
    for (ieq=11; ieq<=18; ieq++)
      fscanf(fl.tmp, " %lf", &Varbar->ptr[ieq]);
    fscanf(fl.tmp, "\n");

    fgets(line, Mline, fl.tmp);
    /* "sFL xFL uFL   sLF xLF uLF \n" */
    fscanf(fl.tmp, "%lf", &Varbar->ptr[19]);
    for (ieq=20; ieq<=Varbar->neq; ieq++)
      fscanf(fl.tmp, " %lf", &Varbar->ptr[ieq]);
    fscanf(fl.tmp, "\n");

  }
  else
  {
    printf("Wrong incond:%c\n", runpar->incond);
    exit(1);
  }

  if (!runpar->sm)
  {
    fprintf(fl.out, "Initial conditions\n");
    fprintf(fl.out, " %8.5lf", Varbar->ptr[1]);
    for (ieq=2; ieq<=Varbar->neq; ieq++) 
      fprintf(fl.out, " %8.5lf", Varbar->ptr[ieq]);
    fprintf(fl.out, "\n");
    fprintf(fl.out, "RR: %lf %lf %lf\n", Varbar->RR[1], Varbar->RR[2],
    Varbar->RR[3]);
    fprintf(fl.out, "LR: %lf %lf %lf\n", Varbar->LR[1], Varbar->LR[2],
    Varbar->LR[3]);
    fprintf(fl.out, "RL: %lf %lf %lf\n", Varbar->RL[1], Varbar->RL[2],
    Varbar->RL[3]);

    fprintf(fl.out, "FR: %lf %lf %lf\n", Varbar->FR[1], Varbar->FR[2],
    Varbar->FR[3]);
    fprintf(fl.out, "RF: %lf %lf %lf\n", Varbar->RF[1], Varbar->RF[2],
    Varbar->RF[3]);
    fprintf(fl.out, "FF: %lf %lf %lf\n", Varbar->FF[1], Varbar->FF[2],
    Varbar->FF[3]);

    fprintf(fl.out, "FL: %lf %lf %lf\n", Varbar->FL[1], Varbar->FL[2],
    Varbar->FL[3]);
    fprintf(fl.out, "LF: %lf %lf %lf\n", Varbar->LF[1], Varbar->LF[2],
    Varbar->LF[3]);

    fprintf(fl.out, "\n");
  }
  fflush(fl.out);
}

/* This function solves the differential equations */
void n_run(Varb_ar *Varbar, net_par *netpar, run_par *runpar,
     int sm, avr_val *av, fl_st fl)
{
  Varb_ar k0, k1, k2, k3, k4, Varc;
  stat_run statrun;
  double time;
  double MLnow, t_ML_pos;
  double MR, MRold, ML, MF;
  double t_up, t_down;
  int it, ieq;
  int MLzero;

  create_variables(&k0, sm, fl);
  create_variables(&k1, sm, fl);
  create_variables(&k2, sm, fl);
  create_variables(&k3, sm, fl);
  create_variables(&k4, sm, fl);
  create_variables(&Varc, sm, fl);

  it=0;
  time=0.0;
  MLzero = 0;
  t_ML_pos = -999.9;
  statrun.mud = 1000;
  statrun.tud = dvector(1, statrun.mud);
  statrun.sud = cvector(1, statrun.mud);
  statrun.nud = 0;

  if ((!runpar->sm) && (runpar->tmcol >= runpar->nt))
    pr_fct(Varbar, netpar, runpar, time, it, fl);

  while ((++it) <= runpar->nt)
  {
    time=it*runpar->deltat;
  
    if (it == runpar->ttrans)
    {
      MR = 1000.0 * MRcal(netpar, runpar, Varbar->RR[1], Varbar->RL[1],
              Varbar->RF[1], fl);
      statrun.MRmin = MR;
      statrun.MRmax = MR; 

      ML = 1000.0 * MLcal(netpar, runpar, Varbar->LR[1], Varbar->LF[1], fl);
      statrun.MLmin = ML;
      statrun.MLmax = ML; 

      MF = 1000.0 * MFcal(netpar, runpar, Varbar->FR[1], Varbar->FL[1],
              Varbar->FF[1], fl);
      statrun.MFmin = MF;
      statrun.MFmax = MF; 
    }

    if (runpar->method == 'r')                   /* Runge-Kutta-4 method */
    {
      one_integration_step(netpar, runpar, Varbar, &k0, &k1,
      0.0,                it, time, &Varc, fl);
      one_integration_step(netpar, runpar, Varbar, &k1, &k2, 
      runpar->deltat/2.0, it, time, &Varc, fl);
      one_integration_step(netpar, runpar, Varbar, &k2, &k3,
      runpar->deltat/2.0, it, time, &Varc, fl);
      one_integration_step(netpar, runpar, Varbar, &k3, &k4,
      runpar->deltat,     it, time, &Varc, fl);
      for (ieq=1; ieq<=Varbar->neq; ieq++)
        Varbar->ptr[ieq] = Varbar->ptr[ieq] + (runpar->deltat / 6.0) *
        (k1.ptr[ieq] + 2.0 * (k2.ptr[ieq] + k3.ptr[ieq]) + k4.ptr[ieq]);
    }
    else if (runpar->method =='t')               /* Runge-Kutta-2 method */
    {
      one_integration_step(netpar, runpar, Varbar, &k0, &k1,
      0.0,                it, time, &Varc, fl);
      one_integration_step(netpar, runpar, Varbar, &k1, &k2,
      runpar->deltat/2.0, it, time, &Varc, fl);
      for (ieq=1; ieq<=Varbar->neq; ieq++)
        Varbar->ptr[ieq] = Varbar->ptr[ieq] + runpar->deltat * k2.ptr[ieq];
    }
    else if (runpar->method =='e')  /* Eulear method */
    {
      one_integration_step(netpar, runpar, Varbar, &k0, &k1,
      0.0,                it, time, &Varc, fl);
      for (ieq=1; ieq<=Varbar->neq; ieq++)
        Varbar->ptr[ieq] = Varbar->ptr[ieq] + runpar->deltat * k1.ptr[ieq];
    }
    else
    {
      printf("wrong method!\n");
      exit(0);
    }

    if ((!runpar->sm) && (it >= runpar->nt-runpar->tmcol))
    {
      if ( !(it%runpar->twrite) )
        pr_fct(Varbar, netpar, runpar, time, it, fl);
    }

    if (it >= runpar->ttrans)
    {
      MR = 1000.0 * MRcal(netpar, runpar, Varbar->RR[1], Varbar->RL[1],
           Varbar->RF[1], fl);      
      ML = 1000.0 * MLcal(netpar, runpar, Varbar->LR[1], Varbar->LF[1], fl);
      MF = 1000.0 * MFcal(netpar, runpar, Varbar->FR[1], Varbar->FL[1],
              Varbar->FF[1], fl);
    }

    if ((it >= runpar->ttrans) && (it < runpar->ttrans + runpar->tupdown))
    {
      if (MR < statrun.MRmin)
      {
        statrun.MRmin = MR;
      }
      else if (MR > statrun.MRmax)
      {
        statrun.MRmax = MR; 
      }

      if (ML < statrun.MLmin)
      {
        statrun.MLmin = ML;
      }
      else if (ML > statrun.MLmax)
      {
        statrun.MLmax = ML; 
      }

      if (MF < statrun.MFmin)
      {
        statrun.MFmin = MF;
      }
      else if (MF > statrun.MFmax)
      {
        statrun.MFmax = MF; 
      }
    }

    if (it == runpar->ttrans + runpar->tupdown )
    {
      statrun.MRthr =  (statrun.MRmin + statrun.MRmax) / 2.0;
      fprintf(fl.out, "MRmin=%lf MRmax=%lf MRthr=%lf\n", statrun.MRmin,
      statrun.MRmax, statrun.MRthr);
    }

    if (statrun.MRmax - statrun.MRmin > 1.0e-3 * statrun.MRthr)
    {
      t_up = -999.9;
      t_down = -999.9;

      if (it >= runpar->ttrans + runpar->tupdown)
      {
        if (MRold < statrun.MRthr)
        {
          if (fabs(MR - statrun.MRthr) < runpar->epsilon)
          {
            t_up = time;
          }
          else if (MR > statrun.MRthr)
          {
            t_up = lininter(MRold, MR, statrun.MRthr, time - runpar->deltat,
            time);
          }

          if (t_up > -999.0)
	  {
            fprintf(fl.out, "up: it=%d MRold=%lf MR=%lf "
            "t=%lf\n", it, MRold, MR, t_up);
            statrun.nud++;
            statrun.tud[statrun.nud] = t_up;
            statrun.sud[statrun.nud] = 'u';
	  }
        }
      
        if (MRold > statrun.MRthr)
        {
          if (fabs(MR - statrun.MRthr) < runpar->epsilon)
          {
            t_down = time;
          }
          else if (MR < statrun.MRthr)
          {
            t_down = lininter(MRold, MR, statrun.MRthr, time - runpar->deltat,
            time);
          }

          if (t_down > -999.0)
	  {
            fprintf(fl.out, "down: it=%d MRold=%lf MR=%lf "
            "t=%lf\n", it, MRold, MR, t_down);
            statrun.nud++;
            statrun.tud[statrun.nud] = t_down;
            statrun.sud[statrun.nud] = 'd';
	  }
        }
      }
    }

    if (it >= runpar->ttrans)
    {
      MRold = MR;
    }

    MLnow = 1000.0 * MLcal(netpar, runpar, Varbar->LR[1], Varbar->LF[1], fl);
    if ((MLzero == 0) && (MLnow > runpar->epsilon))
    {
      MLzero = 1;
      t_ML_pos = time;
    }
  }

  av->MRend = 1000.0 * MRcal(netpar, runpar, Varbar->RR[1], Varbar->RL[1],
              Varbar->RF[1], fl);
  av->MLend = 1000.0 * MLcal(netpar, runpar, Varbar->LR[1], Varbar->LF[1], fl);
  av->MFend = 1000.0 * MFcal(netpar, runpar, Varbar->FR[1], Varbar->FL[1],
              Varbar->FF[1], fl);
  av->t_ML_pos = t_ML_pos;

  if (runpar->ttrans < runpar->nt)
  {
    fprintf(fl.out, "MRmin=%lf MRmax=%lf\n", statrun.MRmin, statrun.MRmax);
    av->MRmin = statrun.MRmin;
    av->MRmax = statrun.MRmax;

    fprintf(fl.out, "MLmin=%lf MLmax=%lf\n", statrun.MLmin, statrun.MLmax);
    av->MLmin = statrun.MLmin;
    av->MLmax = statrun.MLmax;

    fprintf(fl.out, "MFmin=%lf MFmax=%lf\n", statrun.MFmin, statrun.MFmax);
    av->MFmin = statrun.MFmin;
    av->MFmax = statrun.MFmax;
  }
  else
  {
    av->MRmin = -999.9;
    av->MRmax = -999.8;

    av->MLmin = -999.9;
    av->MLmax = -999.8;

    av->MFmin = -999.9;
    av->MFmax = -999.8;
  }

  stat_ud_cal(&statrun, netpar, runpar, fl);

  if (runpar->ttrans + runpar->tupdown < runpar->nt)
  {
    av->TT = statrun.TT;
    av->freq = statrun.freq;
    av->Tup = statrun.Tup;
    av->dutyc = statrun.dutyc;
  }

  if (av->TT > 0)
  {
    av->pattern = 'o';
  }
  else
  {
    if (av->MRend < 1.0e-6)
    {
      av->pattern = 'r';  /* zero R */
    }
    else
    {
      if (av->MLend < 1.0e-6)
      {
        if (av->MFend < 1.0e-6)
        {
          av->pattern = 'z';  /* zero L and F */
        }
        else
        {
          av->pattern = 'l';  /* zero L */
        }
      }
      else
      {
        if (av->MFend < 1.0e-6)
        {
          av->pattern = 'f';  /* zero F */
        }
        else
        {
          av->pattern = 'p';  /* positive L and F */
        }
      }
    }
  }

  free_variables(&k0, sm, fl);
  free_variables(&k1, sm, fl);
  free_variables(&k2, sm, fl);
  free_variables(&k3, sm, fl);
  free_variables(&k4, sm, fl);
  free_variables(&Varc, sm, fl);

  free_dvector(statrun.tud, 1, statrun.mud);
  free_cvector(statrun.sud, 1, statrun.mud);
}

/* This function computes one integration step */
void one_integration_step(net_par *netpar, run_par *runpar, Varb_ar *Varbar,
     Varb_ar *kin, Varb_ar *kout, double delt, int it, double time, 
     Varb_ar *Varc, fl_st fl)
{
  double sRR, xRR, uRR;
  double sLR, xLR, uLR;
  double sRL, xRL, uRL;
  double sFR, xFR, uFR;
  double sRF, xRF, uRF;
  double sFF, xFF, uFF;
  double sFL, xFL, uFL;
  double sLF, xLF, uLF;
  double sfactorRR, sfactorRL, sfactorLR;
  double sfactorFR, sfactorRF, sfactorFF;
  double sfactorFL, sfactorLF;
  double MR, ML, MF;
  int ieq;

  Varb_ar *Varo;   /* Variables for one_integration_step */

  if (delt > runpar->epsilon)
  {
    for (ieq=1; ieq<=Varbar->neq; ieq++)
      Varc->ptr[ieq] = Varbar->ptr[ieq] + delt * kin->ptr[ieq];
    Varo = Varc;
  }
  else
  {
    Varo = Varbar;
  }

  /* Computing the RHS of the ODEs */

  sRR = Varo->RR[1];
  sLR = Varo->LR[1];
  sRL = Varo->RL[1];
  sFR = Varo->FR[1];
  sRF = Varo->RF[1];
  sFF = Varo->FF[1];
  sLF = Varo->LF[1];
  sFL = Varo->FL[1];

  MR = MRcal(netpar, runpar, sRR, sRL, sRF, fl);
  ML = MLcal(netpar, runpar, sLR, sLF, fl);
  MF = MFcal(netpar, runpar, sFR, sFL, sFF, fl);

  one_syn_cal(Varo->RR, kout->RR, MR, netpar->RR, time, it, runpar, fl);
  one_syn_cal(Varo->LR, kout->LR, MR, netpar->LR, time, it, runpar, fl);
  one_syn_cal(Varo->RL, kout->RL, ML, netpar->RL, time, it, runpar, fl);
  one_syn_cal(Varo->FR, kout->FR, MR, netpar->FR, time, it, runpar, fl);
  one_syn_cal(Varo->RF, kout->RF, MF, netpar->RF, time, it, runpar, fl);
  one_syn_cal(Varo->FF, kout->FF, MF, netpar->FF, time, it, runpar, fl);
  one_syn_cal(Varo->FL, kout->FL, ML, netpar->FL, time, it, runpar, fl);
  one_syn_cal(Varo->LF, kout->LF, MF, netpar->LF, time, it, runpar, fl);

  /*
  fprintf(fl.col, "RF: %12.8lf %12.8lf %12.8lf\n", Varo->RF[1], Varo->RF[2],
	  Varo->RF[3]);
   */
}

double MRcal(net_par *netpar, run_par *runpar, double sRR, double sRL, 
       double sRF, fl_st fl)
{
  double IallR, MR;

  IallR = netpar->R.Iapp + netpar->RR.gg * sRR - netpar->RL.gg * sRL - 
          netpar->RF.gg * sRF - netpar->R.theta;

  MR = Mcal(IallR, netpar->R.beta, netpar->R.Mx, netpar->sig, netpar->tghlin,
  fl);

  return(MR);
}

double MLcal(net_par *netpar, run_par *runpar, double sLR, double sLF,
       fl_st fl)
{
  double IallL, ML;

  IallL = netpar->LR.gg * sLR - netpar->LF.gg * sLF - netpar->L.theta;
  
  ML = Mcal(IallL, netpar->L.beta, netpar->L.Mx, netpar->sig, netpar->tghlin,
  fl);

  return(ML);
}

double MFcal(net_par *netpar, run_par *runpar, double sFR, double sFL,
       double sFF, fl_st fl)
{
  double IallF, MF;

  IallF = netpar->F.Iapp + netpar->FR.gg * sFR - netpar->FL.gg * sFL - 
          netpar->FF.gg * sFF - netpar->F.theta;

  MF = Mcal(IallF, netpar->F.beta, netpar->F.Mx, netpar->sig, netpar->tghlin,
  fl);

  return(MF);
}

double Mcal(double Iall, double beta, double Mx, double sig, char tghlin,
       fl_st fl)
{
  double MM;

  if (tghlin == 't')
    MM = Mx * tghfunc(beta * Iall / Mx);
  else if (tghlin == 'l')
    MM = linthr(beta * Iall);
  else if (tghlin == 's')
    MM = beta * Iall / (1 + exp(-sig *Iall));
  else
  {
    printf("tghlin=%c should be t or l!\n", tghlin);
    exit(0);
  }

  return(MM);
}

/* This function calculate the synaptic update term for one synapse   */
void one_syn_cal(double *Varo_syn, double *kout_syn, double MM, syn_par synpar,
     double time, int it, run_par *runpar, fl_st fl)
{
  double ssV, xxV, uuV;

  ssV = Varo_syn[1];
  xxV = Varo_syn[2];
  uuV = Varo_syn[3];

  /* sfactor = ( ((netpar->sfact) == ('y')) ? (1.0 - sRR) : (1.0) ); */

  kout_syn[1] =-ssV / synpar.taus +
                synpar.kf * uuV * xxV * MM;

  kout_syn[2] = ( ((synpar.taur) >= runpar->epsilon) ?
		  ((1.0 - xxV) / synpar.taur - uuV * xxV * MM) :
                (0.0) );

  kout_syn[3] = ( ((synpar.tauf) >= runpar->epsilon) ?
                  ((synpar.UV - uuV) / synpar.tauf + 
		    synpar.UV * (1.0 - uuV) * MM) :
                (0.0) );

  /*
  if (it == 1)
  {
    fprintf(fl.col, "ssV=%lf xxV=%lf uuV=%lf\n", ssV, xxV, uuV);
    fprintf(fl.col, "MM=%lf taur=%lf UV=%lf\n", MM, synpar.taur, synpar.UV);
    fprintf(fl.col, "ssK=%lf xxK=%lf uuK=%lf\n", kout_syn[1], kout_syn[2],
    kout_syn[3]);
  }
  */
}

/* This function writes the simulation results on fl.col   */
void pr_fct(Varb_ar *Varbar, net_par *netpar, run_par *runpar, double time,
     int it, fl_st fl)
{
  double MR, ML, MF;
  double VV;
  int ieq, isyn;

  fprintf(fl.col, "%12.5lf", time);

  MR = MRcal(netpar, runpar, Varbar->RR[1], Varbar->RL[1], Varbar->RF[1], fl);
  /* fprintf(fl.col, "\nR: I=%lf M=%lf\n", netpar->R.Iapp, MR); */
  ML = MLcal(netpar, runpar, Varbar->LR[1], Varbar->LF[1], fl);
  MF = MFcal(netpar, runpar, Varbar->FR[1], Varbar->FL[1], Varbar->FF[1], fl);
  /* fprintf(fl.col, "\nF: I=%lf M=%lf\n", netpar->F.Iapp, MF); */


  fprintf(fl.col, "%13.6lf %13.6lf %13.6lf", 1000.0 * MR, 1000.0 * ML, 
  1000.0 * MF);

  fprintf(fl.col, " | ");

  for (isyn=1; isyn<=8; isyn++)
  {
    for (ieq=1; ieq<=3; ieq++)
    {
      VV = *(Varbar->ptr + (isyn-1) * 3 + ieq);
      fprintf(fl.col, " %lf", VV);
    }
    fprintf(fl.col, " | ");
  }

  if (!runpar->sm)
  {
    if (it == runpar->nt)
    {
      fprintf(fl.out, "\nFinal synaptic variables\n");
      for (isyn=1; isyn<=8; isyn++)
      {
        for (ieq=1; ieq<=3; ieq++)
	{
          VV = *(Varbar->ptr + (isyn-1) * 3 + ieq);
          fprintf(fl.out, " %lf", VV);
	}
        fprintf(fl.out, "\n");
      }
    }
  }

  /*
  fprintf(fl.col, "%12.5lf", Varbar->ptr[1]);
  for (ieq=2; ieq<=Varbar->neq; ieq++)
    fprintf(fl.col, "%12.8lf", Varbar->ptr[ieq]);
  */

  fprintf(fl.col, "\n");
}

/* This function computes the statistics of transition from up to down */
/* states and back.                                                    */
void stat_ud_cal(stat_run *statrun, net_par *netpar, run_par *runpar,
     fl_st fl)
{
  double Ta;
  int kud, lud, iud;

  fprintf(fl.out, "\nnud=%d\n", statrun->nud);
  for (statrun->iud=1; statrun->iud <= statrun->nud; statrun->iud++)
  {
    fprintf(fl.out, "iud=%d sud=%c tud=%lf\n", statrun->iud,
    statrun->sud[statrun->iud], statrun->tud[statrun->iud]);
  }

  /* Computing TT and freq */

  lud = statrun->nud / 2;
  fprintf(fl.out, "lud=%d\n", lud);

  kud = 2 * ((statrun->nud - 1) / 2) + 1;
  fprintf(fl.out, "\nkud=%d\n", kud);

  if (kud <= 2)
  {
    statrun->TT = -999,9;
    statrun->freq = -9.999;
  }
  else 
  {
    statrun->TT = 2.0 * (statrun->tud[kud] - statrun->tud[1]) / (kud - 1);
    statrun->freq = 1000.0 / statrun->TT;
  }
  fprintf(fl.out, "TT=%lf freq=%lf\n", statrun->TT, statrun->freq);


  if (kud <= 3)
  {
    statrun->Tup = -99.99;
    statrun->dutyc = -9.99;
  }
  else
  {
    Ta = 0.0;
    for (iud=1; iud<=lud; iud++)
    {
      Ta += statrun->tud[2*iud] - statrun->tud[2*iud-1];
    }
    Ta /= lud;

    if (statrun->sud[1] == 'u')
    {
      statrun->Tup = Ta;
      statrun->dutyc = statrun->Tup / statrun->TT;
    }
    else if (statrun->sud[1] == 'd')
    {
      statrun->Tup = statrun->TT - Ta;
      statrun->dutyc = statrun->Tup / statrun->TT;
    }
    else
    {
      printf("sud should be either u or d!\n");
      exit(0);
    }
  }

  fprintf(fl.out, "Ta=%lf dutyc=%lf\n", statrun->Tup, statrun->dutyc);
}

/* Linear interpolation */
double lininter(double x1, double x2, double xc, double y1, double y2)
{
  double linter;

  linter = ((xc-x1)*y2+(x2-xc)*y1) / (x2-x1) ;
  return(linter);
}
