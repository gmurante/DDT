/*
* @file
* This file is part of the developer version of GADGET3 and contains
* the license conditions for its usage.
*
* @author GADGET-development team, led by Volker Springel and Klaus Dolag.
*
* @section LICENSE
* Copyright (c) 2016, Volker Springel, Klaus Dolag, and all contributing authors
* (see change logs). All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Received source code may be modified and used as convenient.
*
* 2. Redistributions of source code or in binary form is only possible with
*    explicit agreement of the copyright holders.
*
* 3. Redistributions of source code must retain the above copyright notice,
*    this list of conditions and the following disclaimer.
*
* 4. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
*
* 5. Neither the name of the copyright holder nor the names of its
*    contributors may be used to endorse or promote products derived from this
*    software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

#include "allvars.h"
#include "proto.h"


/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */

/*! This routine contains the main simulation loop that iterates over
 * single timesteps. The loop terminates when the cpu-time limit is
 * reached, when a `stop' file is found in the output directory, or
 * when the simulation ends because we arrived at TimeMax.
 */

void run()
{
  char bfr[10000];
  int n=10000, i ,j;
  FILE *f;


  if(ThisTask==0)
    {
      printf("\n\n Domain Decomposition and Tree Building DONE\n"); fflush(stdout);
    }

  CPU_Step[CPU_MISC] += measure_time();
  

  /* writing off domain decomposition */

  sprintf(bfr,"domaindecomp%03d.dat",ThisTask);
  f=fopen(bfr,"w");
  for(i=0; i<NumPart; i++)
    fprintf(f,"%f %f %f\n",P[i].Pos[0],P[i].Pos[1],P[i].Pos[2]);
  fclose(f);

  if(ThisTask==0)
    printf(" Finding neighbours...\n");

  /* printing neighbours of n random particles */
  for(i=0; i<n; i++)
    {
      do {
	j = rand()/((float)RAND_MAX)*NumPart;
      } while( P[i].Type!=1 && P[i].Type!=0);
      look_around(j, 50.0);
	
    }

  write_cpu_log(); 


  return;
}




/*! This routine writes for every synchronisation point in the timeline information to two log-files:
 * In FdInfo, we just list the timesteps that have been done, while in
 * FdTimebins we inform about the distribution of particles over the timebins, and which timebins are active on this step.
 * code is stored.
 */
void output_log_messages(void)
{
  double z;
  int i;
  long long tot_count[TIMEBINS];
  long long tot_count_sph[TIMEBINS];
  long long tot_cumulative[TIMEBINS];

  sumup_large_ints(TIMEBINS, TimeBinCount, tot_count);
  sumup_large_ints(TIMEBINS, TimeBinCountSph, tot_count_sph);

  if(ThisTask == 0)
    {
      if(All.ComovingIntegrationOn)
	{
	  z = 1.0 / (All.Time) - 1;
	  fprintf(FdInfo, "\nSync-Point %d, Time: %g, Redshift: %g, Nf = %d%09d, Systemstep: %g, Dloga: %g\n",
		  All.NumCurrentTiStep, All.Time, z,
		  (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
		  All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
	  printf("\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep,
		 All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
	  fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n",
		  All.NumCurrentTiStep, All.Time, z, All.TimeStep,
		  log(All.Time) - log(All.Time - All.TimeStep));
	  fflush(FdInfo);
	}
      else
	{
	  fprintf(FdInfo, "\nSync-Point %d, Time: %g, Nf = %d%09d, Systemstep: %g\n", All.NumCurrentTiStep,
		  All.Time, (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000),
		  All.TimeStep);
	  printf("\nSync-Point %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time, All.TimeStep);
	  fprintf(FdTimebin, "\nSync-Point %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time,
		  All.TimeStep);
	  fflush(FdInfo);
	}

      for(i = 1, tot_cumulative[0] = tot_count[0]; i < TIMEBINS; i++)
	tot_cumulative[i] = tot_count[i] + tot_cumulative[i - 1];

      /* GM: timebins.txt not currently enabled

      for(i = 0; i < TIMEBINS; i++)
	{
	  for(j = 0, sum = 0; j < All.CPU_TimeBinCountMeasurements[i]; j++)
	    sum += All.CPU_TimeBinMeasurements[i][j];
	  if(All.CPU_TimeBinCountMeasurements[i])
	    avg_CPU_TimeBin[i] = sum / All.CPU_TimeBinCountMeasurements[i];
	  else
	    avg_CPU_TimeBin[i] = 0;
	}

      for(i = All.HighestOccupiedTimeBin, weight = 1, sum = 0; i >= 0 && tot_count[i] > 0; i--, weight *= 2)
	{
	  if(weight > 1)
	    corr_weight = weight / 2;
	  else
	    corr_weight = weight;

	  frac_CPU_TimeBin[i] = corr_weight * avg_CPU_TimeBin[i];
	  sum += frac_CPU_TimeBin[i];
	}

      for(i = All.HighestOccupiedTimeBin; i >= 0 && tot_count[i] > 0; i--)
	{
	  if(sum)
	    frac_CPU_TimeBin[i] /= sum;
	}


      printf
	("Occupied timebins: non-cells     cells       dt                 cumulative A D    avg-time  cpu-frac\n");
      fprintf(FdTimebin,
	      "Occupied timebins: non-cells     cells       dt                 cumulative A D    avg-time  cpu-frac\n");
      for(i = TIMEBINS - 1, tot = tot_sph = 0; i >= 0; i--)
	if(tot_count_sph[i] > 0 || tot_count[i] > 0)
	  {
	    printf(" %c  bin=%2d      %10llu  %10llu   %16.12f       %10llu %c %c  %10.2f    %5.1f%%\n",
		   TimeBinActive[i] ? 'X' : ' ',
		   i, tot_count[i] - tot_count_sph[i], tot_count_sph[i],
		   i > 0 ? (((integertime) 1) << i) * All.Timebase_interval : 0.0, tot_cumulative[i],
		   (i == All.HighestActiveTimeBin) ? '<' : ' ',
		   (tot_cumulative[i] > All.TreeDomainUpdateFrequency * All.TotNumPart) ? '*' : ' ',
		   avg_CPU_TimeBin[i], 100.0 * frac_CPU_TimeBin[i]);

	    fprintf(FdTimebin,
		    " %c  bin=%2d      %10llu  %10llu   %16.12f       %10llu %c %c  %10.2f    %5.1f%%\n",
		    TimeBinActive[i] ? 'X' : ' ', i, tot_count[i] - tot_count_sph[i], tot_count_sph[i],
		    i > 0 ? (((integertime) 1) << i) * All.Timebase_interval : 0.0, tot_cumulative[i],
		    (i == All.HighestActiveTimeBin) ? '<' : ' ',
		    (tot_cumulative[i] > All.TreeDomainUpdateFrequency * All.TotNumPart) ? '*' : ' ',
		    avg_CPU_TimeBin[i], 100.0 * frac_CPU_TimeBin[i]);

	    if(TimeBinActive[i])
	      {
		tot += tot_count[i];
		tot_sph += tot_count_sph[i];
	      }
	  }
      printf("               ------------------------\n");
      fprintf(FdTimebin, "               ------------------------\n");

	{
	  printf("Total active:   %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);
	  fprintf(FdTimebin, "Total active:   %10llu  %10llu    Sum: %10llu\n", tot - tot_sph, tot_sph, tot);

	}
      fprintf(FdTimebin, "\n");
      fflush(FdTimebin);
      */
    }
}


void write_cpu_log(void)
{
  double max_CPU_Step[CPU_PARTS], avg_CPU_Step[CPU_PARTS], t0, t1, tsum;
  int i;

  CPU_Step[CPU_MISC] += measure_time();

  for(i = 1, CPU_Step[0] = 0; i < CPU_PARTS4SUM; i++)
    CPU_Step[0] += CPU_Step[i];

  MPI_Reduce(CPU_Step, max_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_MAX, 0, MYMPI_COMM_WORLD);
  MPI_Reduce(CPU_Step, avg_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_SUM, 0, MYMPI_COMM_WORLD);


  if(ThisTask == 0)
    {
      for(i = 0; i < CPU_PARTS; i++)
	avg_CPU_Step[i] /= NTask;

      put_symbol(0.0, 1.0, '#');

      for(i = 1, tsum = 0.0; i < CPU_PARTS4SUM; i++)
	{
	  if(max_CPU_Step[i] > 0)
	    {
	      t0 = tsum;
	      t1 = tsum + avg_CPU_Step[i] * (avg_CPU_Step[i] / max_CPU_Step[i]);
	      put_symbol(t0 / avg_CPU_Step[0], t1 / avg_CPU_Step[0], CPU_Symbol[i]);
	      tsum += t1 - t0;

	      t0 = tsum;
	      t1 = tsum + avg_CPU_Step[i] * ((max_CPU_Step[i] - avg_CPU_Step[i]) / max_CPU_Step[i]);
	      put_symbol(t0 / avg_CPU_Step[0], t1 / avg_CPU_Step[0], CPU_SymbolImbalance[i]);
	      tsum += t1 - t0;
	    }
	}

      put_symbol(tsum / max_CPU_Step[0], 1.0, '-');

      /* GM: FdBalance not currently enabled in begrun.c
      fprintf(FdBalance, "Step=%7d  sec=%10.3f  Nf=%2d%09d  %s\n", All.NumCurrentTiStep, max_CPU_Step[0],
	      (int) (GlobNumForceUpdate / 1000000000), (int) (GlobNumForceUpdate % 1000000000), CPU_String);
      fflush(FdBalance);
      */

      if(All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin] == NUMBER_OF_MEASUREMENTS_TO_RECORD)
	{
	  All.CPU_TimeBinCountMeasurements[All.HighestActiveTimeBin]--;
	  memmove(&All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][0],
		  &All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][1],
		  (NUMBER_OF_MEASUREMENTS_TO_RECORD - 1) * sizeof(double));
	}

      All.CPU_TimeBinMeasurements[All.HighestActiveTimeBin][All.CPU_TimeBinCountMeasurements
							    [All.HighestActiveTimeBin]++] = max_CPU_Step[0];
    }

  CPUThisRun += CPU_Step[0];

  for(i = 0; i < CPU_PARTS; i++)
    CPU_Step[i] = 0;

  if(ThisTask == 0)
    {
      for(i = 0; i < CPU_PARTS; i++)
	All.CPU_Sum[i] += avg_CPU_Step[i];

      fprintf(FdCPU, "Step %d, Time: %g, CPUs: %d\n", All.NumCurrentTiStep, All.Time, NTask);
      fprintf(FdCPU,
	      "total                      %10.2f  %5.1f%%\n"
	      "domain                     %10.2f  %5.1f%%\n"
	      "peano                      %10.2f  %5.1f%%\n" 
	      "treebuild                  %10.2f  %5.1f%%\n"
	      "treewalk                   %10.2f  %5.1f%%\n"
	      "i/o                        %10.2f  %5.1f%%\n"
	      "misc                       %10.2f  %5.1f%%\n",	      
	      All.CPU_Sum[CPU_ALL], 100.0,

	      All.CPU_Sum[CPU_DOMAIN], All.CPU_Sum[CPU_DOMAIN]/All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_PEANO], All.CPU_Sum[CPU_PEANO]/All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEBUILD], All.CPU_Sum[CPU_TREEBUILD]/All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEWALK], All.CPU_Sum[CPU_TREEWALK]/All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SNAPSHOT], All.CPU_Sum[CPU_SNAPSHOT]/All.CPU_Sum[CPU_ALL] * 100,

	      All.CPU_Sum[CPU_MISC], (All.CPU_Sum[CPU_MISC]) / All.CPU_Sum[CPU_ALL] * 100);
      fprintf(FdCPU, "\n");
      fflush(FdCPU);


      printf( "\n\n Timings: \n");
      printf( "total                      %10.2f  %5.1f%%\n"
	      "domain                     %10.2f  %5.1f%%\n"
	      "peano                      %10.2f  %5.1f%%\n" 
	      "treebuild                  %10.2f  %5.1f%%\n"
	      "treewalk                   %10.2f  %5.1f%%\n"
	      "i/o                        %10.2f  %5.1f%%\n\n"
	      "ngbcomm                    %10.2f  %5.1f%%\n"
	      "ngbwait                    %10.2f  %5.1f%%\n"
	      "ngbcomput                  %10.2f  %5.1f%%\n\n"	     
	      "misc                       %10.2f  %5.1f%%\n",	      
	      All.CPU_Sum[CPU_ALL], 100.0,

	      All.CPU_Sum[CPU_DOMAIN], All.CPU_Sum[CPU_DOMAIN]/All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_PEANO], All.CPU_Sum[CPU_PEANO]/All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEBUILD], All.CPU_Sum[CPU_TREEBUILD]/All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TREEWALK], All.CPU_Sum[CPU_TREEWALK]/All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_SNAPSHOT], All.CPU_Sum[CPU_SNAPSHOT]/All.CPU_Sum[CPU_ALL] * 100,

	      All.CPU_Sum[CPU_COMM], All.CPU_Sum[CPU_COMM]/All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_COMP], All.CPU_Sum[CPU_COMP]/All.CPU_Sum[CPU_ALL] * 100,
	      All.CPU_Sum[CPU_TIMEWAIT], All.CPU_Sum[CPU_TIMEWAIT]/All.CPU_Sum[CPU_ALL] * 100,



	      All.CPU_Sum[CPU_MISC], (All.CPU_Sum[CPU_MISC]) / All.CPU_Sum[CPU_ALL] * 100);
      printf("\n");
      fflush(stdout);

    }
}


void put_symbol(double t0, double t1, char c)
{
  int i, j;

  i = (int) (t0 * CPU_STRING_LEN + 0.5);
  j = (int) (t1 * CPU_STRING_LEN);

  if(i < 0)
    i = 0;
  if(j < 0)
    j = 0;
  if(i >= CPU_STRING_LEN)
    i = CPU_STRING_LEN;
  if(j >= CPU_STRING_LEN)
    j = CPU_STRING_LEN;

  while(i <= j)
    CPU_String[i++] = c;

  CPU_String[CPU_STRING_LEN] = 0;
}




void check_particles_info(const char *func, const char *file, int linenr)
{
  int i, k, vok = 0, pok = 0, vsph = 0, mok = 0, vstar = 0;
  double vv;

  MPI_Barrier(MYMPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("Checking particle data (function %s in file %s at line %d) ...\n", func, file, linenr);

  for(i = 0; i < NumPart; i++)
    {

      for(k = 0; k < 3; k++)
	{
	  if(P[i].Vel[k] > -1e8 && P[i].Vel[k] < 1e8)
	    {
	      vv = sqrt(P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
	      if(vv > 15000)
		{
		  printf
		    ("task=%d: WARNING: Large velocity for particle %d ID %llu v[%d]=%g, renormalizing it !!\n",
		     ThisTask, i, (unsigned long long) P[i].ID, k, vv);
		  fflush(stdout);
		  P[i].Vel[0] = P[i].Vel[0] / vv * 10000;
		  P[i].Vel[1] = P[i].Vel[1] / vv * 10000;
		  P[i].Vel[2] = P[i].Vel[2] / vv * 10000;
		}
	      vok++;
	    }
	  else
	    {
	      printf
		("task=%d:  strange value in velocity in for particle %d ID %llu , type=%d, mass=%g, v[%d]=%g\n",
		 ThisTask, i, (unsigned long long) P[i].ID, P[i].Type, P[i].Mass, k, P[i].Vel[k]);
	      if(P[i].Type ==0)
		printf("        vred[%d]=%g  a_hydro[%d]=%g\n",k,SphP[i].VelPred[k],k,SphP[i].a.HydroAccel[k]);
	      fflush(stdout);
	      endrun(712401);
	    }

	  if(P[i].Pos[k] > -10000 && P[i].Pos[k] < All.BoxSize + 10000)
	    pok++;
	  else
	    {
	      printf("task=%d:  strange value in position in for particle %d ID %llu x[%d]=%g\n", ThisTask, i,
		     (unsigned long long) P[i].ID, k, P[i].Pos[k]);
	      fflush(stdout);
	      endrun(712402);
	    }
	}

      double massDMpart;
      
      if(All.MassTable[1] > 0)
	massDMpart = All.MassTable[1];
      else
        if(All.MassTable[0] > 0)
          massDMpart = All.MassTable[0]*10;
        else
          massDMpart = P[0].Mass*10;


      if((P[i].Mass > massDMpart/5000 && P[i].Mass < massDMpart*10) || P[i].Type == 2 || P[i].Type == 3 || P[i].Type == 5)
	mok++;
      else
	{
	  printf("task=%d:  strange value in Mass in for particle %d ID %llu Type %d mass=%g\n", ThisTask, i,
		 (unsigned long long) P[i].ID, P[i].Type, P[i].Mass);
	  fflush(stdout);
	  if(P[i].Mass > 0)
	    endrun(712403);
	}

      if(P[i].Type == 0)
	{
	  if(SphP[i].d.Density >= 0 && SphP[i].d.Density < 1e10)
	    vsph++;
	  else
	    printf("task=%d: Gas Particle id=%llu,m=%e strange Density value: %e\n",
		   ThisTask, (unsigned long long) P[i].ID, P[i].Mass, SphP[i].d.Density);

	  if((SphP[i].EntropyPred > -1e20 && SphP[i].EntropyPred < 1e20))
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu, m=%e strange predicted Entropy value in hydro: %e\n",
		   ThisTask, (unsigned long long) P[i].ID, P[i].Mass, SphP[i].EntropyPred);

	  if((SphP[i].Entropy > -1e20 && SphP[i].Entropy < 1e20))
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu, m=%e strange Entropy value in hydro: %e\n",
		   ThisTask, (unsigned long long) P[i].ID, P[i].Mass, SphP[i].Entropy);

	  if((SphP[i].e.DtEntropy > -1e20 && SphP[i].e.DtEntropy < 1e20))
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu, m=%e strange DtEntropy value in hydro: %e\n",
		   ThisTask, (unsigned long long) P[i].ID, P[i].Mass, SphP[i].e.DtEntropy);

	  if((SphP[i].Pressure > -1e20 && SphP[i].Pressure < 1e20))
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu,m=%e strange Pressure value in hydro: %e,%e\n",
		   ThisTask, (unsigned long long) P[i].ID, P[i].Mass, SphP[i].EntropyPred, SphP[i].Pressure);

	  if(SphP[i].a.HydroAccel[0] > -1e20 && SphP[i].a.HydroAccel[0] < 1e20 &&
	     SphP[i].a.HydroAccel[1] > -1e20 && SphP[i].a.HydroAccel[1] < 1e20 &&
	     SphP[i].a.HydroAccel[2] > -1e20 && SphP[i].a.HydroAccel[2] < 1e20)
	    vsph++;
	  else
	    printf("task=%d: Particle id=%llu,m=%e strange acceleration value in hydro: %e,%e,%e\n",
		   ThisTask, (unsigned long long) P[i].ID, P[i].Mass,
		   SphP[i].a.HydroAccel[0], SphP[i].a.HydroAccel[1], SphP[i].a.HydroAccel[2]);
	}

      if(P[i].Type > 5 || P[i].Type < 0)
	{
	  printf("task=%d:  P[i=%d].Type=%d\n", ThisTask, i, P[i].Type);
	  endrun(712411);
	}


    }

  if(ThisTask == 0)
    printf
      ("Positions, Velocities, Masses and  SphVelocities fine for (%d,%d,%d,%d) of %d/%d/%d cases on task 0 (and %d stars)...\n\n",
       pok, vok, mok, vsph / 6, NumPart * 3, NumPart, N_gas, vstar);

  MPI_Barrier(MYMPI_COMM_WORLD);
}
