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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>

#include "allvars.h"
#include "proto.h"



/*! \file init.c
 *  \brief code for initialisation of a simulation from initial conditions
 */


/*! This function reads the initial conditions, and allocates storage for the
 *  tree(s). Various variables of the particle data are initialised and An
 *  intial domain decomposition is performed. If SPH particles are present,
 *  the inial SPH smoothing lengths are determined.
 */

/* GM: stripping violently down this one! 
   only leaving the parsing of parameter file and essential
   initializations */


void init(void)
{
  int i, j;
  double t0_ics, t1_ics;


  All.Time = All.TimeBegin;
  /*
  set_cosmo_factors_for_current_time();
  */

  /* GM: change this, no restart>0 allowed now */

  if(RestartFlag >0)
    {
      if(ThisTask == 0)
	printf("This mini-app does not allow checkpointing/restarting\n or other advanced functionalities\n");
      endrun(0);
    }
  if(All.ComovingIntegrationOn==1)
    {
      if(ThisTask == 0)
	printf("This mini-app does not allow the use of comoving integration\n");
      endrun(0);
    }


  t0_ics = second();
  switch (All.ICFormat)
    {
    case 1:
    case 2:
      //    case 3:
      //    case 4:
      /* GM: unused at the moment
      if(RestartFlag >= 2 && RestartSnapNum >= 0)
	{

	  char fname[1000];

	  if(All.NumFilesPerSnapshot > 1)
	    sprintf(fname, "%s/snapdir_%03d/%s_%03d", All.OutputDir, RestartSnapNum, All.SnapshotFileBase,
		    RestartSnapNum);
	  else
	    sprintf(fname, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, RestartSnapNum);
	  read_ic(fname);

	}
      else
	{ */
	  read_ic(All.InitCondFile);
	  /*	} */
      break;

    default:
    case 3:
    case 4:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }
  t1_ics = second();

  if(ThisTask == 0)
    printf("reading ICs took %g sec\n", timediff(t0_ics, t1_ics));

  All.Time = All.TimeBegin;
  /*
  set_cosmo_factors_for_current_time();
  */



  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      All.Ti_Current = 0;
      //      a3 = All.Time * All.Time * All.Time;
      //      atime = All.Time;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      All.Ti_Current = 0;
      //      a3 = 1;
      //      atime = 1;
    }


  set_softenings();

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;
  if(RestartFlag == 2)
    {
      if(RestartSnapNum < 0)
	{
	  char *underscore = strrchr(All.InitCondFile, '_');
	  if(!underscore)
	    {
	      char buf[1000];
	      sprintf(buf, "Your input file '%s' lacks an underscore. Cannot infer next snapshot number.\n",
		      All.InitCondFile);
	      terminate(buf);
	    }
	  else
	    {
	      All.SnapshotFileCount = atoi(underscore + 1) + 1;
	    }
	}
      else
	{
	  All.SnapshotFileCount = RestartSnapNum + 1;
	  All.SnapshotFileCount--;
	}
    }

  All.TotNumOfForces = 0;

  All.TopNodeAllocFactor = 0.008;
  All.TreeAllocFactor = 0.25;

  /*
  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();
  */

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;

  for(i = 0; i < GRAVCOSTLEVELS; i++)
    All.LevelToTimeBin[i] = 0;

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < GRAVCOSTLEVELS; j++)
      P[i].GravCost[j] = 0;

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
    }


  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].g.GravAccel[j] = 0;
	}
      P[i].Ti_begstep = 0;
      P[i].Ti_current = 0;
      P[i].TimeBin = 0;

      if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS)
	P[i].OldAcc = 0;	/* Do not zero in 2lpt case as masses are stored here */



    }

  for(i = 0; i < TIMEBINS; i++)
    TimeBinActive[i] = 1;

  /*
  reconstruct_timebins();
  */


  /* GM: I comment out SPH initialization for the moment */
  /*  for(i = 0; i < N_gas; i++)	/initialize sph_properties 
    {
      SphP[i].EntropyPred = SphP[i].Entropy;

      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].a.HydroAccel[j] = 0;
	}

      SphP[i].e.DtEntropy = 0;

      if(RestartFlag == 0)
	{
	  PPP[i].Hsml = 0;

	  SphP[i].d.Density = -1;
	  SphP[i].v.DivVel = 0;
	}

    }
*/

  Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */

  TreeReconstructFlag = 1;


  domain_Decomposition(0, 0);	/* do initial domain decomposition (gives equal numbers of particles) */

  set_softenings();

  /* will build tree */
  ngb_treebuild();

  All.Ti_Current = 0;

  /*
  if(RestartFlag != 3 && RestartFlag != 5)
    setup_smoothinglengths();
  */



  /* at this point, the entropy variable actually contains the
   * internal energy, read in from the initial conditions file.
   * Once the density has been computed, we can convert to entropy.
   */

  /* GM: I comment out SPH initialization for the moment */
  /*for(i = 0; i < N_gas; i++)	/initialize sph_properties 
    {
      if(header.flag_entropy_instead_u == 0)
	{
	  if(ThisTask == 0 && i == 0)
	    printf("Converting u -> entropy !\n");

	  SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].d.Density / a3, GAMMA_MINUS1);
	  SphP[i].EntropyPred = SphP[i].Entropy;
	}

      SphP[i].e.DtEntropy = 0;

      SphP[i].v.DivVel = 0;


    }
  */

}



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
/* GM: commented out for the moment (will be used with sph) 
void setup_smoothinglengths(void)
{
  int i, no, p;

  if(RestartFlag == 0)
    {
      for(i = 0; i < N_gas; i++)
	{
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }

	  PPP[i].Hsml =
	    pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
	  if(All.SofteningTable[0] != 0 && PPP[i].Hsml > 200.0 * All.SofteningTable[0])
	    PPP[i].Hsml = All.SofteningTable[0];
	}


    }


  density();

}
*/


void test_id_uniqueness(void)
{
  int i;
  double t0, t1;
  MyIDType *ids, *ids_first;

  if(ThisTask == 0)
    {
      printf("Testing ID uniqueness...\n");
      fflush(stdout);
    }

  if(NumPart == 0)
    {
      printf("need at least one particle per cpu\n");
      endrun(8);
    }

  t0 = second();

  ids = (MyIDType *) mymalloc("ids", NumPart * sizeof(MyIDType));
  ids_first = (MyIDType *) mymalloc("ids_first", NTask * sizeof(MyIDType));

  for(i = 0; i < NumPart; i++)
    ids[i] = P[i].ID;

  parallel_sort(ids, NumPart, sizeof(MyIDType), compare_IDs);

  for(i = 1; i < NumPart; i++)
    if(ids[i] == ids[i - 1])
      {
#ifdef LONGIDS
	printf("non-unique ID=%d%09d found on task=%d (i=%d NumPart=%d)\n",
	       (int) (ids[i] / 1000000000), (int) (ids[i] % 1000000000), ThisTask, i, NumPart);

#else
	printf("non-unique ID=%d found on task=%d   (i=%d NumPart=%d)\n", (int) ids[i], ThisTask, i, NumPart);
#endif
	endrun(12);
      }

  MPI_Allgather(&ids[0], sizeof(MyIDType), MPI_BYTE, ids_first, sizeof(MyIDType), MPI_BYTE, MYMPI_COMM_WORLD);

  if(ThisTask < NTask - 1)
    if(ids[NumPart - 1] == ids_first[ThisTask + 1])
      {
	printf("non-unique ID=%d found on task=%d\n", (int) ids[NumPart - 1], ThisTask);
	endrun(13);
      }

  myfree(ids_first);
  myfree(ids);

  t1 = second();

  if(ThisTask == 0)
    {
      printf("success.  took=%g sec\n", timediff(t0, t1));
      fflush(stdout);
    }
}

int compare_IDs(const void *a, const void *b)
{
  if(*((MyIDType *) a) < *((MyIDType *) b))
    return -1;

  if(*((MyIDType *) a) > *((MyIDType *) b))
    return +1;

  return 0;
}


/* GM: this was in gravtree.c 
   note that comoving integration is currently forbidden */
/*! This function sets the (comoving) softening length of all particle
 *  types in the table All.SofteningTable[...].  We check that the physical
 *  softening length is bounded by the Softening-MaxPhys values.
 */
void set_softenings(void)
{
  int i;

  if(All.ComovingIntegrationOn)
    {
      if(All.SofteningGas * All.Time > All.SofteningGasMaxPhys)
	All.SofteningTable[0] = All.SofteningGasMaxPhys / All.Time;
      else
	All.SofteningTable[0] = All.SofteningGas;

      if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys)
	All.SofteningTable[1] = All.SofteningHaloMaxPhys / All.Time;
      else
	All.SofteningTable[1] = All.SofteningHalo;

      if(All.SofteningDisk * All.Time > All.SofteningDiskMaxPhys)
	All.SofteningTable[2] = All.SofteningDiskMaxPhys / All.Time;
      else
	All.SofteningTable[2] = All.SofteningDisk;

      if(All.SofteningBulge * All.Time > All.SofteningBulgeMaxPhys)
	All.SofteningTable[3] = All.SofteningBulgeMaxPhys / All.Time;
      else
	All.SofteningTable[3] = All.SofteningBulge;

      if(All.SofteningStars * All.Time > All.SofteningStarsMaxPhys)
	All.SofteningTable[4] = All.SofteningStarsMaxPhys / All.Time;
      else
	All.SofteningTable[4] = All.SofteningStars;

      if(All.SofteningBndry * All.Time > All.SofteningBndryMaxPhys)
	All.SofteningTable[5] = All.SofteningBndryMaxPhys / All.Time;
      else
	All.SofteningTable[5] = All.SofteningBndry;
    }
  else
    {
      All.SofteningTable[0] = All.SofteningGas;
      All.SofteningTable[1] = All.SofteningHalo;
      All.SofteningTable[2] = All.SofteningDisk;
      All.SofteningTable[3] = All.SofteningBulge;
      All.SofteningTable[4] = All.SofteningStars;
      All.SofteningTable[5] = All.SofteningBndry;
    }

  for(i = 0; i < 6; i++)
    All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

  /* GM: not used now
  All.MinGasHsml = All.MinGasHsmlFractional * All.ForceSoftening[0];
  */
}
