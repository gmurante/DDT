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
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include "allvars.h"
#include "utilities.h"
#include "proto.h"


/*! \file io.c
 *  \brief Output of a snapshot file to disk.
 */

static int n_type[6];
static long long ntot_type_all[6];

#ifdef LT_SEvDbg
double metals[2 * LT_NMetP + 3], tot_metals[2 * LT_NMetP + 3];
#endif

#ifdef LT_SEv_INFO_DETAILS
double DetailsWSum[LT_NMetP], DetailsWoSum[LT_NMetP];
double my_weight_sum_details, My_weight_sum_details;
unsigned int *my_weight_sum_details_N;
#endif

#define min(x, y) (((x) < (y)) ? (x) : (y))

/* GM: the I/O has to be COMPLETELY rewritten */


static int n_info;

/*! This function writes a snapshot of the particle distribution to one or
 * several files using Gadget's default file format.  If
 * NumFilesPerSnapshot>1, the snapshot is distributed into several files,
 * which are written simultaneously. Each file contains data from a group of
 * processors of size roughly NTask/NumFilesPerSnapshot.
 */
void savepositions(int num)
{
  size_t bytes;
  char buf[500];
  int n, filenr, gr, ngroups, masterTask, lastTask;

  CPU_Step[CPU_MISC] += measure_time();



  if(DumpFlag == 1)
    {
      if(ThisTask == 0)
	printf("\nwriting snapshot file #%d... \n", num);


      if(!(CommBuffer = mymalloc("CommBuffer", bytes = All.BufferSize * 1024 * 1024)))
	{
	  printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(2);
	}


      if(NTask < All.NumFilesPerSnapshot)
	{
	  if(ThisTask == 0)
	    printf
	      ("Fatal error.\nNumber of processors must be larger or equal than All.NumFilesPerSnapshot.\n");
	  endrun(0);
	}
      if(All.SnapFormat < 1 || All.SnapFormat > 3)
	{
	  if(ThisTask == 0)
	    printf("Unsupported File-Format\n");
	  endrun(0);
	}
      if(All.SnapFormat == 3)
	{
	  if(ThisTask == 0)
	    printf("Code wasn't compiled with HDF5 support enabled!\n");
	  endrun(0);
	}


      /* determine global and local particle numbers */
      for(n = 0; n < 6; n++)
	n_type[n] = 0;

      for(n = 0; n < NumPart; n++)
	n_type[P[n].Type]++;

      sumup_large_ints(6, n_type, ntot_type_all);



      /* assign processors to output files */
      distribute_file(All.NumFilesPerSnapshot, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

      if(All.NumFilesPerSnapshot > 1)
	{
	  if(ThisTask == 0)
	    {
	      sprintf(buf, "%s/snapdir_%03d", All.OutputDir, num);
	      mkdir(buf, 02755);
	    }
	  MPI_Barrier(MYMPI_COMM_WORLD);
	}

      if(All.NumFilesPerSnapshot > 1)
	sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num, filenr);
      else
	sprintf(buf, "%s%s_%03d", All.OutputDir, All.SnapshotFileBase, num);


      ngroups = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
      if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel))
	ngroups++;

      for(gr = 0; gr < ngroups; gr++)
	{
	  if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
	    {
	      write_file(buf, masterTask, lastTask);

	    }
	  MPI_Barrier(MYMPI_COMM_WORLD);
	}

      myfree(CommBuffer);


      if(ThisTask == 0)
	printf("done with snapshot.\n");

      All.Ti_lastoutput = All.Ti_Current;

      CPU_Step[CPU_SNAPSHOT] += measure_time();

    }



}



/*! This function fills the write buffer with particle data. New output blocks can in
 *  principle be added here.
 */
/* GM: warning: I here ONLY leave standard data (pos, vel, mass, id, rho, entropy...) */
void fill_write_buffer(enum iofields blocknr, int *startindex, int pc, int type)
{
  int n, k, pindex;
  MyOutputFloat *fp;
  MyIDType *ip;


#ifdef PERIODIC
  MyFloat boxSize;
#endif

  double a3inv = 1, fac1, fac2;


  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      fac1 = 1 / (All.Time * All.Time);
      fac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
    }
  else
    {
      a3inv = fac1 = fac2 = 1;
    }


  fp = (MyOutputFloat *) CommBuffer;
  ip = (MyIDType *) CommBuffer;

  pindex = *startindex;

  switch (blocknr)
    {
    case IO_POS:		/* positions */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Pos[k];
#ifdef PERIODIC
		boxSize = All.BoxSize;
#ifdef LONG_X
		if(k == 0)
		  boxSize = All.BoxSize * LONG_X;
#endif
#ifdef LONG_Y
		if(k == 1)
		  boxSize = All.BoxSize * LONG_Y;
#endif
#ifdef LONG_Z
		if(k == 2)
		  boxSize = All.BoxSize * LONG_Z;
#endif
		while(fp[k] < 0)
		  fp[k] += boxSize;
		while(fp[k] >= boxSize)
		  fp[k] -= boxSize;
#ifdef DS_SHIFT_BOX
		fp[k] -= boxSize / 2;
#endif
#endif
	      }
	    n++;
	    fp += 3;
	  }
      break;

    case IO_VEL:		/* velocities */
      for(n = 0; n < pc; pindex++) /* GM: WARNING, velocities are currently saved as they are (no kick) */
	if(P[pindex].Type == type)
	  {
//	    dt_step = (P[pindex].TimeBin ? (((integertime) 1) << P[pindex].TimeBin) : 0);
//
//	    if(All.ComovingIntegrationOn)
//	      {
//		dt_gravkick =
//		  get_gravkick_factor(P[pindex].Ti_begstep, All.Ti_Current) -
//		  get_gravkick_factor(P[pindex].Ti_begstep, P[pindex].Ti_begstep + dt_step / 2);
//		dt_hydrokick =
//		  get_hydrokick_factor(P[pindex].Ti_begstep, All.Ti_Current) -
//		  get_hydrokick_factor(P[pindex].Ti_begstep, P[pindex].Ti_begstep + dt_step / 2);
//	      }
//	    else
//	      dt_gravkick = dt_hydrokick =
//		(All.Ti_Current - (P[pindex].Ti_begstep + dt_step / 2)) * All.Timebase_interval;
//
	    for(k = 0; k < 3; k++)
	      {
		fp[k] = P[pindex].Vel[k];

//		fp[k] += P[pindex].g.GravAccel[k] * dt_gravkick;
//		if(P[pindex].Type == 0)
//		  fp[k] += SphP[pindex].a.HydroAccel[k] * dt_hydrokick;
	      }

	    for(k = 0; k < 3; k++)
	      fp[k] *= sqrt(a3inv);

	    n++;
	    fp += 3;
	  }
      break;

    case IO_ID:		/* particle ID */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *ip++ = P[pindex].ID;
	    n++;
	  }
      break;

    case IO_MASS:		/* particle mass */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = P[pindex].Mass;
	    n++;
	  }
      break;

    case IO_U:			/* internal energy */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {

	    *fp++ =
	      DMAX(All.MinEgySpec,
		   SphP[pindex].Entropy / GAMMA_MINUS1 * pow(SphP[pindex].d.Density * a3inv, GAMMA_MINUS1));
	    n++;
	  }
      break;

    case IO_RHO:		/* density */
      for(n = 0; n < pc; pindex++)
	if(P[pindex].Type == type)
	  {
	    *fp++ = SphP[pindex].d.Density;
	    n++;
	  }
      break;


    case IO_LASTENTRY:
      endrun(213);
      break;
    }

  *startindex = pindex;
}




/*! This function tells the size of one data entry in each of the blocks
 *  defined for the output file.
 */
int get_bytes_per_blockelement(enum iofields blocknr, int mode)
{
  int bytes_per_blockelement = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
      if(mode)
	bytes_per_blockelement = 3 * sizeof(MyInputFloat);
      else
	bytes_per_blockelement = 3 * sizeof(MyOutputFloat);
      break;

    case IO_ID:
      bytes_per_blockelement = sizeof(MyIDType);
      break;


    case IO_MASS:
    case IO_U:
    case IO_RHO:
      if(mode)
	bytes_per_blockelement = sizeof(MyInputFloat);
      else
	bytes_per_blockelement = sizeof(MyOutputFloat);
      break;


    case IO_LASTENTRY:
      endrun(214);
      break;
    }

  return bytes_per_blockelement;
}

int get_datatype_in_block(enum iofields blocknr)
{
  int typekey;

  switch (blocknr)
    {
    case IO_ID:
#ifdef LONGIDS
      typekey = 2;		/* native long long */
#else
      typekey = 0;		/* native int */
#endif
      break;


    default:
      typekey = 1;		/* native MyOutputFloat */
      break;
    }

  return typekey;
}



int get_values_per_blockelement(enum iofields blocknr)
{
  int values = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
      values = 3;
      break;

    case IO_ID:
    case IO_MASS:
    case IO_U:
    case IO_RHO:
      values = 1;
      break;

    case IO_LASTENTRY:
      endrun(215);
      break;
    }
  return values;
}




/*! This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array.
 */
int get_particles_in_block(enum iofields blocknr, int *typelist)
{
  int i, nall, ntot_withmasses, ngas;

  nall = 0;
  ntot_withmasses = 0;

  for(i = 0; i < 6; i++)
    {
      typelist[i] = 0;

      if(header.npart[i] > 0)
	{
	  nall += header.npart[i];
	  typelist[i] = 1;
	}

      if(All.MassTable[i] == 0)
	ntot_withmasses += header.npart[i];
    }

  ngas = header.npart[0];

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ID:
      return nall;
      break;

    case IO_MASS:
      for(i = 0; i < 6; i++)
	{
	  typelist[i] = 0;
	  if(All.MassTable[i] == 0 && header.npart[i] > 0)
	    typelist[i] = 1;
	}
      return ntot_withmasses;
      break;

    case IO_U:
    case IO_RHO:
      for(i = 1; i < 6; i++)
	typelist[i] = 0;

      return ngas;
      break;

    case IO_LASTENTRY:
      endrun(216);
      break;
    }

  endrun(212);
  return 0;
}



/*! This function tells whether a block in the output file is present or not.
 */
int blockpresent(enum iofields blocknr)
{
  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ID:
    case IO_MASS:
    case IO_U:
    case IO_RHO:
      return 1;			/* always present */

    case IO_LASTENTRY:
      return 0;			/* will not occur */
    }


  return 0;			/* default: not present */
}




/*! This function associates a short 4-character block name with each block number.
 *  This is stored in front of each block for snapshot FileFormat=2.
 */
void get_Tab_IO_Label(enum iofields blocknr, char *label)
{
  switch (blocknr)
    {
    case IO_POS:
      strncpy(label, "POS ", 4);
      break;
    case IO_VEL:
      strncpy(label, "VEL ", 4);
      break;
    case IO_ID:
      strncpy(label, "ID  ", 4);
      break;
    case IO_MASS:
      strncpy(label, "MASS", 4);
      break;
    case IO_U:
      strncpy(label, "U   ", 4);
      break;
    case IO_RHO:
      strncpy(label, "RHO ", 4);
      break;
 
    case IO_LASTENTRY:
      endrun(217);
      break;
    }
}


void get_dataset_name(enum iofields blocknr, char *buf)
{
  strcpy(buf, "default");

  switch (blocknr)
    {
    case IO_POS:
      strcpy(buf, "Coordinates");
      break;
    case IO_VEL:
      strcpy(buf, "Velocities");
      break;
    case IO_ID:
      strcpy(buf, "ParticleIDs");
      break;
    case IO_MASS:
      strcpy(buf, "Masses");
      break;
    case IO_U:
      strcpy(buf, "InternalEnergy");
      break;
    case IO_RHO:
      strcpy(buf, "Density");
      break;

    case IO_LASTENTRY:
      endrun(218);
      break;
    }
}



/*! This function writes a snapshot file containing the data from processors
 *  'writeTask' to 'lastTask'. 'writeTask' is the one that actually writes.
 *  Each snapshot file contains a header first, then particle positions,
 *  velocities and ID's.  Then particle masses are written for those particle
 *  types with zero entry in MassTable.  After that, first the internal
 *  energies u, and then the density is written for the SPH particles.  If
 *  cooling is enabled, mean molecular weight and neutral hydrogen abundance
 *  are written for the gas particles. This is followed by the SPH smoothing
 *  length and further blocks of information, depending on included physics
 *  and compile-time flags.
 */
void write_file(char *fname, int writeTask, int lastTask)
{
  int type, bytes_per_blockelement, npart, nextblock, typelist[6];
  int n_for_this_task, n, p, pc, offset = 0, task, i;
  size_t blockmaxlen;
  int ntot_type[6], nn[6];
  enum iofields blocknr;
  char label[8];
  int bnr;
  int blksize;
  MPI_Status status;
  FILE *fd = 0;





#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}


  /* determine particle numbers of each type in file */

  if(ThisTask == writeTask)
    {
      for(n = 0; n < 6; n++)
	ntot_type[n] = n_type[n];

      for(task = writeTask + 1; task <= lastTask; task++)
	{
	  MPI_Recv(&nn[0], 6, MPI_INT, task, TAG_LOCALN, MYMPI_COMM_WORLD, &status);
	  for(n = 0; n < 6; n++)
	    ntot_type[n] += nn[n];
	}

      for(task = writeTask + 1; task <= lastTask; task++)
	MPI_Send(&ntot_type[0], 6, MPI_INT, task, TAG_N, MYMPI_COMM_WORLD);
    }
  else
    {
      MPI_Send(&n_type[0], 6, MPI_INT, writeTask, TAG_LOCALN, MYMPI_COMM_WORLD);
      MPI_Recv(&ntot_type[0], 6, MPI_INT, writeTask, TAG_N, MYMPI_COMM_WORLD, &status);
    }


  /* fill file header */

  for(n = 0; n < 6; n++)
    {
      header.npart[n] = ntot_type[n];
      header.npartTotal[n] = (unsigned int) ntot_type_all[n];
      header.npartTotalHighWord[n] = (unsigned int) (ntot_type_all[n] >> 32);
    }


  if(header.flag_ic_info == FLAG_SECOND_ORDER_ICS)
    header.flag_ic_info = FLAG_EVOLVED_2LPT;

  if(header.flag_ic_info == FLAG_ZELDOVICH_ICS)
    header.flag_ic_info = FLAG_EVOLVED_ZELDOVICH;

  if(header.flag_ic_info == FLAG_NORMALICS_2LPT)
    header.flag_ic_info = FLAG_EVOLVED_2LPT;

  if(header.flag_ic_info == 0 && All.ComovingIntegrationOn != 0)
    header.flag_ic_info = FLAG_EVOLVED_ZELDOVICH;

  for(n = 0; n < 6; n++)
    header.mass[n] = All.MassTable[n];



  header.time = All.Time;

  if(All.ComovingIntegrationOn)
    header.redshift = 1.0 / All.Time - 1;
  else
    header.redshift = 0;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;
  header.flag_entropy_instead_u = 0;


  header.num_files = All.NumFilesPerSnapshot;
  header.BoxSize = All.BoxSize;
  header.Omega0 = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;

#ifdef OUTPUT_IN_DOUBLEPRECISION
  header.flag_doubleprecision = 1;
#else
  header.flag_doubleprecision = 0;
#endif


  /* open file and write header */


  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == 3)
	{
	}
      else
	{
	  if(!(fd = fopen(fname, "w")))
	    {
	      printf("can't open file `%s' for writing snapshot.\n", fname);
	      endrun(123);
	    }

	  if(All.SnapFormat == 2)
	    {
	      blksize = sizeof(int) + 4 * sizeof(char);
	      SKIP;
	      my_fwrite((void *) "HEAD", sizeof(char), 4, fd);
	      nextblock = sizeof(header) + 2 * sizeof(int);
	      my_fwrite(&nextblock, sizeof(int), 1, fd);
	      SKIP;
	    }

	  blksize = sizeof(header);
	  SKIP;
	  my_fwrite(&header, sizeof(header), 1, fd);
	  SKIP;
	}

    }

  if((All.SnapFormat == 1 || All.SnapFormat == 2) && ThisTask == writeTask)
    {
      n_info = 0;
      InfoBlock = (struct info_block *) mymalloc("InfoBlock", sizeof(struct info_block) * 1000);

      for(bnr = 0; bnr < 1000; bnr++)
	{
	  blocknr = (enum iofields) bnr;


	  if(blocknr == IO_LASTENTRY)
	    break;

	  if(blockpresent(blocknr))
	    {
	      bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 0);
	      npart = get_particles_in_block(blocknr, &typelist[0]);

	      if(npart > 0)
		{
		  for(type = 0; type < 6; type++)
		    InfoBlock[n_info].is_present[type] = typelist[type];
		  InfoBlock[n_info].ndim = get_values_per_blockelement(blocknr);
		  get_Tab_IO_Label(blocknr, label);
		  for(i = 0; i < 4; i++)
		    InfoBlock[n_info].label[i] = label[i];
		  switch (get_datatype_in_block(blocknr))
		    {
		    case 0:
		      if(InfoBlock[n_info].ndim <= 1)
			strncpy(InfoBlock[n_info].type, "LONG    ", 8);
		      else
			strncpy(InfoBlock[n_info].type, "LONGN   ", 8);
		      break;
		    case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
		      if(InfoBlock[n_info].ndim <= 1)
			strncpy(InfoBlock[n_info].type, "DOUBLE  ", 8);
		      else
			strncpy(InfoBlock[n_info].type, "DOUBLEN ", 8);
		      break;
#else
		      if(InfoBlock[n_info].ndim <= 1)
			strncpy(InfoBlock[n_info].type, "FLOAT   ", 8);
		      else
			strncpy(InfoBlock[n_info].type, "FLOATN  ", 8);
		      break;
#endif
		    case 2:
		      if(InfoBlock[n_info].ndim <= 1)
			strncpy(InfoBlock[n_info].type, "LLONG   ", 8);
		      else
			strncpy(InfoBlock[n_info].type, "LLONGN  ", 8);
		      break;
		    }
		  n_info++;
		}
	    }
	}
    }

  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum iofields) bnr;


      if(blocknr == IO_LASTENTRY)
	break;

      if(blockpresent(blocknr))
	{
	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 0);

	  blockmaxlen = (size_t) ((All.BufferSize * 1024 * 1024) / bytes_per_blockelement);

	  npart = get_particles_in_block(blocknr, &typelist[0]);

	  if(npart > 0)
	    {
	      if(ThisTask == 0)
		{
		  char buf[1000];

		  get_dataset_name(blocknr, buf);
		  printf("writing block %d (%s)...\n", bnr, buf);
		}

	      if(ThisTask == writeTask)
		{

		  if(All.SnapFormat == 1 || All.SnapFormat == 2)
		    {

		    }
		}

	      for(type = 0; type < 6; type++)
		{
		  if(typelist[type])
		    {

		      for(task = writeTask, offset = 0; task <= lastTask; task++)
			{
			  if(task == ThisTask)
			    {
			      n_for_this_task = n_type[type];

			      for(p = writeTask; p <= lastTask; p++)
				if(p != ThisTask)
				  MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK,
					   MYMPI_COMM_WORLD);
			    }
			  else
			    MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MYMPI_COMM_WORLD,
				     &status);

			  while(n_for_this_task > 0)
			    {
			      pc = n_for_this_task;

			      if(pc > (int) blockmaxlen)
				pc = blockmaxlen;

			      if(ThisTask == task)
				fill_write_buffer(blocknr, &offset, pc, type);

			      if(ThisTask == writeTask && task != writeTask)
				MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task,
					 TAG_PDATA, MYMPI_COMM_WORLD, &status);

			      if(ThisTask != writeTask && task == ThisTask)
				MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask,
					  TAG_PDATA, MYMPI_COMM_WORLD);

			      if(ThisTask == writeTask)
				{
				  if(All.SnapFormat == 3)
				    {
				    }
				  else
				    {
				      my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
				    }
				}

			      n_for_this_task -= pc;
			    }
			}

		    }
		}

	      if(ThisTask == writeTask)
		{
		  if(All.SnapFormat == 1 || All.SnapFormat == 2)
		    {
		      SKIP;
		    }
		}
	    }
	}
    }


  if((All.SnapFormat == 1 || All.SnapFormat == 2) && ThisTask == writeTask)
    {
      if(All.SnapFormat == 2)
	{
	  blksize = sizeof(int) + 4 * sizeof(char);
	  SKIP;
	  my_fwrite((void *) "INFO", sizeof(char), 4, fd);
	  nextblock = n_info * sizeof(struct info_block) + 2 * sizeof(int);
	  my_fwrite(&nextblock, sizeof(int), 1, fd);
	  SKIP;
	}
      blksize = n_info * sizeof(struct info_block);
      SKIP;
      my_fwrite(InfoBlock, sizeof(struct info_block), n_info, fd);
      SKIP;
      myfree(InfoBlock);
    }


  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == 3)
	{
	}
      else
	{
	  fclose(fd);
	}
    }
}





/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if(size * nmemb > 0)
    {
      if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
	{
	  printf("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
	  fflush(stdout);
	  endrun(777);
	}
    }
  else
    nwritten = 0;

  return nwritten;
}


/*! This catches I/O errors occuring for fread(). In this case we
 *  better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if(size * nmemb == 0)
    return 0;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      if(feof(stream))
	printf("I/O error (fread) on task=%d has occured: end of file\n", ThisTask);
      else
	printf("I/O error (fread) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      fflush(stdout);
      endrun(778);
    }
  return nread;
}

/** This is so we don't have to preface every printf call with the if
    statement to only make it print once. */
void mpi_printf(const char *fmt, ...)
{
  if(ThisTask == 0)
    {
      va_list l;
      va_start(l, fmt);
      vprintf(fmt, l);
      fflush(stdout);
      va_end(l);
    }
}
