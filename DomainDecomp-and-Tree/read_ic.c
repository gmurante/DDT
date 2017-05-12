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
#include <string.h>


#include "allvars.h"
#include "proto.h"

/* This function reads initial conditions that are in the default file format
 * of Gadget, i.e. snapshot files can be used as input files.  However, when a
 * snapshot file is used as input, not all the information in the header is
 * used: THE STARTING TIME NEEDS TO BE SET IN THE PARAMETERFILE.
 * Alternatively, the code can be started with restartflag==2, then snapshots
 * from the code can be used as initial conditions-files without having to
 * change the parameterfile.  For gas particles, only the internal energy is
 * read, the density and mean molecular weight will be recomputed by the code.
 * When InitGasTemp>0 is given, the gas temperature will be initialzed to this
 * value assuming a mean colecular weight either corresponding to complete
 * neutrality, or full ionization.
 */

#ifdef AUTO_SWAP_ENDIAN_READIC
int swap_file = 8;
#endif




void read_ic(char *fname)
{
  int i, num_files, rest_files, ngroups, gr, filenr, masterTask, lastTask, groupMaster;
  double u_init, molecular_weight;
  char buf[500];

  CPU_Step[CPU_MISC] += measure_time();


  NumPart = 0;
  N_gas = 0;
  All.TotNumPart = 0;

  num_files = find_files(fname);


  rest_files = num_files;

  if(ThisTask == 0)
    printf("num_files = %i\n", num_files);


  while(rest_files > NTask)
    {
      sprintf(buf, "%s.%d", fname, ThisTask + (rest_files - NTask));
      if(All.ICFormat == 3)
	sprintf(buf, "%s.%d.hdf5", fname, ThisTask + (rest_files - NTask));

      ngroups = NTask / All.NumFilesWrittenInParallel;
      if((NTask % All.NumFilesWrittenInParallel))
	ngroups++;
      groupMaster = (ThisTask / ngroups) * ngroups;

      for(gr = 0; gr < ngroups; gr++)
	{
	  if(ThisTask == (groupMaster + gr))	/* ok, it's this processor's turn */
	  MPI_Barrier(MYMPI_COMM_WORLD);
	}

      rest_files -= NTask;
    }


  if(rest_files > 0)
    {
      distribute_file(rest_files, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

      if(num_files > 1)
	{
	  sprintf(buf, "%s.%d", fname, filenr);
	  if(All.ICFormat == 3)
	    sprintf(buf, "%s.%d.hdf5", fname, filenr);
	}
      else
	{
	  sprintf(buf, "%s", fname);
	  if(All.ICFormat == 3)
	    sprintf(buf, "%s.hdf5", fname);
	}

      ngroups = rest_files / All.NumFilesWrittenInParallel;
      if((rest_files % All.NumFilesWrittenInParallel))
	ngroups++;


      for(gr = 0; gr < ngroups; gr++)
	{

	  if((filenr / num_files) == gr)	/* ok, it's this processor's turn */
	    read_file(buf, masterTask, lastTask);
	  MPI_Barrier(MYMPI_COMM_WORLD);
	}
    }



  myfree(CommBuffer);


  if(header.flag_ic_info != FLAG_SECOND_ORDER_ICS)
    {
      /* this makes sure that masses are initialized in the case that the mass-block
         is empty for this particle type */
      for(i = 0; i < NumPart; i++)
	{
	  if(All.MassTable[P[i].Type] != 0)
	    P[i].Mass = All.MassTable[P[i].Type];
	}
    }




  u_init = (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.InitGasTemp;
  u_init *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;	/* unit conversion */

  if(All.InitGasTemp > 1.0e4)	/* assuming FULL ionization */
    molecular_weight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
  else				/* assuming NEUTRAL GAS */
    molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);

  u_init /= molecular_weight;

  All.InitGasU = u_init;


  if(RestartFlag == 0)
    {
      if(All.InitGasTemp > 0)
	{
	  for(i = 0; i < N_gas; i++)
	    {
	      if(ThisTask == 0 && i == 0 && SphP[i].Entropy == 0)
		printf("Initializing u from InitGasTemp !\n");

	      if(SphP[i].Entropy == 0)
		SphP[i].Entropy = All.InitGasU;

	      /* Note: the coversion to entropy will be done in the function init(),
	         after the densities have been computed */
	    }
	}


    }

  for(i = 0; i < N_gas; i++)
    SphP[i].Entropy = DMAX(All.MinEgySpec, SphP[i].Entropy);


  MPI_Barrier(MYMPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("reading done.\n");
      fflush(stdout);
    }

  if(ThisTask == 0)
    {
      printf("Total number of particles :  %d%09d\n\n",
	     (int) (All.TotNumPart / 1000000000), (int) (All.TotNumPart % 1000000000));
      fflush(stdout);
    }

  CPU_Step[CPU_SNAPSHOT] += measure_time();
}


/*! This function reads out the buffer that was filled with particle data.
 */
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type)
{
  int n, k;
  MyInputFloat *fp;
  MyIDType *ip;
#ifdef AUTO_SWAP_ENDIAN_READIC
  int vt, vpb;
  char *cp;
#endif

  fp = (MyInputFloat *) CommBuffer;
  ip = (MyIDType *) CommBuffer;

#ifdef AUTO_SWAP_ENDIAN_READIC
  if(blocknr != IO_LASTENTRY)
    {
      cp = (char *) CommBuffer;
      vt = get_datatype_in_block(blocknr);
      vpb = get_values_per_blockelement(blocknr);
      if(vt == 2)
	swap_Nbyte(cp, pc * vpb, 8);
      else
	{
#ifdef INPUT_IN_DOUBLEPRECISION
	  if(vt == 1)
	    swap_Nbyte(cp, pc * vpb, 8);
	  else
#endif
	    swap_Nbyte(cp, pc * vpb, 4);
	}
    }
#endif


  switch (blocknr)
    {
    case IO_POS:		/* positions */
      
      for(n = 0; n < pc; n++)
	{
	  for(k = 0; k < 3; k++)
	    P[offset + n].Pos[k] = *fp++;
	}

      for(n = 0; n < pc; n++)
	{
	  P[offset + n].Type = type;	/* initialize type here as well */

	}
      break;

    case IO_VEL:		/* velocities */
      for(n = 0; n < pc; n++)
	for(k = 0; k < 3; k++)
	  P[offset + n].Vel[k] = *fp++;
      break;

    case IO_ID:		/* particle ID */
      for(n = 0; n < pc; n++)
	{
	  P[offset + n].ID = *ip++;
	}
      break;
      
    case IO_MASS:		/* particle mass */
      for(n = 0; n < pc; n++)
	P[offset + n].Mass = *fp++;
      break;

    case IO_U:			/* temperature */
      for(n = 0; n < pc; n++)
	{
	  SphP[offset + n].Entropy = *fp++;
	}
      break;


    case IO_RHO:		/* density */
      for(n = 0; n < pc; n++)
	SphP[offset + n].d.Density = *fp++;
      break;


    case IO_LASTENTRY:
      endrun(220);
      break;
    }
}



/*! This function reads a snapshot file and distributes the data it contains
 *  to tasks 'readTask' to 'lastTask'.
 */
void read_file(char *fname, int readTask, int lastTask)
{
  size_t blockmaxlen;
  int i, n_in_file, n_for_this_task, ntask, pc, offset = 0, task;
  int blksize1, blksize2;
  MPI_Status status;
  FILE *fd = 0;
  int nall, nread;
  int type, bnr;
  char label[4], buf[500];
  int nstart, bytes_per_blockelement, npart, nextblock, typelist[6];
  enum iofields blocknr;
  size_t bytes;


#define SKIP  {my_fread(&blksize1,sizeof(int),1,fd);}
#define SKIP2  {my_fread(&blksize2,sizeof(int),1,fd);}

  if(ThisTask == readTask)
    {
      if(All.ICFormat == 1 || All.ICFormat == 2)
	{
	  if(!(fd = fopen(fname, "r")))
	    {
	      printf("can't open file `%s' for reading initial conditions.\n", fname);
	      endrun(123);
	    }


	  if(All.ICFormat == 2)
	    {
	      SKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_file = blksize1;
#endif
	      my_fread(&label, sizeof(char), 4, fd);
	      my_fread(&nextblock, sizeof(int), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_Nbyte((char *) &nextblock, 1, 4);
#endif
	      printf("Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2], label[3],
		     nextblock);
	      SKIP2;
	    }

	  SKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
	  if(All.ICFormat == 1)
	    {
	      if(blksize1 != 256)
		swap_file = 1;
	    }
#endif
	  my_fread(&header, sizeof(header), 1, fd);
	  SKIP2;
#ifdef AUTO_SWAP_ENDIAN_READIC
	  swap_Nbyte((char *) &blksize1, 1, 4);
	  swap_Nbyte((char *) &blksize2, 1, 4);
#endif

	  if(blksize1 != 256 || blksize2 != 256)
	    {
	      printf("incorrect header format\n");
	      fflush(stdout);
	      endrun(890);
	      /* Probable error is wrong size of fill[] in header file. Needs to be 256 bytes in total. */
	    }
#ifdef AUTO_SWAP_ENDIAN_READIC
	  swap_header();
#endif
	}



      for(task = readTask + 1; task <= lastTask; task++)
	{
	  MPI_Ssend(&header, sizeof(header), MPI_BYTE, task, TAG_HEADER, MYMPI_COMM_WORLD);
#ifdef AUTO_SWAP_ENDIAN_READIC
	  MPI_Ssend(&swap_file, sizeof(int), MPI_INT, task, TAG_SWAP, MYMPI_COMM_WORLD);
#endif
	}

    }
  else
    {
      MPI_Recv(&header, sizeof(header), MPI_BYTE, readTask, TAG_HEADER, MYMPI_COMM_WORLD, &status);
#ifdef AUTO_SWAP_ENDIAN_READIC
      MPI_Recv(&swap_file, sizeof(int), MPI_INT, readTask, TAG_SWAP, MYMPI_COMM_WORLD, &status);
#endif
    }

#ifdef INPUT_IN_DOUBLEPRECISION
  if(header.flag_doubleprecision == 0)
    {
      if(ThisTask == 0)
	printf
	  ("\nProblem: Code compiled with INPUT_IN_DOUBLEPRECISION, but input files are in single precision!\n");
      endrun(11);
    }
#else
  if(header.flag_doubleprecision)
    {
      if(ThisTask == 0)
	printf
	  ("\nProblem: Code not compiled with INPUT_IN_DOUBLEPRECISION, but input files are in double precision!\n");
      endrun(10);
    }
#endif


  if(All.TotNumPart == 0)
    {
      if(header.num_files <= 1)
	for(i = 0; i < 6; i++)
	  {
	    header.npartTotal[i] = header.npart[i];
	  }

      All.TotN_gas = header.npartTotal[0] + (((long long) header.npartTotalHighWord[0]) << 32);

      for(i = 0, All.TotNumPart = 0; i < 6; i++)
	{
	  All.TotNumPart += header.npartTotal[i];
	  All.TotNumPart += (((long long) header.npartTotalHighWord[i]) << 32);
	}



      for(i = 0; i < 6; i++)
	All.MassTable[i] = header.mass[i];

      All.MaxPart = (int) (All.PartAllocFactor * (All.TotNumPart / NTask));	/* sets the maximum number of particles that may */
      All.MaxPartSph = (int) (All.PartAllocFactor * (All.TotN_gas / NTask));	/* sets the maximum number of particles that may
										   reside on a processor */

      allocate_memory();

      if(!(CommBuffer = mymalloc("CommBuffer", bytes = All.BufferSize * 1024 * 1024)))
	{
	  printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(2);
	}

      if(RestartFlag >= 2)
	{
	  All.Time = All.TimeBegin = header.time;
	  //	  set_cosmo_factors_for_current_time();
	}

    }

  if(ThisTask == readTask)
    {
      for(i = 0, n_in_file = 0; i < 6; i++)
	n_in_file += header.npart[i];

      printf("\nreading file `%s' on task=%d (contains %d particles.)\n"
	     "distributing this file to tasks %d-%d\n"
	     "Type 0 (gas):   %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 1 (halo):  %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 2 (disk):  %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 3 (bulge): %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 4 (stars): %8d  (tot=%6d%09d) masstab=%g\n"
	     "Type 5 (bndry): %8d  (tot=%6d%09d) masstab=%g\n\n", fname, ThisTask, n_in_file, readTask,
	     lastTask, header.npart[0], (int) (header.npartTotal[0] / 1000000000),
	     (int) (header.npartTotal[0] % 1000000000), All.MassTable[0], header.npart[1],
	     (int) (header.npartTotal[1] / 1000000000), (int) (header.npartTotal[1] % 1000000000),
	     All.MassTable[1], header.npart[2], (int) (header.npartTotal[2] / 1000000000),
	     (int) (header.npartTotal[2] % 1000000000), All.MassTable[2], header.npart[3],
	     (int) (header.npartTotal[3] / 1000000000), (int) (header.npartTotal[3] % 1000000000),
	     All.MassTable[3], header.npart[4], (int) (header.npartTotal[4] / 1000000000),
	     (int) (header.npartTotal[4] % 1000000000), All.MassTable[4], header.npart[5],
	     (int) (header.npartTotal[5] / 1000000000), (int) (header.npartTotal[5] % 1000000000),
	     All.MassTable[5]);
      fflush(stdout);
    }


  ntask = lastTask - readTask + 1;


  /* to collect the gas particles all at the beginning (in case several
     snapshot files are read on the current CPU) we move the collisionless
     particles such that a gap of the right size is created */

  for(type = 0, nall = 0; type < 6; type++)
    {
      n_in_file = header.npart[type];

      n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
	n_for_this_task++;


      if(type == 0)
	{
	  if(N_gas + n_for_this_task > All.MaxPartSph)
	    {
	      printf("Not enough space on task=%d for SPH particles (space for %d, need at least %d)\n",
		     ThisTask, All.MaxPartSph, N_gas + n_for_this_task);
	      fflush(stdout);
	      endrun(172);
	    }
	}

      nall += n_for_this_task;
    }

  if(NumPart + nall > All.MaxPart)
    {
      printf("Not enough space on task=%d (space for %d, need at least %d)\n",
	     ThisTask, All.MaxPart, NumPart + nall);
      fflush(stdout);
      endrun(173);
    }

  memmove(&P[N_gas + nall], &P[N_gas], (NumPart - N_gas) * sizeof(struct particle_data));
  nstart = N_gas;



  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum iofields) bnr;

      if(blocknr == IO_LASTENTRY)
	break;


      if(RestartFlag == 5 && blocknr > IO_MASS)	/* if we only do power spectra, we don't need to read other blocks beyond the mass */
	continue;



      if(blockpresent(blocknr))
	{
	  if(RestartFlag == 0 && blocknr > IO_U )
		  continue;	/* ignore all other blocks in initial conditions */


	  if(ThisTask == readTask)
	    {
	      get_dataset_name(blocknr, buf);
	      printf("reading block %d (%s)...\n", bnr, buf);
	      fflush(stdout);
	    }

	  bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 1);

	  blockmaxlen = (size_t) ((All.BufferSize * 1024 * 1024) / bytes_per_blockelement);

	  npart = get_particles_in_block(blocknr, &typelist[0]);

	  if(npart > 0)
	    {
	      if(blocknr != IO_LASTENTRY)
		if(ThisTask == readTask)
		  {

		    if(All.ICFormat == 2)
		      {
			get_Tab_IO_Label(blocknr, label);
			find_block(label, fd);
		      }

		    if(All.ICFormat == 1 || All.ICFormat == 2)
		      SKIP;
		  }

	      for(type = 0, offset = 0, nread = 0; type < 6; type++)
		{
		  n_in_file = header.npart[type];
		  if(typelist[type] == 0)
		    {
		      n_for_this_task = n_in_file / ntask;
		      if((ThisTask - readTask) < (n_in_file % ntask))
			n_for_this_task++;

		      offset += n_for_this_task;
		    }
		  else
		    {
		      for(task = readTask; task <= lastTask; task++)
			{
			  n_for_this_task = n_in_file / ntask;
			  if((task - readTask) < (n_in_file % ntask))
			    n_for_this_task++;

			  if(task == ThisTask)
			    if(NumPart + n_for_this_task > All.MaxPart)
			      {
				printf("too many particles. %d %d %d\n", NumPart, n_for_this_task,
				       All.MaxPart);
				endrun(1313);
			      }

			  do
			    {
			      pc = n_for_this_task;

			      if(pc > (int) blockmaxlen)
				pc = blockmaxlen;

			      if(ThisTask == readTask)
				{
				  if(All.ICFormat == 1 || All.ICFormat == 2)
				    {
				      if(blocknr != IO_LASTENTRY)
					{
					  my_fread(CommBuffer, bytes_per_blockelement, pc, fd);
					  nread += pc;
					}
				      else
					{
					  nread += pc;
					}
				    }
				}

			      if(ThisTask == readTask && task != readTask && pc > 0)
				MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task,
					  TAG_PDATA, MYMPI_COMM_WORLD);

			      if(ThisTask != readTask && task == ThisTask && pc > 0)
				MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, readTask,
					 TAG_PDATA, MYMPI_COMM_WORLD, &status);

			      if(ThisTask == task)
				{
				  empty_read_buffer(blocknr, nstart + offset, pc, type);

				  offset += pc;
				}


			      n_for_this_task -= pc;
			    }
			  while(n_for_this_task > 0);
			}
		    }
		}

	      if(ThisTask == readTask)
		{
		  if(blocknr != IO_LASTENTRY)
		    if(All.ICFormat == 1 || All.ICFormat == 2)
		      {
			SKIP2;


#ifdef AUTO_SWAP_ENDIAN_READIC
			swap_Nbyte((char *) &blksize1, 1, 4);
			swap_Nbyte((char *) &blksize2, 1, 4);
#endif
			if(blksize1 != blksize2)
			  {
			    printf("incorrect block-sizes detected!\n");
			    printf("Task=%d   blocknr=%d  blksize1=%d  blksize2=%d\n", ThisTask, bnr,
				   blksize1, blksize2);
			    if(blocknr == IO_ID)
			      {
				printf
				  ("Possible mismatch of 32bit and 64bit ID's in IC file and GADGET compilation !\n");
			      }
			    fflush(stdout);
			    endrun(1889);
			  }
		      }
		}
	    }
	}
    }



  for(type = 0; type < 6; type++)
    {
      n_in_file = header.npart[type];

      n_for_this_task = n_in_file / ntask;
      if((ThisTask - readTask) < (n_in_file % ntask))
	n_for_this_task++;

      NumPart += n_for_this_task;

      if(type == 0)
	N_gas += n_for_this_task;

    }

  if(ThisTask == readTask)
    {
      if(All.ICFormat == 1 || All.ICFormat == 2)
	fclose(fd);
    }


}



/*! This function determines on how many files a given snapshot is distributed.
 */
int find_files(char *fname)
{
  FILE *fd;
  char buf[200], buf1[200];
  int dummy;

  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);

  if(All.ICFormat == 3)
    {
      sprintf(buf, "%s.%d.hdf5", fname, 0);
      sprintf(buf1, "%s.hdf5", fname);
    }

  if(All.ICFormat == 3)
    {
      if(ThisTask == 0)
	printf("Code wasn't compiled with HDF5 support enabled!\n");
      endrun(0);
    }

  header.num_files = 0;

  if(ThisTask == 0)
    {
      if((fd = fopen(buf, "r")))
	{
	  if(All.ICFormat == 1 || All.ICFormat == 2)
	    {
	      if(All.ICFormat == 2)
		{
		  my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
		  swap_file = dummy;
#endif
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		}

	      my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      if(All.ICFormat == 1)
		{
		  if(dummy == 256)
		    swap_file = 8;
		  else
		    swap_file = 1;
		}
#endif
	      my_fread(&header, sizeof(header), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_header();
#endif
	      my_fread(&dummy, sizeof(dummy), 1, fd);
	    }
	  fclose(fd);

	}
    }

#ifdef AUTO_SWAP_ENDIAN_READIC
  MPI_Bcast(&swap_file, sizeof(int), MPI_INT, 0, MYMPI_COMM_WORLD);
#endif
  MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MYMPI_COMM_WORLD);

  if(header.num_files > 0)
    return header.num_files;

  if(ThisTask == 0)
    {
      if((fd = fopen(buf1, "r")))
	{
	  if(All.ICFormat == 1 || All.ICFormat == 2)
	    {
	      if(All.ICFormat == 2)
		{
		  my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
		  swap_file = dummy;
#endif
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		  my_fread(&dummy, sizeof(dummy), 1, fd);
		}

	      my_fread(&dummy, sizeof(dummy), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      if(All.ICFormat == 1)
		{
		  if(dummy == 256)
		    swap_file = 8;
		  else
		    swap_file = 1;
		}
#endif
	      my_fread(&header, sizeof(header), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	      swap_header();
#endif
	      my_fread(&dummy, sizeof(dummy), 1, fd);
	    }
	  fclose(fd);

	  header.num_files = 1;
	}
    }

#ifdef AUTO_SWAP_ENDIAN_READIC
  MPI_Bcast(&swap_file, sizeof(int), MPI_INT, 0, MYMPI_COMM_WORLD);
#endif
  MPI_Bcast(&header, sizeof(header), MPI_BYTE, 0, MYMPI_COMM_WORLD);

  if(header.num_files > 0)
    return header.num_files;

  if(ThisTask == 0)
    {
      printf("\nCan't find initial conditions file.");
      printf("neither as '%s'\nnor as '%s'\n", buf, buf1);
      fflush(stdout);
    }

  endrun(0);
  return 0;
}



/*! This function assigns a certain number of files to processors, such that
 *  each processor is exactly assigned to one file, and the number of cpus per
 *  file is as homogenous as possible. The number of files may at most be
 *  equal to the number of processors.
 */
void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master,
		     int *last)
{
  int ntask, filesleft, filesright, tasksleft;

  if(nfiles > 1)
    {
      ntask = lasttask - firsttask + 1;

      filesleft = (int) ((((double) (ntask / 2)) / ntask) * nfiles);
      if(filesleft <= 0)
	filesleft = 1;
      if(filesleft >= nfiles)
	filesleft = nfiles - 1;

      filesright = nfiles - filesleft;

      tasksleft = ntask / 2;
      distribute_file(filesleft, firstfile, firsttask, firsttask + tasksleft - 1, filenr, master, last);
      distribute_file(filesright, firstfile + filesleft, firsttask + tasksleft, lasttask, filenr, master,
		      last);
    }
  else
    {
      if(ThisTask >= firsttask && ThisTask <= lasttask)
	{
	  *filenr = firstfile;
	  *master = firsttask;
	  *last = lasttask;
	}
    }
}




#ifdef AUTO_SWAP_ENDIAN_READIC
/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to swap ENDIAN -------------------------------*/
/*-------- char *data:    Pointer to the data ---------------------------------*/
/*-------- int n:         Number of elements to swap --------------------------*/
/*-------- int m:         Size of single element to swap ----------------------*/
/*--------                int,float = 4 ---------------------------------------*/
/*--------                double    = 8 ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
void swap_Nbyte(char *data, int n, int m)
{
  int i, j;
  char old_data[16];

  if(swap_file != 8)
    {
      for(j = 0; j < n; j++)
	{
	  memcpy(&old_data[0], &data[j * m], m);
	  for(i = 0; i < m; i++)
	    {
	      data[j * m + i] = old_data[m - i - 1];
	    }
	}
    }
}

/*------------------------------------------------------------------*/
/*----------- procedure to swap header if needed -------------------*/
/*------------------------------------------------------------------*/

void swap_header()
{
  swap_Nbyte((char *) &header.npart, 6, 4);
  swap_Nbyte((char *) &header.mass, 6, 8);
  swap_Nbyte((char *) &header.time, 1, 8);
  swap_Nbyte((char *) &header.redshift, 1, 8);
  swap_Nbyte((char *) &header.flag_sfr, 1, 4);
  swap_Nbyte((char *) &header.flag_feedback, 1, 4);
  swap_Nbyte((char *) &header.npartTotal, 6, 4);
  swap_Nbyte((char *) &header.flag_cooling, 1, 4);
  swap_Nbyte((char *) &header.num_files, 1, 4);
  swap_Nbyte((char *) &header.BoxSize, 1, 8);
  swap_Nbyte((char *) &header.Omega0, 1, 8);
  swap_Nbyte((char *) &header.OmegaLambda, 1, 8);
  swap_Nbyte((char *) &header.HubbleParam, 1, 8);
  swap_Nbyte((char *) &header.flag_stellarage, 1, 4);
  swap_Nbyte((char *) &header.flag_metals, 1, 4);
  swap_Nbyte((char *) &header.npartTotalHighWord, 6, 4);
  swap_Nbyte((char *) &header.flag_entropy_instead_u, 1, 4);
  swap_Nbyte((char *) &header.flag_doubleprecision, 1, 4);
  swap_Nbyte((char *) &header.flag_ic_info, 1, 4);
  swap_Nbyte((char *) &header.lpt_scalingfactor, 1, 4);
}

#endif

/*---------------------- Routine find a block in a snapfile -------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- char *label:   4 byte identifyer for block -------------------------*/
/*-------- returns length of block found, -------------------------------------*/
/*-------- the file fd points to starting point of block ----------------------*/
/*-----------------------------------------------------------------------------*/
void find_block(char *label, FILE * fd)
{
  unsigned int blocksize = 0, blksize;
  char blocklabel[5] = { "    " };

#define FBSKIP  {my_fread(&blksize,sizeof(int),1,fd);}

  rewind(fd);


  while(!feof(fd) && blocksize == 0)
    {
      FBSKIP;
#ifdef AUTO_SWAP_ENDIAN_READIC
      swap_file = blksize;
      swap_Nbyte((char *) &blksize, 1, 4);
#endif
      if(blksize != 8)
	{
	  printf("Incorrect Format (blksize=%u)!\n", blksize);
	  exit(1891);
	}
      else
	{
	  my_fread(blocklabel, 4 * sizeof(char), 1, fd);
	  my_fread(&blocksize, sizeof(int), 1, fd);
#ifdef AUTO_SWAP_ENDIAN_READIC
	  swap_Nbyte((char *) &blocksize, 1, 4);
#endif
	  /*
	     printf("Searching <%c%c%c%c>, found Block <%s> with %d bytes\n",
	     label[0],label[1],label[2],label[3],blocklabel,blocksize);
	   */
	  FBSKIP;
	  if(strncmp(label, blocklabel, 4) != 0)
	    {
	      fseek(fd, blocksize, 1);
	      blocksize = 0;
	    }
	}
    }
  if(feof(fd))
    {
      printf("Block '%c%c%c%c' not found !\n", label[0], label[1], label[2], label[3]);
      fflush(stdout);
      endrun(1890);
    }
}
