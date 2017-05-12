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


#include "allvars.h"
#include "proto.h"


/*! \file domain.c
 *  \brief code for domain decomposition
 *
 *  This file contains the code for the domain decomposition of the
 *  simulation volume.  The domains are constructed from disjoint subsets
 *  of the leaves of a fiducial top-level tree that covers the full
 *  simulation volume. Domain boundaries hence run along tree-node
 *  divisions of a fiducial global BH tree. As a result of this method, the
 *  tree force are in principle strictly independent of the way the domains
 *  are cut. The domain decomposition can be carried out for an arbitrary
 *  number of CPUs. Individual domains are not cubical, but spatially
 *  coherent since the leaves are traversed in a Peano-Hilbert order and
 *  individual domains form segments along this order.  This also ensures
 *  that each domain has a small surface to volume ratio, which minimizes
 *  communication.
 */


#define REDUC_FAC      0.98


/*! toGo[task*NTask + partner] gives the number of particles in task 'task'
 *  that have to go to task 'partner'
 */
static int *toGo, *toGoSph;
static int *toGet, *toGetSph;
static int *list_NumPart;
static int *list_N_gas;
static int *list_load;
static int *list_loadsph;
static double *list_work;
static double *list_worksph;

/* GM: these are instantiated in restart.c and 
  used when a restart is done, now I instantiate them here
 (but shouldn't be used...)
 extern int old_MaxPart, new_MaxPart; */
int old_MaxPart, new_MaxPart; 


static struct local_topnode_data
{
  peanokey Size;		/*!< number of Peano-Hilbert mesh-cells represented by top-level node */
  peanokey StartKey;		/*!< first Peano-Hilbert key in top-level node */
  long long Count;		/*!< counts the number of particles in this top-level node */
  double Cost;
  double SphCost;
  int Daughter;			/*!< index of first daughter cell (out of 8) of top-level node */
  int Leaf;			/*!< if the node is a leaf, this gives its number when all leaves are traversed in Peano-Hilbert order */
  int Parent;
  int PIndex;			/*!< first particle in node */
}
 *topNodes;			/*!< points to the root node of the top-level tree */

static struct peano_hilbert_data
{
  peanokey key;
  int index;
}
 *mp;



static void domain_insertnode(struct local_topnode_data *treeA, struct local_topnode_data *treeB, int noA,
			      int noB);
static void domain_add_cost(struct local_topnode_data *treeA, int noA, long long count, double cost,
			    double sphcost);



static float *domainWork;	/*!< a table that gives the total "work" due to the particles stored by each processor */
static float *domainWorkSph;	/*!< a table that gives the total "work" due to the particles stored by each processor */
static int *domainCount;	/*!< a table that gives the total number of particles held by each processor */
static int *domainCountSph;	/*!< a table that gives the total number of SPH particles held by each processor */


static int domain_allocated_flag = 0;

static int maxLoad, maxLoadsph;


static double totgravcost, gravcost, totsphcost, sphcost;
static long long totpartcount;

static int UseAllParticles;



/*! This is the main routine for the domain decomposition.  It acts as a
 *  driver routine that allocates various temporary buffers, maps the
 *  particles back onto the periodic box if needed, and then does the
 *  domain decomposition, and a final Peano-Hilbert order of all particles
 *  as a tuning measure.
 */
void domain_Decomposition(int UseAllTimeBins, int SaveKeys)
{
  int i, ret, retsum, diff, highest_bin_to_include;
  size_t bytes, all_bytes;
  double t0, t1;

  UseAllParticles = UseAllTimeBins;

  CPU_Step[CPU_MISC] += measure_time();

  /* GM: this is for getting the right position with multiple timesteps.
     not needed now, but I leavhe it here as a memo 
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_current != All.Ti_Current)
    drift_particle(i, All.Ti_Current); */

  force_treefree();
  domain_free();

  if(old_MaxPart)
    {
      All.MaxPart = new_MaxPart;
      old_MaxPart = 0;
    }


#ifdef PERIODIC
  do_box_wrapping();		/* map the particles back onto the box */
#endif


  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type > 5 || P[i].Type < 0)
	{
	  printf("task=%d:  P[i=%d].Type=%d\n", ThisTask, i, P[i].Type);
	  endrun(112411);
	}
    }


  TreeReconstructFlag = 1;	/* ensures that new tree will be constructed */

  /* we take the closest cost factor */
  if(UseAllParticles)
    highest_bin_to_include = All.HighestOccupiedTimeBin;
  else
    highest_bin_to_include = All.HighestActiveTimeBin;

  for(i = 1, TakeLevel = 0, diff = abs(All.LevelToTimeBin[0] - highest_bin_to_include); i < GRAVCOSTLEVELS;
      i++)
    {
      if(diff > abs(All.LevelToTimeBin[i] - highest_bin_to_include))
	{
	  TakeLevel = i;
	  diff = abs(All.LevelToTimeBin[i] - highest_bin_to_include);
	}
    }

  if(ThisTask == 0)
    {
      printf("domain decomposition... LevelToTimeBin[TakeLevel=%d]=%d  (presently allocated=%g MB)\n",
	     TakeLevel, All.LevelToTimeBin[TakeLevel], AllocatedBytes / (1024.0 * 1024.0));
      fflush(stdout);
    }

  t0 = second();

  do
    {
      domain_allocate();

      all_bytes = 0;

      Key = (peanokey *) mymalloc("domain_key", bytes = (sizeof(peanokey) * All.MaxPart));
      all_bytes += bytes;

      toGo = (int *) mymalloc("toGo", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      toGoSph = (int *) mymalloc("toGoSph", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      toGet = (int *) mymalloc("toGet", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      toGetSph = (int *) mymalloc("toGetSph", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      list_NumPart = (int *) mymalloc("list_NumPart", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      list_N_gas = (int *) mymalloc("list_N_gas", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      list_load = (int *) mymalloc("list_load", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      list_loadsph = (int *) mymalloc("list_loadsph", bytes = (sizeof(int) * NTask));
      all_bytes += bytes;
      list_work = (double *) mymalloc("list_work", bytes = (sizeof(double) * NTask));
      all_bytes += bytes;
      list_worksph = (double *) mymalloc("list_worksph", bytes = (sizeof(double) * NTask));
      all_bytes += bytes;
      domainWork = (float *) mymalloc("domainWork", bytes = (MaxTopNodes * sizeof(float)));
      all_bytes += bytes;
      domainWorkSph = (float *) mymalloc("domainWorkSph", bytes = (MaxTopNodes * sizeof(float)));
      all_bytes += bytes;
      domainCount = (int *) mymalloc("domainCount", bytes = (MaxTopNodes * sizeof(int)));
      all_bytes += bytes;
      domainCountSph = (int *) mymalloc("domainCountSph", bytes = (MaxTopNodes * sizeof(int)));
      all_bytes += bytes;


      topNodes = (struct local_topnode_data *) mymalloc("topNodes", bytes =
							(MaxTopNodes * sizeof(struct local_topnode_data)));
      all_bytes += bytes;

      if(ThisTask == 0)
	{
	  printf
	    ("use of %g MB of temporary storage for domain decomposition... (presently allocated=%g MB)\n",
	     all_bytes / (1024.0 * 1024.0), AllocatedBytes / (1024.0 * 1024.0));
	  fflush(stdout);
	}

      maxLoad = (int) (All.MaxPart * REDUC_FAC);
      maxLoadsph = (int) (All.MaxPartSph * REDUC_FAC);

      report_memory_usage(&HighMark_domain, "DOMAIN");

      ret = domain_decompose();

      /* copy what we need for the topnodes */
      for(i = 0; i < NTopnodes; i++)
	{
	  TopNodes[i].StartKey = topNodes[i].StartKey;
	  TopNodes[i].Size = topNodes[i].Size;
	  TopNodes[i].Daughter = topNodes[i].Daughter;
	  TopNodes[i].Leaf = topNodes[i].Leaf;
	}

      myfree(topNodes);

      myfree(domainCountSph);
      myfree(domainCount);
      myfree(domainWorkSph);
      myfree(domainWork);
      myfree(list_worksph);
      myfree(list_work);
      myfree(list_loadsph);
      myfree(list_load);
      myfree(list_N_gas);
      myfree(list_NumPart);
      myfree(toGetSph);
      myfree(toGet);
      myfree(toGoSph);
      myfree(toGo);


      MPI_Allreduce(&ret, &retsum, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
      if(retsum)
	{
	  myfree(Key);

	  domain_free();

	  if(ThisTask == 0)
	    printf("Increasing TopNodeAllocFactor=%g  ", All.TopNodeAllocFactor);

	  All.TopNodeAllocFactor *= 1.3;

	  if(ThisTask == 0)
	    {
	      printf("new value=%g\n", All.TopNodeAllocFactor);
	      fflush(stdout);
	    }

	  if(All.TopNodeAllocFactor > 1000)
	    {
	      if(ThisTask == 0)
		printf("something seems to be going seriously wrong here. Stopping.\n");
	      fflush(stdout);
	      endrun(781);
	    }
	}
    }
  while(retsum);

  t1 = second();

  if(ThisTask == 0)
    {
      printf("domain decomposition done. (took %g sec)\n", timediff(t0, t1));
      fflush(stdout);
    }

  CPU_Step[CPU_DOMAIN] += measure_time();

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type > 5 || P[i].Type < 0)
	{
	  printf("task=%d:  P[i=%d].Type=%d\n", ThisTask, i, P[i].Type);
	  endrun(111111);
	}


    }

#ifdef PEANOHILBERT
  peano_hilbert_order();
  
  CPU_Step[CPU_PEANO] += measure_time();
#endif

  myfree(Key);

  memmove(TopNodes + NTopnodes, DomainTask, NTopnodes * sizeof(int));

  TopNodes = (struct topnode_data *) myrealloc(TopNodes, bytes =
					       (NTopnodes * sizeof(struct topnode_data) +
						NTopnodes * sizeof(int)));
  if(ThisTask == 0)
    printf("Freed %g MByte in top-level domain structure\n",
	   (MaxTopNodes - NTopnodes) * sizeof(struct topnode_data) / (1024.0 * 1024.0));

  DomainTask = (int *) (TopNodes + NTopnodes);

  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);

  CPU_Step[CPU_DOMAIN] += measure_time();

  /* GM: this is needed with multiple timestep, unused now, I leave it here as a memo 
     reconstruct_timebins(); */
}

/*! This function allocates all the stuff that will be required for the tree-construction/walk later on */
void domain_allocate(void)
{
  size_t bytes, all_bytes = 0;

  MaxTopNodes = (int) (All.TopNodeAllocFactor * All.MaxPart + 1);

  DomainStartList = (int *) mymalloc("DomainStartList", bytes = (NTask * MULTIPLEDOMAINS * sizeof(int)));
  all_bytes += bytes;

  DomainEndList = (int *) mymalloc("DomainEndList", bytes = (NTask * MULTIPLEDOMAINS * sizeof(int)));
  all_bytes += bytes;

  TopNodes = (struct topnode_data *) mymalloc("TopNodes", bytes =
					      (MaxTopNodes * sizeof(struct topnode_data) +
					       MaxTopNodes * sizeof(int)));
  all_bytes += bytes;

  DomainTask = (int *) (TopNodes + MaxTopNodes);

  if(ThisTask == 0)
    printf("Allocated %g MByte for top-level domain structure\n", all_bytes / (1024.0 * 1024.0));

  domain_allocated_flag = 1;
}

void domain_free(void)
{
  if(domain_allocated_flag)
    {
      myfree(TopNodes);
      myfree(DomainEndList);
      myfree(DomainStartList);
      domain_allocated_flag = 0;
    }
}

static struct topnode_data *save_TopNodes;
static int *save_DomainStartList, *save_DomainEndList;

void domain_free_trick(void)
{
  if(domain_allocated_flag)
    {
      save_TopNodes = TopNodes;
      save_DomainEndList = DomainEndList;
      save_DomainStartList = DomainStartList;
      domain_allocated_flag = 0;
    }
  else
    endrun(131231);
}

void domain_allocate_trick(void)
{
  domain_allocated_flag = 1;

  TopNodes = save_TopNodes;
  DomainEndList = save_DomainEndList;
  DomainStartList = save_DomainStartList;
}




double domain_particle_costfactor(int i)
{
  return 0.1 + P[i].GravCost[TakeLevel];
}


/*! This function carries out the actual domain decomposition for all
 *  particle types. It will try to balance the work-load for each domain,
 *  as estimated based on the P[i]-GravCost values.  The decomposition will
 *  respect the maximum allowed memory-imbalance given by the value of
 *  PartAllocFactor.
 */
int domain_decompose(void)
{
  int i, status;
  long long sumtogo, sumload, sumloadsph;
  int maxload, maxloadsph, multipledomains = MULTIPLEDOMAINS;
  double sumwork, maxwork, sumworksph, maxworksph;


  for(i = 0; i < 6; i++)
    NtypeLocal[i] = 0;

  gravcost = sphcost = 0;


  for(i = 0; i < NumPart; i++)
    {
      NtypeLocal[P[i].Type]++;
      
      gravcost += domain_particle_costfactor(i);
      if(P[i].Type == 0)
	if(TimeBinActive[P[i].TimeBin] || UseAllParticles)
	  sphcost += 1.0;
    }

  /* because Ntype[] is of type `long long', we cannot do a simple
   * MPI_Allreduce() to sum the total particle numbers 
   */
  sumup_large_ints(6, NtypeLocal, Ntype);

  for(i = 0, totpartcount = 0; i < 6; i++)
    totpartcount += Ntype[i];

  MPI_Allreduce(&gravcost, &totgravcost, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(&sphcost, &totsphcost, 1, MPI_DOUBLE, MPI_SUM, MYMPI_COMM_WORLD);

#ifndef UNEQUALSOFTENINGS
  for(i = 0; i < 6; i++)
    if(Ntype[i] > 0)
      break;

  int j;

  for(j = i + 1; j < 6; j++)
    {
      if(Ntype[j] > 0)
	if(All.SofteningTable[j] != All.SofteningTable[i])
	  {
	    if(ThisTask == 0)
	      {
		fprintf(stdout, "Code was not compiled with UNEQUALSOFTENINGS, but some of the\n");
		fprintf(stdout, "softening lengths are unequal nevertheless.\n");
		fprintf(stdout, "This is not allowed.\n");
	      }
	    endrun(0);
	  }
    }
#endif


  /* determine global dimensions of domain grid */
  domain_findExtent();

  if(domain_determineTopTree())
    return 1;

  /* find the split of the domain grid */
  domain_findSplit_work_balanced(multipledomains * NTask, NTopleaves);
  domain_assign_load_or_work_balanced(1, multipledomains);
  
  status = domain_check_memory_bound(multipledomains);
  
  if(status != 0)		/* the optimum balanced solution violates memory constraint, let's try something different */
    {
      if(ThisTask == 0)
	printf
	  ("Note: the domain decomposition is suboptimum because the ceiling for memory-imbalance is reached\n");
      
#ifndef OVERRIDE_STOP_FOR_SUBOPTIMUM_DOMAINS
      if(ThisTask == 0)
	{
	  printf("We better stop. (OVERRIDE_STOP_FOR_SUBOPTIMUM_DOMAINS is not set)\n");
	}
      endrun(0);
#endif
      
      domain_findSplit_load_balanced(multipledomains * NTask, NTopleaves);
      domain_assign_load_or_work_balanced(0, multipledomains);
      
      status = domain_check_memory_bound(multipledomains);
      
      if(status != 0)
	{
	  if(ThisTask == 0)
	    printf("No domain decomposition that stays within memory bounds is possible.\n");
	  endrun(0);
	}
    }
  

  if(ThisTask == 0)
    {
      sumload = maxload = sumloadsph = maxloadsph = 0;
      sumwork = sumworksph = maxwork = maxworksph = 0;

      for(i = 0; i < NTask; i++)
	{
	  sumload += list_load[i];
	  sumloadsph += list_loadsph[i];
	  sumwork += list_work[i];
	  sumworksph += list_worksph[i];

	  if(list_load[i] > maxload)
	    maxload = list_load[i];

	  if(list_loadsph[i] > maxloadsph)
	    maxloadsph = list_loadsph[i];



	  if(list_work[i] > maxwork)
	    maxwork = list_work[i];

	  if(list_worksph[i] > maxworksph)
	    maxworksph = list_worksph[i];
	}

      printf("gravity work-load balance=%g   memory-balance=%g   SPH work-load balance=%g\n",
	     maxwork / (sumwork / NTask), maxload / (((double) sumload) / NTask),
	     maxworksph / ((sumworksph + 1.0e-30) / NTask));
    }

  for(i = 0; i < NumPart; i++)
    {

      int no = 0;

      while(topNodes[no].Daughter >= 0)
	no = topNodes[no].Daughter + (Key[i] - topNodes[no].StartKey) / (topNodes[no].Size / 8);

      no = topNodes[no].Leaf;

      int task = DomainTask[no];

      if(task != ThisTask)
	P[i].Type |= 32;
    }


  int iter = 0, ret;
  size_t exchange_limit;

  do
    {
      exchange_limit = FreeBytes - NTask * (24 * sizeof(int) + 16 * sizeof(MPI_Request));

#ifdef MPISENDRECV_CHECKSUM
      exchange_limit -= All.BufferSize * 1024 * 1024;
#endif
      if(exchange_limit <= 0)
	{
	  printf("task=%d: exchange_limit=%d\n", ThisTask, (int) exchange_limit);
	  endrun(1223);
	}

      /* determine for each cpu how many particles have to be shifted to other cpus */
      ret = domain_countToGo(exchange_limit);

      for(i = 0, sumtogo = 0; i < NTask; i++)
	sumtogo += toGo[i];

      sumup_longs(1, &sumtogo, &sumtogo);

      if(ThisTask == 0)
	{
	  printf("iter=%d exchange of %d%09d particles (ret=%d)\n", iter,
		 (int) (sumtogo / 1000000000), (int) (sumtogo % 1000000000), ret);
	  fflush(stdout);
	}

      domain_exchange();

      iter++;
    }
  while(ret > 0);

  return 0;
}

int domain_check_memory_bound(int multipledomains)
{
  int ta, m, i;
  int load, sphload, max_load, max_sphload;
  double work, worksph;


  max_load = max_sphload = 0;

  for(ta = 0; ta < NTask; ta++)
    {
      load = sphload = 0;
      work = worksph = 0;

      for(m = 0; m < multipledomains; m++)
	for(i = DomainStartList[ta * multipledomains + m]; i <= DomainEndList[ta * multipledomains + m]; i++)
	  {
	    load += domainCount[i];
	    sphload += domainCountSph[i];
	    work += domainWork[i];
	    worksph += domainWorkSph[i];
	  }

      list_load[ta] = load;
      list_loadsph[ta] = sphload;
      list_work[ta] = work;
      list_worksph[ta] = worksph;

      if(load > max_load)
	max_load = load;
      if(sphload > max_sphload)
	max_sphload = sphload;
    }


  if(max_load > maxLoad)
    {
      if(ThisTask == 0)
	{
	  printf("desired memory imbalance=%g  (limit=%d, needed=%d)\n",
		 (max_load * All.PartAllocFactor) / maxLoad, maxLoad, max_load);
	}

      domain_print();

      return 1;
    }

  if(max_sphload > maxLoadsph)
    {
      if(ThisTask == 0)
	{
	  printf("desired memory imbalance=%g  (SPH) (limit=%d, needed=%d)\n",
		 (max_sphload * All.PartAllocFactor) / maxLoadsph, maxLoadsph, max_sphload);
	}

      domain_print();

      return 1;
    }

  return 0;
}

void domain_print(void)
{
  int ta;

  if(ThisTask == 0)
    {
      printf("load loadsph work worksph");
      printf("\n");
    }
  MPI_Barrier(MYMPI_COMM_WORLD);

  for(ta = 0; ta < NTask; ta++)
    {
      if(ThisTask == ta)
	{
	  printf("%d ", ThisTask);
	  printf("%d %d %g %g ", list_load[ta], list_loadsph[ta], list_work[ta], list_worksph[ta]);
	  printf("\n");
	  fflush(stdout);
	}
      MPI_Barrier(MYMPI_COMM_WORLD);
    }
}

void domain_exchange(void)
{
  int count_togo = 0, count_togo_sph = 0, count_get = 0, count_get_sph = 0;
  int *count, *count_sph, *offset, *offset_sph;
  int *count_recv, *count_recv_sph, *offset_recv, *offset_recv_sph;
  int i, n, ngrp, no, target;
  struct particle_data *partBuf;
  struct sph_particle_data *sphBuf;
  peanokey *keyBuf;

  count = (int *) mymalloc("count", NTask * sizeof(int));
  count_sph = (int *) mymalloc("count_sph", NTask * sizeof(int));
  offset = (int *) mymalloc("offset", NTask * sizeof(int));
  offset_sph = (int *) mymalloc("offset_sph", NTask * sizeof(int));

  count_recv = (int *) mymalloc("count_recv", NTask * sizeof(int));
  count_recv_sph = (int *) mymalloc("count_recv_sph", NTask * sizeof(int));
  offset_recv = (int *) mymalloc("offset_recv", NTask * sizeof(int));
  offset_recv_sph = (int *) mymalloc("offset_recv_sph", NTask * sizeof(int));




  int prec_offset, prec_count;
  int *decrease;

  decrease = (int *) mymalloc("decrease", NTask * sizeof(int));

  for(i = 1, offset_sph[0] = 0, decrease[0] = 0; i < NTask; i++)
    {
      offset_sph[i] = offset_sph[i - 1] + toGoSph[i - 1];
      decrease[i] = toGoSph[i - 1];
    }

  prec_offset = offset_sph[NTask - 1] + toGoSph[NTask - 1];



  offset[0] = prec_offset;
  for(i = 1; i < NTask; i++)
    offset[i] = offset[i - 1] + (toGo[i - 1] - decrease[i]);

  myfree(decrease);

  for(i = 0; i < NTask; i++)
    {
      count_togo += toGo[i];
      count_togo_sph += toGoSph[i];

      count_get += toGet[i];
      count_get_sph += toGetSph[i];

    }

  partBuf = (struct particle_data *) mymalloc("partBuf", count_togo * sizeof(struct particle_data));
  sphBuf = (struct sph_particle_data *) mymalloc("sphBuf", count_togo_sph * sizeof(struct sph_particle_data));
  keyBuf = (peanokey *) mymalloc("keyBuf", count_togo * sizeof(peanokey));

  for(i = 0; i < NTask; i++)
    count[i] = count_sph[i] = 0;


  for(n = 0; n < NumPart; n++)
    {
      if((P[n].Type & (32 + 16)) == (32 + 16))
	{
	  P[n].Type &= 15;

	  no = 0;

	  while(topNodes[no].Daughter >= 0)
	    no = topNodes[no].Daughter + (Key[n] - topNodes[no].StartKey) / (topNodes[no].Size / 8);

	  no = topNodes[no].Leaf;

	  target = DomainTask[no];

	  if(P[n].Type == 0)
	    {
	      partBuf[offset_sph[target] + count_sph[target]] = P[n];
	      keyBuf[offset_sph[target] + count_sph[target]] = Key[n];
	      sphBuf[offset_sph[target] + count_sph[target]] = SphP[n];
	      count_sph[target]++;
	    }

	  else
	    {
	      partBuf[offset[target] + count[target]] = P[n];
	      keyBuf[offset[target] + count[target]] = Key[n];
	      count[target]++;
	    }


	  if(P[n].Type == 0)
	    {
	      P[n] = P[N_gas - 1];
	      SphP[n] = SphP[N_gas - 1];
	      Key[n] = Key[N_gas - 1];

	      P[N_gas - 1] = P[NumPart - 1];
	      Key[N_gas - 1] = Key[NumPart - 1];

	      NumPart--;
	      N_gas--;
	      n--;
	    }


	  else
	    {
	      P[n] = P[NumPart - 1];
	      Key[n] = Key[NumPart - 1];
	      NumPart--;
	      n--;
	    }
	}
    }

  int count_totget;

  count_totget = count_get_sph;

  if(count_totget)
    {
      memmove(P + N_gas + count_totget, P + N_gas, (NumPart - N_gas) * sizeof(struct particle_data));
      memmove(Key + N_gas + count_totget, Key + N_gas, (NumPart - N_gas) * sizeof(peanokey));
    }




  for(i = 0; i < NTask; i++)
    {
      count_recv_sph[i] = toGetSph[i];
      count_recv[i] = toGet[i] - toGetSph[i];
    }


  for(i = 1, offset_recv_sph[0] = N_gas; i < NTask; i++)
    offset_recv_sph[i] = offset_recv_sph[i - 1] + count_recv_sph[i - 1];
  prec_count = N_gas + count_get_sph;


  offset_recv[0] = NumPart - N_gas + prec_count;

  for(i = 1; i < NTask; i++)
    offset_recv[i] = offset_recv[i - 1] + count_recv[i - 1];

#ifndef NO_ISEND_IRECV_IN_DOMAIN

  int n_requests = 0, max_requests = 10;
  MPI_Request *requests;


  requests = (MPI_Request *) mymalloc("requests", max_requests * NTask * sizeof(MPI_Request));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  if(count_recv_sph[target] > 0)
	    {
	      MPI_Irecv(P + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct particle_data),
			MPI_BYTE, target, TAG_PDATA_SPH, MYMPI_COMM_WORLD, &requests[n_requests++]);

	      MPI_Irecv(Key + offset_recv_sph[target], count_recv_sph[target] * sizeof(peanokey),
			MPI_BYTE, target, TAG_KEY_SPH, MYMPI_COMM_WORLD, &requests[n_requests++]);

	      MPI_Irecv(SphP + offset_recv_sph[target],
			count_recv_sph[target] * sizeof(struct sph_particle_data), MPI_BYTE, target,
			TAG_SPHDATA, MYMPI_COMM_WORLD, &requests[n_requests++]);
	    }




	  if(count_recv[target] > 0)
	    {
	      MPI_Irecv(P + offset_recv[target], count_recv[target] * sizeof(struct particle_data),
			MPI_BYTE, target, TAG_PDATA, MYMPI_COMM_WORLD, &requests[n_requests++]);

	      MPI_Irecv(Key + offset_recv[target], count_recv[target] * sizeof(peanokey),
			MPI_BYTE, target, TAG_KEY, MYMPI_COMM_WORLD, &requests[n_requests++]);
	    }
	}
    }


  MPI_Barrier(MYMPI_COMM_WORLD);	/* not really necessary, but this will guarantee that all receives are
					   posted before the sends, which helps the stability of MPI on 
					   bluegene, and perhaps some mpich1-clusters */

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  if(count_sph[target] > 0)
	    {
	      MPI_Isend(partBuf + offset_sph[target], count_sph[target] * sizeof(struct particle_data),
			MPI_BYTE, target, TAG_PDATA_SPH, MYMPI_COMM_WORLD, &requests[n_requests++]);

	      MPI_Isend(keyBuf + offset_sph[target], count_sph[target] * sizeof(peanokey),
			MPI_BYTE, target, TAG_KEY_SPH, MYMPI_COMM_WORLD, &requests[n_requests++]);

	      MPI_Isend(sphBuf + offset_sph[target], count_sph[target] * sizeof(struct sph_particle_data),
			MPI_BYTE, target, TAG_SPHDATA, MYMPI_COMM_WORLD, &requests[n_requests++]);
	    }




	  if(count[target] > 0)
	    {
	      MPI_Isend(partBuf + offset[target], count[target] * sizeof(struct particle_data),
			MPI_BYTE, target, TAG_PDATA, MYMPI_COMM_WORLD, &requests[n_requests++]);

	      MPI_Isend(keyBuf + offset[target], count[target] * sizeof(peanokey),
			MPI_BYTE, target, TAG_KEY, MYMPI_COMM_WORLD, &requests[n_requests++]);
	    }
	}
    }

  MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);

  if(n_requests > max_requests * NTask)
    {
      printf("Not enough memory reserved for requests: %d > %d !\n", n_requests, max_requests * NTask);
      endrun(52097);
    }

  myfree(requests);
#else

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  if(count_sph[target] > 0 || count_recv_sph[target] > 0)
	    {
	      MPI_Sendrecv(partBuf + offset_sph[target], count_sph[target] * sizeof(struct particle_data),
			   MPI_BYTE, target, TAG_PDATA_SPH,
			   P + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct particle_data),
			   MPI_BYTE, target, TAG_PDATA_SPH, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);

	      MPI_Sendrecv(sphBuf + offset_sph[target], count_sph[target] * sizeof(struct sph_particle_data),
			   MPI_BYTE, target, TAG_SPHDATA,
			   SphP + offset_recv_sph[target],
			   count_recv_sph[target] * sizeof(struct sph_particle_data), MPI_BYTE, target,
			   TAG_SPHDATA, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);

	      MPI_Sendrecv(keyBuf + offset_sph[target], count_sph[target] * sizeof(peanokey),
			   MPI_BYTE, target, TAG_KEY_SPH,
			   Key + offset_recv_sph[target], count_recv_sph[target] * sizeof(peanokey),
			   MPI_BYTE, target, TAG_KEY_SPH, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }



	  if(count[target] > 0 || count_recv[target] > 0)
	    {
	      MPI_Sendrecv(partBuf + offset[target], count[target] * sizeof(struct particle_data),
			   MPI_BYTE, target, TAG_PDATA,
			   P + offset_recv[target], count_recv[target] * sizeof(struct particle_data),
			   MPI_BYTE, target, TAG_PDATA, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);

	      MPI_Sendrecv(keyBuf + offset[target], count[target] * sizeof(peanokey),
			   MPI_BYTE, target, TAG_KEY,
			   Key + offset_recv[target], count_recv[target] * sizeof(peanokey),
			   MPI_BYTE, target, TAG_KEY, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }


#endif

  NumPart += count_get;
  N_gas += count_get_sph;



  if(NumPart > All.MaxPart)
    {
      printf("Task=%d NumPart=%d All.MaxPart=%d\n", ThisTask, NumPart, All.MaxPart);
      endrun(787878);
    }

  if(N_gas > All.MaxPartSph)
    endrun(787879);


  myfree(keyBuf);
  myfree(sphBuf);
  myfree(partBuf);


  myfree(offset_recv_sph);
  myfree(offset_recv);
  myfree(count_recv_sph);
  myfree(count_recv);

  myfree(offset_sph);
  myfree(offset);
  myfree(count_sph);
  myfree(count);
}



void domain_findSplit_work_balanced(int ncpu, int ndomain)
{
  int i, start, end;
  double work, worksph, workavg, work_before, workavg_before;
  double load, fac_work, fac_load, fac_worksph;

  for(i = 0, work = load = worksph = 0; i < ndomain; i++)
    {
      work += domainWork[i];
      load += domainCount[i];
      worksph += domainWorkSph[i];
    }

  if(worksph > 0)
    {
      /* in this case we give equal weight to gravitational work-load, SPH work load, and particle load.
       */
      fac_work = 0.333 / work;
      fac_load = 0.333 / load;
      fac_worksph = 0.333 / worksph;
    }
  else
    {
      /* in this case we give equal weight to gravitational work-load and particle load.
       * The final pieces should have at most imbalance 2.0 in either of the two
       */
      fac_work = 0.5 / work;
      fac_load = 0.5 / load;
      fac_worksph = 0.0;
    }

  workavg = 1.0 / ncpu;

  work_before = workavg_before = 0;

  start = 0;

  for(i = 0; i < ncpu; i++)
    {
      work = 0;
      end = start;

      work += fac_work * domainWork[end] + fac_load * domainCount[end] + fac_worksph * domainWorkSph[end];

      while((work + work_before < workavg + workavg_before) || (i == ncpu - 1 && end < ndomain - 1))
	{
	  if((ndomain - end) > (ncpu - i))
	    end++;
	  else
	    break;

	  work += fac_work * domainWork[end] + fac_load * domainCount[end] + fac_worksph * domainWorkSph[end];
	}

      DomainStartList[i] = start;
      DomainEndList[i] = end;

      work_before += work;
      workavg_before += workavg;
      start = end + 1;
    }
}

static struct domain_segments_data
{
  int task, start, end;
  double work;
  double load;
  double load_activesph;
  double normalized_load;
}
 *domainAssign;


struct queue_data
{
  int first, last;
  int *next;
  int *previous;
  double *value;
}
queues[3];

struct tasklist_data
{
  double work;
  double load;
  double load_activesph;
  int count;
}
 *tasklist;

int domain_sort_task(const void *a, const void *b)
{
  if(((struct domain_segments_data *) a)->task < (((struct domain_segments_data *) b)->task))
    return -1;

  if(((struct domain_segments_data *) a)->task > (((struct domain_segments_data *) b)->task))
    return +1;

  return 0;
}

int domain_sort_load(const void *a, const void *b)
{
  if(((struct domain_segments_data *) a)->normalized_load >
     (((struct domain_segments_data *) b)->normalized_load))
    return -1;

  if(((struct domain_segments_data *) a)->normalized_load <
     (((struct domain_segments_data *) b)->normalized_load))
    return +1;

  return 0;
}

void domain_assign_load_or_work_balanced(int mode, int multipledomains)
{
  double target_work_balance, target_load_balance, target_load_activesph_balance;
  double value, target_max_balance, best_balance;
  double tot_work, tot_load, tot_loadactivesph;

  int best_queue, target, next, prev;
  int i, n, q, ta;

  domainAssign = (struct domain_segments_data *) mymalloc("domainAssign",
							  multipledomains * NTask *
							  sizeof(struct domain_segments_data));

  tasklist = (struct tasklist_data *) mymalloc("tasklist", NTask * sizeof(struct tasklist_data));

  for(ta = 0; ta < NTask; ta++)
    {
      tasklist[ta].work = 0;
      tasklist[ta].load = 0;
      tasklist[ta].load_activesph = 0;
      tasklist[ta].count = 0;
    }

  tot_work = 0;
  tot_load = 0;
  tot_loadactivesph = 0;

  for(n = 0; n < multipledomains * NTask; n++)
    {
      domainAssign[n].start = DomainStartList[n];
      domainAssign[n].end = DomainEndList[n];
      domainAssign[n].work = 0;
      domainAssign[n].load = 0;
      domainAssign[n].load_activesph = 0;

      for(i = DomainStartList[n]; i <= DomainEndList[n]; i++)
	{
	  domainAssign[n].work += domainWork[i];
	  domainAssign[n].load += domainCount[i];
	  domainAssign[n].load_activesph += domainWorkSph[i];
	}

      if(mode == 0)
	domainAssign[n].work = 1;

      tot_work += domainAssign[n].work;
      tot_load += domainAssign[n].load;
      tot_loadactivesph += domainAssign[n].load_activesph;
    }

  for(n = 0; n < multipledomains * NTask; n++)
    {
      domainAssign[n].normalized_load =
	domainAssign[n].work / (tot_work + 1.0e-30) +
	domainAssign[n].load_activesph / (tot_loadactivesph + 1.0e-30);
    }

  qsort(domainAssign, multipledomains * NTask, sizeof(struct domain_segments_data), domain_sort_load);

/* initialize three queues */
  for(q = 0; q < 3; q++)
    {
      queues[q].next = (int *) mymalloc("queues[q].next", NTask * sizeof(int));
      queues[q].previous = (int *) mymalloc("queues[q].previous", NTask * sizeof(int));
      queues[q].value = (double *) mymalloc("queues[q].value", NTask * sizeof(double));

      for(ta = 0; ta < NTask; ta++)
	{
	  queues[q].next[ta] = ta + 1;
	  queues[q].previous[ta] = ta - 1;
	  queues[q].value[ta] = 0;
	}
      queues[q].previous[0] = -1;
      queues[q].next[NTask - 1] = -1;
      queues[q].first = 0;
      queues[q].last = NTask - 1;
    }

  for(n = 0; n < multipledomains * NTask; n++)
    {
      /* need to decide, which of the tasks that has the lowest load in one of the three quantities is best */
      for(q = 0, best_balance = 1.0e30, best_queue = 0; q < 3; q++)
	{
	  target = queues[q].first;

	  while(tasklist[target].count == multipledomains)
	    target = queues[q].next[target];

	  target_work_balance = (domainAssign[n].work + tasklist[target].work) / (tot_work + 1.0e-30);
	  target_load_balance = (domainAssign[n].load + tasklist[target].load) / (tot_load + 1.0e-30);
	  target_load_activesph_balance = (domainAssign[n].load_activesph + tasklist[target].load_activesph)
	    / (tot_loadactivesph + 1.0e-30);

	  target_max_balance = target_work_balance;
	  if(target_max_balance < target_load_balance)
	    target_max_balance = target_load_balance;
	  if(target_max_balance < target_load_activesph_balance)
	    target_max_balance = target_load_activesph_balance;

	  if(target_max_balance < best_balance)
	    {
	      best_balance = target_max_balance;
	      best_queue = q;
	    }
	}

      /* Now we now the best queue, and hence the best target task. Assign this piece to this task */
      target = queues[best_queue].first;

      while(tasklist[target].count == multipledomains)
	target = queues[best_queue].next[target];

      domainAssign[n].task = target;
      tasklist[target].work += domainAssign[n].work;
      tasklist[target].load += domainAssign[n].load;
      tasklist[target].load_activesph += domainAssign[n].load_activesph;
      tasklist[target].count++;

      /* now we need to remove the element 'target' from the 3 queue's and reinsert it */
      for(q = 0; q < 3; q++)
	{
	  switch (q)
	    {
	    case 0:
	      value = tasklist[target].work;
	      break;
	    case 1:
	      value = tasklist[target].load;
	      break;
	    case 2:
	      value = tasklist[target].load_activesph;
	      break;
	    default:
	      value = 0;
	      break;
	    }

	  /* now remove the element target */
	  prev = queues[q].previous[target];
	  next = queues[q].next[target];

	  if(prev >= 0)		/* previous exists */
	    queues[q].next[prev] = next;
	  else
	    queues[q].first = next;	/* we remove the head of the queue */


	  if(next >= 0)		/* next exists */
	    queues[q].previous[next] = prev;
	  else
	    queues[q].last = prev;	/* we remove the end of the queue */

	  /* now we insert the element again, in an ordered fashion, starting from the end of the queue */
	  if(queues[q].last >= 0)
	    {
	      ta = queues[q].last;

	      while(value < queues[q].value[ta])
		{
		  ta = queues[q].previous[ta];
		  if(ta < 0)
		    break;
		}

	      if(ta < 0)	/* we insert the element as the first element */
		{
		  queues[q].next[target] = queues[q].first;
		  queues[q].previous[queues[q].first] = target;
		  queues[q].first = target;
		}
	      else
		{
		  /* insert behind ta */
		  queues[q].next[target] = queues[q].next[ta];
		  if(queues[q].next[ta] >= 0)
		    queues[q].previous[queues[q].next[ta]] = target;
		  else
		    queues[q].last = target;	/* we insert a new last element */
		  queues[q].previous[target] = ta;
		  queues[q].next[ta] = target;
		}
	    }
	  else
	    {
	      /* queue was empty */
	      queues[q].previous[target] = queues[q].next[target] = -1;
	      queues[q].first = queues[q].last = target;
	    }

	  queues[q].value[target] = value;
	}
    }

  qsort(domainAssign, multipledomains * NTask, sizeof(struct domain_segments_data), domain_sort_task);


  for(n = 0; n < multipledomains * NTask; n++)
    {
      DomainStartList[n] = domainAssign[n].start;
      DomainEndList[n] = domainAssign[n].end;

      for(i = DomainStartList[n]; i <= DomainEndList[n]; i++)
	DomainTask[i] = domainAssign[n].task;
    }

/* free the queues */
  for(q = 2; q >= 0; q--)
    {
      myfree(queues[q].value);
      myfree(queues[q].previous);
      myfree(queues[q].next);
    }

  myfree(tasklist);

  myfree(domainAssign);
}


void domain_findSplit_load_balanced(int ncpu, int ndomain)
{
  int i, start, end;
  double load, loadavg, load_before, loadavg_before, fac_load, fac;

  for(i = 0, load = 0; i < ndomain; i++)
    {
      load += domainCount[i];
    }

  fac = 1;
  fac = 1. / fac;

  fac_load = fac / load;

  loadavg = 1.0 / ncpu;

  load_before = loadavg_before = 0;

  start = 0;

  for(i = 0; i < ncpu; i++)
    {
      load = 0;
      end = start;

      load += fac_load * domainCount[end];
      while((load + load_before < loadavg + loadavg_before) || (i == ncpu - 1 && end < ndomain - 1))
	{
	  if((ndomain - end) > (ncpu - i))
	    end++;
	  else
	    break;

	  load += fac_load * domainCount[end];
	}

      DomainStartList[i] = start;
      DomainEndList[i] = end;

      load_before += load;
      loadavg_before += loadavg;
      start = end + 1;
    }
}








/*! This function determines how many particles that are currently stored
 *  on the local CPU have to be moved off according to the domain
 *  decomposition.
 */
int domain_countToGo(size_t nlimit)
{
  int n, no, ret, retsum, mpi_overflow_flag;
  size_t mpi_Particle_limit = 2147483647 / sizeof(struct particle_data);
  size_t mpi_SphParticle_limit = 2147483647 / sizeof(struct sph_particle_data);
  size_t package;

  for(n = 0; n < NTask; n++)
    {
      toGo[n] = 0;
      toGoSph[n] = 0;
    }

  package = (sizeof(struct particle_data) + sizeof(struct sph_particle_data) + sizeof(peanokey));

  if(package >= nlimit)
    endrun(212);


  for(n = 0; n < NumPart && package < nlimit; n++)
    {
      mpi_overflow_flag = 0;

      if(P[n].Type & 32)
	{
	  no = 0;

	  while(topNodes[no].Daughter >= 0)
	    no = topNodes[no].Daughter + (Key[n] - topNodes[no].StartKey) / (topNodes[no].Size / 8);

	  no = topNodes[no].Leaf;

	  if(DomainTask[no] != ThisTask)
	    {
	      toGo[DomainTask[no]] += 1;
	      nlimit -= sizeof(struct particle_data) + sizeof(peanokey);

	      if(toGo[DomainTask[no]] >= mpi_Particle_limit)
		mpi_overflow_flag = 1;

	      if((P[n].Type & 15) == 0)
		{
		  toGoSph[DomainTask[no]] += 1;
		  nlimit -= sizeof(struct sph_particle_data);

		  if(toGoSph[DomainTask[no]] >= mpi_SphParticle_limit)
		    mpi_overflow_flag = 1;
		}

	      P[n].Type |= 16;	/* flag this particle for export */
	      if(mpi_overflow_flag == 1)
		{
		  printf("Task %d: Applying size limit for MPI \n", ThisTask);
		  break;
		}
	    }
	}
    }

  MPI_Alltoall(toGo, 1, MPI_INT, toGet, 1, MPI_INT, MYMPI_COMM_WORLD);
  MPI_Alltoall(toGoSph, 1, MPI_INT, toGetSph, 1, MPI_INT, MYMPI_COMM_WORLD);

  if(package < nlimit && mpi_overflow_flag == 0)
    ret = 0;
  else
    ret = 1;


  MPI_Allreduce(&ret, &retsum, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

  if(retsum)
    {
      /* in this case, we are not guaranteed that the temporary state after
         the partial exchange will actually observe the particle limits on all
         processors... we need to test this explicitly and rework the exchange
         such that this is guaranteed. This is actually a rather non-trivial
         constraint. */

      MPI_Allgather(&NumPart, 1, MPI_INT, list_NumPart, 1, MPI_INT, MYMPI_COMM_WORLD);
      MPI_Allgather(&N_gas, 1, MPI_INT, list_N_gas, 1, MPI_INT, MYMPI_COMM_WORLD);

      int flag, flagsum, ntoomany, ta, i, target;
      int count_togo, count_toget, count_togo_sph, count_toget_sph;


      do
	{
	  flagsum = 0;

	  do
	    {
	      flag = 0;

	      for(ta = 0; ta < NTask; ta++)
		{
		  if(ta == ThisTask)
		    {
		      count_togo = count_toget = 0;
		      count_togo_sph = count_toget_sph = 0;
		      for(i = 0; i < NTask; i++)
			{
			  count_togo += toGo[i];
			  count_toget += toGet[i];
			  count_togo_sph += toGoSph[i];
			  count_toget_sph += toGetSph[i];
			}
		    }
		  MPI_Bcast(&count_togo, 1, MPI_INT, ta, MYMPI_COMM_WORLD);
		  MPI_Bcast(&count_toget, 1, MPI_INT, ta, MYMPI_COMM_WORLD);
		  MPI_Bcast(&count_togo_sph, 1, MPI_INT, ta, MYMPI_COMM_WORLD);
		  MPI_Bcast(&count_toget_sph, 1, MPI_INT, ta, MYMPI_COMM_WORLD);


		  int ifntoomany;

		  ntoomany = list_N_gas[ta] + count_toget_sph - count_togo_sph - All.MaxPartSph;
		  ifntoomany = (ntoomany > 0);


		  if(ifntoomany)
		    {
		      if(ThisTask == 0)
			{
			  if(ntoomany > 0)
			    {
			      printf
				("exchange needs to be modified because I can't receive %d SPH-particles on task=%d\n",
				 ntoomany, ta);
			      if(flagsum > 25)
				printf
				  ("list_N_gas[ta=%d]=%d  count_toget_sph=%d count_togo_sph=%d MaxPartSPH=%d\n",
				   ta, list_N_gas[ta], count_toget_sph, count_togo_sph, All.MaxPartSph);
			      fflush(stdout);
			    }


			}


		      flag = 1;
		      i = flagsum % NTask;

		      while(ifntoomany)
			{
			  if(i == ThisTask)
			    {
			      if(toGoSph[ta] > 0)
				if(ntoomany > 0)
				  {
				    toGoSph[ta]--;
				    count_toget_sph--;
				    count_toget--;
				    ntoomany--;
				  }

			    }

			  MPI_Bcast(&ntoomany, 1, MPI_INT, i, MYMPI_COMM_WORLD);
			  MPI_Bcast(&count_toget, 1, MPI_INT, i, MYMPI_COMM_WORLD);
			  MPI_Bcast(&count_toget_sph, 1, MPI_INT, i, MYMPI_COMM_WORLD);
			  i++;
			  if(i >= NTask)
			    i = 0;

			  ifntoomany = (ntoomany > 0);

			}
		    }

		  if((ntoomany = list_NumPart[ta] + count_toget - count_togo - All.MaxPart) > 0)
		    {
		      if(ThisTask == 0)
			{
			  printf
			    ("exchange needs to be modified because I can't receive %d particles on task=%d\n",
			     ntoomany, ta);
			  if(flagsum > 25)
			    printf("list_NumPart[ta=%d]=%d  count_toget=%d count_togo=%d MaxPart=%d\n",
				   ta, list_NumPart[ta], count_toget, count_togo, All.MaxPart);
			  fflush(stdout);
			}

		      flag = 1;
		      i = flagsum % NTask;
		      while(ntoomany)
			{
			  if(i == ThisTask)
			    {
			      if(toGo[ta] > 0)
				{
				  toGo[ta]--;
				  count_toget--;
				  ntoomany--;
				}
			    }

			  MPI_Bcast(&ntoomany, 1, MPI_INT, i, MYMPI_COMM_WORLD);
			  MPI_Bcast(&count_toget, 1, MPI_INT, i, MYMPI_COMM_WORLD);

			  i++;
			  if(i >= NTask)
			    i = 0;
			}
		    }
		}
	      flagsum += flag;

	      if(ThisTask == 0)
		{
		  printf("flagsum = %d\n", flagsum);
		  fflush(stdout);
		  if(flagsum > 100)
		    endrun(1013);
		}
	    }
	  while(flag);

	  if(flagsum)
	    {
	      int *local_toGo, *local_toGoSph;

	      local_toGo = (int *) mymalloc("	      local_toGo", NTask * sizeof(int));
	      local_toGoSph = (int *) mymalloc("	      local_toGoSph", NTask * sizeof(int));



	      for(n = 0; n < NTask; n++)
		{
		  local_toGo[n] = 0;
		  local_toGoSph[n] = 0;
		}

	      for(n = 0; n < NumPart; n++)
		{
		  if(P[n].Type & 32)
		    {
		      P[n].Type &= (15 + 32);	/* clear 16 */

		      no = 0;

		      while(topNodes[no].Daughter >= 0)
			no =
			  topNodes[no].Daughter + (Key[n] - topNodes[no].StartKey) / (topNodes[no].Size / 8);

		      no = topNodes[no].Leaf;

		      target = DomainTask[no];

		      if((P[n].Type & 15) == 0)
			{
			  if(local_toGoSph[target] < toGoSph[target] && local_toGo[target] < toGo[target])
			    {
			      local_toGo[target] += 1;
			      local_toGoSph[target] += 1;
			      P[n].Type |= 16;
			    }
			}


		      else
			{
			  if(local_toGo[target] < toGo[target])
			    {
			      local_toGo[target] += 1;
			      P[n].Type |= 16;
			    }
			}
		    }
		}

	      for(n = 0; n < NTask; n++)
		{
		  toGo[n] = local_toGo[n];
		  toGoSph[n] = local_toGoSph[n];
		}

	      MPI_Alltoall(toGo, 1, MPI_INT, toGet, 1, MPI_INT, MYMPI_COMM_WORLD);
	      MPI_Alltoall(toGoSph, 1, MPI_INT, toGetSph, 1, MPI_INT, MYMPI_COMM_WORLD);

	      myfree(local_toGoSph);
	      myfree(local_toGo);
	    }
	}
      while(flagsum);

      return 1;
    }
  else
    return 0;
}





/*! This function walks the global top tree in order to establish the
 *  number of leaves it has. These leaves are distributed to different
 *  processors.
 */
void domain_walktoptree(int no)
{
  int i;

  if(topNodes[no].Daughter == -1)
    {
      topNodes[no].Leaf = NTopleaves;
      NTopleaves++;
    }
  else
    {
      for(i = 0; i < 8; i++)
	domain_walktoptree(topNodes[no].Daughter + i);
    }
}


int domain_compare_key(const void *a, const void *b)
{
  if(((struct peano_hilbert_data *) a)->key < (((struct peano_hilbert_data *) b)->key))
    return -1;

  if(((struct peano_hilbert_data *) a)->key > (((struct peano_hilbert_data *) b)->key))
    return +1;

  return 0;
}




int domain_check_for_local_refine(int i, double countlimit, double costlimit)
{
  int j, p, sub, flag = 0;




  if(topNodes[i].Parent >= 0)
    {
      if(topNodes[i].Count > 0.8 * topNodes[topNodes[i].Parent].Count ||
	 topNodes[i].Cost > 0.8 * topNodes[topNodes[i].Parent].Cost)
	flag = 1;
    }

  if((topNodes[i].Count > countlimit || topNodes[i].Cost > costlimit || flag == 1) && topNodes[i].Size >= 8)
    {
      if(topNodes[i].Size >= 8)
	{
	  if((NTopnodes + 8) <= MaxTopNodes)
	    {
	      topNodes[i].Daughter = NTopnodes;

	      for(j = 0; j < 8; j++)
		{
		  sub = topNodes[i].Daughter + j;
		  topNodes[sub].Daughter = -1;
		  topNodes[sub].Parent = i;
		  topNodes[sub].Size = (topNodes[i].Size >> 3);
		  topNodes[sub].StartKey = topNodes[i].StartKey + j * topNodes[sub].Size;
		  topNodes[sub].PIndex = topNodes[i].PIndex;
		  topNodes[sub].Count = 0;
		  topNodes[sub].Cost = 0;
		  /*topNodes[sub].SphCost = 0; *//* added because of arithmetic exception , needed? */

		}

	      NTopnodes += 8;

	      sub = topNodes[i].Daughter;

	      for(p = topNodes[i].PIndex, j = 0; p < topNodes[i].PIndex + topNodes[i].Count; p++)
		{
		  if(j < 7)
		    while(mp[p].key >= topNodes[sub + 1].StartKey)
		      {
			j++;
			sub++;
			topNodes[sub].PIndex = p;
			if(j >= 7)
			  break;
		      }

		  topNodes[sub].Cost += domain_particle_costfactor(mp[p].index);
		  topNodes[sub].Count++;
		}

	      for(j = 0; j < 8; j++)
		{
		  sub = topNodes[i].Daughter + j;

		  if(domain_check_for_local_refine(sub, countlimit, costlimit))
		    return 1;
		}
	    }
	  else
	    return 1;
	}
    }

  return 0;
}


int domain_recursively_combine_topTree(int start, int ncpu)
{
  int i, nleft, nright, errflag = 0;
  int recvTask, ntopnodes_import;
  int master_left, master_right;
  struct local_topnode_data *topNodes_import = 0, *topNodes_temp;

  nleft = ncpu / 2;
  nright = ncpu - nleft;

  if(ncpu > 2)
    {
      errflag += domain_recursively_combine_topTree(start, nleft);
      errflag += domain_recursively_combine_topTree(start + nleft, nright);
    }

  if(ncpu >= 2)
    {
      master_left = start;
      master_right = start + nleft;
      if(master_left == master_right)
	endrun(123);

      if(ThisTask == master_left || ThisTask == master_right)
	{
	  if(ThisTask == master_left)
	    recvTask = master_right;
	  else
	    recvTask = master_left;

	  /* inform each other about the length of the trees */
	  MPI_Sendrecv(&NTopnodes, 1, MPI_INT, recvTask, TAG_GRAV_A,
		       &ntopnodes_import, 1, MPI_INT, recvTask, TAG_GRAV_A, MYMPI_COMM_WORLD,
		       MPI_STATUS_IGNORE);


	  topNodes_import =
	    (struct local_topnode_data *) mymalloc("topNodes_import",
						   IMAX(ntopnodes_import,
							NTopnodes) * sizeof(struct local_topnode_data));

	  /* exchange the trees */
	  MPI_Sendrecv(topNodes,
		       NTopnodes * sizeof(struct local_topnode_data), MPI_BYTE,
		       recvTask, TAG_GRAV_B,
		       topNodes_import,
		       ntopnodes_import * sizeof(struct local_topnode_data), MPI_BYTE,
		       recvTask, TAG_GRAV_B, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

      if(ThisTask == master_left)
	{
	  for(recvTask = master_left + 1; recvTask < master_left + nleft; recvTask++)
	    {
	      MPI_Send(&ntopnodes_import, 1, MPI_INT, recvTask, TAG_GRAV_A, MYMPI_COMM_WORLD);
	      MPI_Send(topNodes_import,
		       ntopnodes_import * sizeof(struct local_topnode_data), MPI_BYTE,
		       recvTask, TAG_GRAV_B, MYMPI_COMM_WORLD);
	    }
	}

      if(ThisTask == master_right)
	{
	  for(recvTask = master_right + 1; recvTask < master_right + nright; recvTask++)
	    {
	      MPI_Send(&ntopnodes_import, 1, MPI_INT, recvTask, TAG_GRAV_A, MYMPI_COMM_WORLD);
	      MPI_Send(topNodes_import,
		       ntopnodes_import * sizeof(struct local_topnode_data), MPI_BYTE,
		       recvTask, TAG_GRAV_B, MYMPI_COMM_WORLD);
	    }
	}

      if(ThisTask > master_left && ThisTask < master_left + nleft)
	{
	  MPI_Recv(&ntopnodes_import, 1, MPI_INT, master_left, TAG_GRAV_A, MYMPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);

	  topNodes_import =
	    (struct local_topnode_data *) mymalloc("topNodes_import",
						   IMAX(ntopnodes_import,
							NTopnodes) * sizeof(struct local_topnode_data));

	  MPI_Recv(topNodes_import,
		   ntopnodes_import * sizeof(struct local_topnode_data), MPI_BYTE,
		   master_left, TAG_GRAV_B, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);

	}


      if(ThisTask > master_right && ThisTask < master_right + nright)
	{
	  MPI_Recv(&ntopnodes_import, 1, MPI_INT, master_right, TAG_GRAV_A, MYMPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);

	  topNodes_import =
	    (struct local_topnode_data *) mymalloc("topNodes_import",
						   IMAX(ntopnodes_import,
							NTopnodes) * sizeof(struct local_topnode_data));

	  MPI_Recv(topNodes_import,
		   ntopnodes_import * sizeof(struct local_topnode_data), MPI_BYTE,
		   master_right, TAG_GRAV_B, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

      if(ThisTask >= master_left && ThisTask < master_left + nleft)
	{
	  /* swap the two trees so that result will be equal on all cpus */

	  topNodes_temp =
	    (struct local_topnode_data *) mymalloc("topNodes_temp",
						   NTopnodes * sizeof(struct local_topnode_data));
	  memcpy(topNodes_temp, topNodes, NTopnodes * sizeof(struct local_topnode_data));
	  memcpy(topNodes, topNodes_import, ntopnodes_import * sizeof(struct local_topnode_data));
	  memcpy(topNodes_import, topNodes_temp, NTopnodes * sizeof(struct local_topnode_data));
	  myfree(topNodes_temp);
	  i = NTopnodes;
	  NTopnodes = ntopnodes_import;
	  ntopnodes_import = i;
	}

      if(ThisTask >= start && ThisTask < start + ncpu)
	{
	  if(errflag == 0)
	    {
	      if((NTopnodes + ntopnodes_import) <= MaxTopNodes)
		{
		  domain_insertnode(topNodes, topNodes_import, 0, 0);
		}
	      else
		{
		  errflag += 1;
		}
	    }

	  myfree(topNodes_import);
	}
    }

  return errflag;
}



/*! This function constructs the global top-level tree node that is used
 *  for the domain decomposition. This is done by considering the string of
 *  Peano-Hilbert keys for all particles, which is recursively chopped off
 *  in pieces of eight segments until each segment holds at most a certain
 *  number of particles.
 */
int domain_determineTopTree(void)
{
  int i, count, j, sub, ngrp;
  int recvTask, sendTask, ntopnodes_import, errflag, errsum;
  struct local_topnode_data *topNodes_import, *topNodes_temp;
  double costlimit, countlimit;
  MPI_Status status;
  int multipledomains = MULTIPLEDOMAINS;

  mp = (struct peano_hilbert_data *) mymalloc("mp", sizeof(struct peano_hilbert_data) * NumPart);

  for(i = 0; i < NumPart; i++)
    {
      Key[i] = peano_hilbert_key((int) ((P[i].Pos[0] - DomainCorner[0]) * DomainFac),
				 (int) ((P[i].Pos[1] - DomainCorner[1]) * DomainFac),
				 (int) ((P[i].Pos[2] - DomainCorner[2]) * DomainFac), BITS_PER_DIMENSION);

    }

  for(i = 0, count = 0; i < NumPart; i++)
    {
      mp[count].key = Key[i];
      mp[count].index = i;
      count++;
    }


  qsort(mp, count, sizeof(struct peano_hilbert_data), domain_compare_key);

  NTopnodes = 1;
  topNodes[0].Daughter = -1;
  topNodes[0].Parent = -1;
  topNodes[0].Size = PEANOCELLS;
  topNodes[0].StartKey = 0;
  topNodes[0].PIndex = 0;
  topNodes[0].Count = count;
  topNodes[0].Cost = gravcost;

  costlimit = totgravcost / (TOPNODEFACTOR * multipledomains * NTask);
  countlimit = totpartcount / (TOPNODEFACTOR * multipledomains * NTask);

  errflag = domain_check_for_local_refine(0, countlimit, costlimit);

  myfree(mp);

  MPI_Allreduce(&errflag, &errsum, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  if(errsum)
    {
      if(ThisTask == 0)
	printf
	  ("We are out of Topnodes. We'll try to repeat with a higher value than All.TopNodeAllocFactor=%g\n",
	   All.TopNodeAllocFactor);
      fflush(stdout);

      return errsum;
    }


  /* we now need to exchange tree parts and combine them as needed */

  if(NTask == (1 << PTask))	/* the following algoritm only works for power of 2 */
    {
      for(ngrp = 1, errflag = 0; ngrp < (1 << PTask); ngrp <<= 1)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      /* inform each other about the length of the trees */
	      MPI_Sendrecv(&NTopnodes, 1, MPI_INT, recvTask, TAG_GRAV_A,
			   &ntopnodes_import, 1, MPI_INT, recvTask, TAG_GRAV_A, MYMPI_COMM_WORLD, &status);


	      topNodes_import =
		(struct local_topnode_data *) mymalloc("topNodes_import",
						       IMAX(ntopnodes_import,
							    NTopnodes) * sizeof(struct local_topnode_data));

	      /* exchange the trees */
	      MPI_Sendrecv(topNodes,
			   NTopnodes * sizeof(struct local_topnode_data), MPI_BYTE,
			   recvTask, TAG_GRAV_B,
			   topNodes_import,
			   ntopnodes_import * sizeof(struct local_topnode_data), MPI_BYTE,
			   recvTask, TAG_GRAV_B, MYMPI_COMM_WORLD, &status);

	      if(sendTask > recvTask)	/* swap the two trees so that result will be equal on all cpus */
		{
		  topNodes_temp =
		    (struct local_topnode_data *) mymalloc("topNodes_temp",
							   NTopnodes * sizeof(struct local_topnode_data));
		  memcpy(topNodes_temp, topNodes, NTopnodes * sizeof(struct local_topnode_data));
		  memcpy(topNodes, topNodes_import, ntopnodes_import * sizeof(struct local_topnode_data));
		  memcpy(topNodes_import, topNodes_temp, NTopnodes * sizeof(struct local_topnode_data));
		  myfree(topNodes_temp);
		  i = NTopnodes;
		  NTopnodes = ntopnodes_import;
		  ntopnodes_import = i;
		}


	      if(errflag == 0)
		{
		  if((NTopnodes + ntopnodes_import) <= MaxTopNodes)
		    {
		      domain_insertnode(topNodes, topNodes_import, 0, 0);
		    }
		  else
		    {
		      errflag = 1;
		    }
		}

	      myfree(topNodes_import);
	    }
	}
    }
  else
    {
      errflag = domain_recursively_combine_topTree(0, NTask);
    }

  MPI_Allreduce(&errflag, &errsum, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

  if(errsum)
    {
      if(ThisTask == 0)
	printf("can't combine trees due to lack of storage. Will try again.\n");
      return errsum;
    }

  /* now let's see whether we should still append more nodes, based on the estimated cumulative cost/count in each cell */

  if(ThisTask == 0)
    printf("Before=%d\n", NTopnodes);

  for(i = 0, errflag = 0; i < NTopnodes; i++)
    {
      if(topNodes[i].Daughter < 0)
	if(topNodes[i].Count > countlimit || topNodes[i].Cost > costlimit)	/* ok, let's add nodes if we can */
	  if(topNodes[i].Size > 1)
	    {
	      if((NTopnodes + 8) <= MaxTopNodes)
		{
		  topNodes[i].Daughter = NTopnodes;

		  for(j = 0; j < 8; j++)
		    {
		      sub = topNodes[i].Daughter + j;
		      topNodes[sub].Size = (topNodes[i].Size >> 3);
		      topNodes[sub].Count = topNodes[i].Count / 8;
		      topNodes[sub].Cost = topNodes[i].Cost / 8;
		      topNodes[sub].Daughter = -1;
		      topNodes[sub].Parent = i;
		      topNodes[sub].StartKey = topNodes[i].StartKey + j * topNodes[sub].Size;
		    }

		  NTopnodes += 8;
		}
	      else
		{
		  errflag = 1;
		  break;
		}
	    }
    }

  MPI_Allreduce(&errflag, &errsum, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  if(errsum)
    return errsum;

  if(ThisTask == 0)
    printf("After=%d\n", NTopnodes);

  /* count toplevel leaves */
  domain_sumCost();

  if(NTopleaves < multipledomains * NTask)
    endrun(112);

  return 0;
}



void domain_sumCost(void)
{
  int i, n, no;
  float *local_domainWork;
  float *local_domainWorkSph;
  int *local_domainCount;
  int *local_domainCountSph;


  local_domainWork = (float *) mymalloc("local_domainWork", NTopnodes * sizeof(float));
  local_domainWorkSph = (float *) mymalloc("local_domainWorkSph", NTopnodes * sizeof(float));
  local_domainCount = (int *) mymalloc("local_domainCount", NTopnodes * sizeof(int));
  local_domainCountSph = (int *) mymalloc("local_domainCountSph", NTopnodes * sizeof(int));



  NTopleaves = 0;
  domain_walktoptree(0);

  for(i = 0; i < NTopleaves; i++)
    {
      local_domainWork[i] = 0;
      local_domainWorkSph[i] = 0;
      local_domainCount[i] = 0;
      local_domainCountSph[i] = 0;
    }

  if(ThisTask == 0)
    printf("NTopleaves= %d  NTopnodes=%d (space for %d)\n", NTopleaves, NTopnodes, MaxTopNodes);

  for(n = 0; n < NumPart; n++)
    {

      no = 0;

      while(topNodes[no].Daughter >= 0)
	no = topNodes[no].Daughter + (Key[n] - topNodes[no].StartKey) / (topNodes[no].Size >> 3);

      no = topNodes[no].Leaf;

      local_domainWork[no] += (float) domain_particle_costfactor(n);

      local_domainCount[no] += 1;
      if(P[n].Type == 0)
	{
	  local_domainCountSph[no] += 1;
	  if(TimeBinActive[P[n].TimeBin] || UseAllParticles)
	    local_domainWorkSph[no] += 1.0;
	}
    }

  MPI_Allreduce(local_domainWork, domainWork, NTopleaves, MPI_FLOAT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(local_domainWorkSph, domainWorkSph, NTopleaves, MPI_FLOAT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(local_domainCount, domainCount, NTopleaves, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
  MPI_Allreduce(local_domainCountSph, domainCountSph, NTopleaves, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);

  myfree(local_domainCountSph);
  myfree(local_domainCount);
  myfree(local_domainWorkSph);
  myfree(local_domainWork);


}


/*! This routine finds the extent of the global domain grid.
 */
void domain_findExtent(void)
{
  int i, j;
  double len, xmin[3], xmax[3], xmin_glob[3], xmax_glob[3];


  /* determine local extension */
  for(j = 0; j < 3; j++)
    {
      xmin[j] = MAX_REAL_NUMBER;
      xmax[j] = -MAX_REAL_NUMBER;
    }

    for(i = 0; i < NumPart; i++)
      {

	for(j = 0; j < 3; j++)
	  {
	    if(xmin[j] > P[i].Pos[j])
	      xmin[j] = P[i].Pos[j];

	    if(xmax[j] < P[i].Pos[j])
	      xmax[j] = P[i].Pos[j];
	  }
      }

  MPI_Allreduce(xmin, xmin_glob, 3, MPI_DOUBLE, MPI_MIN, MYMPI_COMM_WORLD);
  MPI_Allreduce(xmax, xmax_glob, 3, MPI_DOUBLE, MPI_MAX, MYMPI_COMM_WORLD);

  len = 0;
  for(j = 0; j < 3; j++)
    if(xmax_glob[j] - xmin_glob[j] > len)
      len = xmax_glob[j] - xmin_glob[j];

  len *= 1.001;

  for(j = 0; j < 3; j++)
    {
      DomainCenter[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]);
      DomainCorner[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]) - 0.5 * len;
    }


  DomainLen = len;
  DomainFac = 1.0 / len * (((peanokey) 1) << (BITS_PER_DIMENSION));
}




void domain_add_cost(struct local_topnode_data *treeA, int noA, long long count, double cost, double sphcost)
{
  int i, sub;
  long long countA, countB;

  countB = count / 8;
  countA = count - 7 * countB;

  cost = cost / 8;
  sphcost = sphcost / 8;

  for(i = 0; i < 8; i++)
    {
      sub = treeA[noA].Daughter + i;

      if(i == 0)
	count = countA;
      else
	count = countB;

      treeA[sub].Count += count;
      treeA[sub].Cost += cost;
      treeA[sub].SphCost += sphcost;

      if(treeA[sub].Daughter >= 0)
	domain_add_cost(treeA, sub, count, cost, sphcost);
    }
}


void domain_insertnode(struct local_topnode_data *treeA, struct local_topnode_data *treeB, int noA, int noB)
{
  int j, sub;
  long long count, countA, countB;
  double cost, costA, costB;

  if(treeB[noB].Size < treeA[noA].Size)
    {
      if(treeA[noA].Daughter < 0)
	{
	  if((NTopnodes + 8) <= MaxTopNodes)
	    {
	      count = treeA[noA].Count - treeB[treeB[noB].Parent].Count;
	      countB = count / 8;
	      countA = count - 7 * countB;

	      cost = treeA[noA].Cost - treeB[treeB[noB].Parent].Cost;
	      costB = cost / 8;
	      costA = cost - 7 * costB;

	      treeA[noA].Daughter = NTopnodes;
	      for(j = 0; j < 8; j++)
		{
		  if(j == 0)
		    {
		      count = countA;
		      cost = costA;
		    }
		  else
		    {
		      count = countB;
		      cost = costB;
		    }

		  sub = treeA[noA].Daughter + j;
		  topNodes[sub].Size = (treeA[noA].Size >> 3);
		  topNodes[sub].Count = count;
		  topNodes[sub].Cost = cost;
		  topNodes[sub].Daughter = -1;
		  topNodes[sub].Parent = noA;
		  topNodes[sub].StartKey = treeA[noA].StartKey + j * treeA[sub].Size;
		}
	      NTopnodes += 8;
	    }
	  else
	    endrun(88);
	}

      sub = treeA[noA].Daughter + (treeB[noB].StartKey - treeA[noA].StartKey) / (treeA[noA].Size >> 3);
      domain_insertnode(treeA, treeB, sub, noB);
    }
  else if(treeB[noB].Size == treeA[noA].Size)
    {
      treeA[noA].Count += treeB[noB].Count;
      treeA[noA].Cost += treeB[noB].Cost;

      if(treeB[noB].Daughter >= 0)
	{
	  for(j = 0; j < 8; j++)
	    {
	      sub = treeB[noB].Daughter + j;
	      domain_insertnode(treeA, treeB, noA, sub);
	    }
	}
      else
	{
	  if(treeA[noA].Daughter >= 0)
	    domain_add_cost(treeA, noA, treeB[noB].Count, treeB[noB].Cost, treeB[noB].SphCost);
	}
    }
  else
    endrun(89);
}



static void msort_domain_with_tmp(struct peano_hilbert_data *b, size_t n, struct peano_hilbert_data *t)
{
  struct peano_hilbert_data *tmp;
  struct peano_hilbert_data *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  msort_domain_with_tmp(b1, n1, t);
  msort_domain_with_tmp(b2, n2, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(b1->key <= b2->key)
	{
	  --n1;
	  *tmp++ = *b1++;
	}
      else
	{
	  --n2;
	  *tmp++ = *b2++;
	}
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * sizeof(struct peano_hilbert_data));

  memcpy(b, t, (n - n2) * sizeof(struct peano_hilbert_data));
}

void mysort_domain(void *b, size_t n, size_t s)
{
  const size_t size = n * s;
  struct peano_hilbert_data *tmp;

  tmp = (struct peano_hilbert_data *) mymalloc("tmp", size);

  msort_domain_with_tmp((struct peano_hilbert_data *) b, n, tmp);

  myfree(tmp);
}


/* GM note: this was in predict.c */
/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef PERIODIC
void do_box_wrapping(void)
{
  int i, j;
  double boxsize[3];

  for(j = 0; j < 3; j++)
    boxsize[j] = All.BoxSize;

#ifdef LONG_X
  boxsize[0] *= LONG_X;
#endif
#ifdef LONG_Y
  boxsize[1] *= LONG_Y;
#endif
#ifdef LONG_Z
  boxsize[2] *= LONG_Z;
#endif

  for(i = 0; i < NumPart; i++)
      for(j = 0; j < 3; j++)
	{
	  while(P[i].Pos[j] < 0)
	    P[i].Pos[j] += boxsize[j];

	  while(P[i].Pos[j] >= boxsize[j])
	    P[i].Pos[j] -= boxsize[j];
	}
}
#endif
