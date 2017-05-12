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
#include <time.h>

#include "allvars.h"
#include "proto.h"



/*! \file forcetree.c
 *  \brief gravitational tree and code for Ewald correction
 *
 *  This file contains the computation of the gravitational force by means
 *  of a tree. The type of tree implemented is a geometrical oct-tree,
 *  starting from a cube encompassing all particles. This cube is
 *  automatically found in the domain decomposition, which also splits up
 *  the global "top-level" tree along node boundaries, moving the particles
 *  of different parts of the tree to separate processors. Tree nodes can
 *  be dynamically updated in drift/kick operations to avoid having to
 *  reconstruct the tree every timestep.
 */

/*! auxialiary variable used to set-up non-recursive walk */
static int last;



/*! length of lock-up table for short-range force kernel in TreePM algorithm */
#define NTAB 1000
/*! variables for short-range lookup table */
static float shortrange_table[NTAB], shortrange_table_potential[NTAB];
/*! toggles after first tree-memory allocation, has only influence on log-files */
static int first_flag = 0;

static int tree_allocated_flag = 0;



/*! This function is a driver routine for constructing the gravitational
 *  oct-tree, which is done by calling a small number of other functions.
 */
int force_treebuild(int npart, struct unbind_data *mp)
{
  int flag;


  do
    {
      Numnodestree = force_treebuild_single(npart, mp);

      MPI_Allreduce(&Numnodestree, &flag, 1, MPI_INT, MPI_MIN, MYMPI_COMM_WORLD);
      if(flag == -1)
	{
	  force_treefree();

	  if(ThisTask == 0)
	    printf("Increasing TreeAllocFactor=%g", All.TreeAllocFactor);

	  All.TreeAllocFactor *= 1.15;

	  if(ThisTask == 0)
	    printf("new value=%g\n", All.TreeAllocFactor);

	  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart) + NTopnodes, All.MaxPart);
	}
    }
  while(flag == -1);

  force_flag_localnodes();

  force_exchange_pseudodata();

  force_treeupdate_pseudos(All.MaxPart);

  TimeOfLastTreeConstruction = All.Time;

  return Numnodestree;
}



/*! Constructs the gravitational oct-tree.  
 *
 *  The index convention for accessing tree nodes is the following: the
 *  indices 0...NumPart-1 reference single particles, the indices
 *  All.MaxPart.... All.MaxPart+nodes-1 reference tree nodes. `Nodes_base'
 *  points to the first tree node, while `nodes' is shifted such that
 *  nodes[All.MaxPart] gives the first tree node. Finally, node indices
 *  with values 'All.MaxPart + MaxNodes' and larger indicate "pseudo
 *  particles", i.e. multipole moments of top-level nodes that lie on
 *  different CPUs. If such a node needs to be opened, the corresponding
 *  particle must be exported to that CPU. The 'Extnodes' structure
 *  parallels that of 'Nodes'. Its information is only needed for the SPH
 *  part of the computation. (The data is split onto these two structures
 *  as a tuning measure.  If it is merged into 'Nodes' a somewhat bigger
 *  size of the nodes also for gravity would result, which would reduce
 *  cache utilization slightly.
 */
int force_treebuild_single(int npart, struct unbind_data *mp)
{
  int i, j, k, subnode = 0, shift, parent, numnodes, rep;
  int nfree, th, nn, no;
  struct NODE *nfreep;
  MyFloat lenhalf;
  peanokey key, morton, th_key, *morton_list;


  /* create an empty root node  */
  nfree = All.MaxPart;		/* index of first free node */
  nfreep = &Nodes[nfree];	/* select first node */

  nfreep->len = DomainLen;
  for(j = 0; j < 3; j++)
    nfreep->center[j] = DomainCenter[j];
  for(j = 0; j < 8; j++)
    nfreep->u.suns[j] = -1;


  numnodes = 1;
  nfreep++;
  nfree++;

  /* create a set of empty nodes corresponding to the top-level domain
   * grid. We need to generate these nodes first to make sure that we have a
   * complete top-level tree which allows the easy insertion of the
   * pseudo-particles at the right place 
   */

  force_create_empty_nodes(All.MaxPart, 0, 1, 0, 0, 0, &numnodes, &nfree);

  /* if a high-resolution region in a global tree is used, we need to generate
   * an additional set empty nodes to make sure that we have a complete
   * top-level tree for the high-resolution inset
   */

  nfreep = &Nodes[nfree];
  parent = -1;			/* note: will not be used below before it is changed */

  morton_list = (peanokey *) mymalloc("morton_list", NumPart * sizeof(peanokey));

  /* now we insert all particles */
  for(k = 0; k < npart; k++)
    {
      if(mp)
	i = mp[k].index;
      else
	i = k;


      rep = 0;

      key = peano_and_morton_key((int) ((P[i].Pos[0] - DomainCorner[0]) * DomainFac),
				 (int) ((P[i].Pos[1] - DomainCorner[1]) * DomainFac),
				 (int) ((P[i].Pos[2] - DomainCorner[2]) * DomainFac), BITS_PER_DIMENSION,
				 &morton);
      morton_list[i] = morton;

      shift = 3 * (BITS_PER_DIMENSION - 1);

      no = 0;
      while(TopNodes[no].Daughter >= 0)
	{
	  no = TopNodes[no].Daughter + (key - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);
	  shift -= 3;
	  rep++;
	}

      no = TopNodes[no].Leaf;
      th = DomainNodeIndex[no];

      while(1)
	{
	  if(th >= All.MaxPart)	/* we are dealing with an internal node */
	    {
	      if(shift >= 0)
		{
		  subnode = ((morton >> shift) & 7);
		}
	      else
		{
		  subnode = 0;
		  if(P[i].Pos[0] > Nodes[th].center[0])
		    subnode += 1;
		  if(P[i].Pos[1] > Nodes[th].center[1])
		    subnode += 2;
		  if(P[i].Pos[2] > Nodes[th].center[2])
		    subnode += 4;
		}

	      if(Nodes[th].len < 1.0e-3 * All.ForceSoftening[P[i].Type])
		{
		  /* seems like we're dealing with particles at identical (or extremely close)
		   * locations. Randomize subnode index to allow tree construction. Note: Multipole moments
		   * of tree are still correct, but this will only happen well below gravitational softening
		   * length-scale anyway.
		   */
		  subnode = (int) (8.0 * get_random_number((P[i].ID + rep) % (RNDTABLE + (rep & 3))));

		  if(subnode >= 8)
		    subnode = 7;
		}

	      nn = Nodes[th].u.suns[subnode];

	      shift -= 3;

	      if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
		{
		  parent = th;
		  th = nn;
		  rep++;
		}
	      else
		{
		  /* here we have found an empty slot where we can attach
		   * the new particle as a leaf.
		   */
		  Nodes[th].u.suns[subnode] = i;
		  break;	/* done for this particle */
		}
	    }
	  else
	    {
	      /* We try to insert into a leaf with a single particle.  Need
	       * to generate a new internal node at this point.
	       */
	      Nodes[parent].u.suns[subnode] = nfree;

	      nfreep->len = 0.5 * Nodes[parent].len;
	      lenhalf = 0.25 * Nodes[parent].len;

	      if(subnode & 1)
		nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
	      else
		nfreep->center[0] = Nodes[parent].center[0] - lenhalf;

	      if(subnode & 2)
		nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
	      else
		nfreep->center[1] = Nodes[parent].center[1] - lenhalf;

	      if(subnode & 4)
		nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
	      else
		nfreep->center[2] = Nodes[parent].center[2] - lenhalf;

	      nfreep->u.suns[0] = -1;
	      nfreep->u.suns[1] = -1;
	      nfreep->u.suns[2] = -1;
	      nfreep->u.suns[3] = -1;
	      nfreep->u.suns[4] = -1;
	      nfreep->u.suns[5] = -1;
	      nfreep->u.suns[6] = -1;
	      nfreep->u.suns[7] = -1;

	      if(shift >= 0)
		{
		  th_key = morton_list[th];
		  subnode = ((th_key >> shift) & 7);
		}
	      else
		{
		  subnode = 0;
		  if(P[th].Pos[0] > nfreep->center[0])
		    subnode += 1;
		  if(P[th].Pos[1] > nfreep->center[1])
		    subnode += 2;
		  if(P[th].Pos[2] > nfreep->center[2])
		    subnode += 4;
		}

	      if(nfreep->len < 1.0e-3 * All.ForceSoftening[P[th].Type])
		{
		  /* seems like we're dealing with particles at identical (or extremely close)
		   * locations. Randomize subnode index to allow tree construction. Note: Multipole moments
		   * of tree are still correct, but this will only happen well below gravitational softening
		   * length-scale anyway.
		   */
		  subnode = (int) (8.0 * get_random_number((P[th].ID + rep) % (RNDTABLE + (rep & 3))));

		  if(subnode >= 8)
		    subnode = 7;
		}
	      nfreep->u.suns[subnode] = th;

	      th = nfree;	/* resume trying to insert the new particle at
				 * the newly created internal node
				 */

	      numnodes++;
	      nfree++;
	      nfreep++;

	      if((numnodes) >= MaxNodes)
		{
		  printf("task %d: maximum number %d of tree-nodes reached for particle %d.\n", ThisTask,
			 MaxNodes, i);

		  if(All.TreeAllocFactor > 5.0)
		    {
		      printf
			("task %d: looks like a serious problem for particle %d, stopping with particle dump.\n",
			 ThisTask, i);
		      dump_particles();
		      endrun(1);
		    }
		  else
		    {
		      myfree(morton_list);
		      return -1;
		    }
		}
	    }
	}
    }

  myfree(morton_list);


  /* insert the pseudo particles that represent the mass distribution of other domains */
  force_insert_pseudo_particles();





  /* now compute the multipole moments recursively */
  last = -1;

  force_update_node_recursive(All.MaxPart, -1, -1, 0);

  if(last >= All.MaxPart)
    {
      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
	Nextnode[last - MaxNodes] = -1;
      else
	Nodes[last].u.d.nextnode = -1;
    }
  else
    Nextnode[last] = -1;

  return numnodes;
}



/*! This function recursively creates a set of empty tree nodes which
 *  corresponds to the top-level tree for the domain grid. This is done to
 *  ensure that this top-level tree is always "complete" so that we can easily
 *  associate the pseudo-particles of other CPUs with tree-nodes at a given
 *  level in the tree, even when the particle population is so sparse that
 *  some of these nodes are actually empty.
*/
void force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, int *nodecount,
			      int *nextfree)
{
  int i, j, k, n, sub, count;
  MyFloat lenhalf;

  if(TopNodes[topnode].Daughter >= 0)
    {
      for(i = 0; i < 2; i++)
	for(j = 0; j < 2; j++)
	  for(k = 0; k < 2; k++)
	    {
	      sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

	      count = i + 2 * j + 4 * k;

	      Nodes[no].u.suns[count] = *nextfree;

	      lenhalf = 0.25 * Nodes[no].len;
	      Nodes[*nextfree].len = 0.5 * Nodes[no].len;
	      Nodes[*nextfree].center[0] = Nodes[no].center[0] + (2 * i - 1) * lenhalf;
	      Nodes[*nextfree].center[1] = Nodes[no].center[1] + (2 * j - 1) * lenhalf;
	      Nodes[*nextfree].center[2] = Nodes[no].center[2] + (2 * k - 1) * lenhalf;

	      for(n = 0; n < 8; n++)
		Nodes[*nextfree].u.suns[n] = -1;

	      if(TopNodes[TopNodes[topnode].Daughter + sub].Daughter == -1)
		DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf] = *nextfree;

	      *nextfree = *nextfree + 1;
	      *nodecount = *nodecount + 1;

	      if((*nodecount) >= MaxNodes || (*nodecount) >= MaxTopNodes)
		{
		  printf("task %d: maximum number MaxNodes=%d of tree-nodes reached."
			 "MaxTopNodes=%d NTopnodes=%d NTopleaves=%d nodecount=%d\n",
			 ThisTask, MaxNodes, MaxTopNodes, NTopnodes, NTopleaves, *nodecount);
		  printf("in create empty nodes\n");
		  dump_particles();
		  endrun(11);
		}

	      force_create_empty_nodes(*nextfree - 1, TopNodes[topnode].Daughter + sub,
				       bits + 1, 2 * x + i, 2 * y + j, 2 * z + k, nodecount, nextfree);
	    }
    }
}



/*! this function inserts pseudo-particles which will represent the mass
 *  distribution of the other CPUs. Initially, the mass of the
 *  pseudo-particles is set to zero, and their coordinate is set to the
 *  center of the domain-cell they correspond to. These quantities will be
 *  updated later on.
 */
void force_insert_pseudo_particles(void)
{
  int i, index;

  for(i = 0; i < NTopleaves; i++)
    {
      index = DomainNodeIndex[i];

      if(DomainTask[i] != ThisTask)
	Nodes[index].u.suns[0] = All.MaxPart + MaxNodes + i;
    }
}


/*! this routine determines the multipole moments for a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Nodes[] structure in the sequence of this tree-walk.
 *
 *  Note that the bitflags-variable for each node is used to store in the
 *  lowest bits some special information: Bit 0 flags whether the node
 *  belongs to the top-level tree corresponding to the domain
 *  decomposition, while Bit 1 signals whether the top-level node is
 *  dependent on local mass.
 * 
 *  If UNEQUALSOFTENINGS is set, bits 2-4 give the particle type with
 *  the maximum softening among the particles in the node, and bit 5
 *  flags whether the node contains any particles with lower softening
 *  than that.  
 */
void force_update_node_recursive(int no, int sib, int father, float LocSoft)
{
  int j, jj, k, p, pp, nextsib, suns[8], count_particles, multiple_flag;
  MyFloat hmax, vmax, v;
  MyFloat divVmax;
  MyFloat s[3], vs[3], mass;
  struct particle_data *pa;


#ifdef UNEQUALSOFTENINGS
  int maxsofttype, current_maxsofttype, diffsoftflag;
#endif


  if(no >= All.MaxPart && no < All.MaxPart + MaxNodes)	/* internal node */
    {
      for(j = 0; j < 8; j++)
	suns[j] = Nodes[no].u.suns[j];	/* this "backup" is necessary because the nextnode entry will overwrite one element (union!) */
      if(last >= 0)
	{
	  if(last >= All.MaxPart)
	    {
	      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
		Nextnode[last - MaxNodes] = no;
	      else
		Nodes[last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[last] = no;
	}

      last = no;

      mass = 0;
      s[0] = 0;
      s[1] = 0;
      s[2] = 0;
      vs[0] = 0;
      vs[1] = 0;
      vs[2] = 0;
      hmax = 0;
      vmax = 0;
      count_particles = 0;

#ifdef UNEQUALSOFTENINGS
      maxsofttype = 7;
      diffsoftflag = 0;
#endif


      for(j = 0; j < 8; j++)
	{
	  if((p = suns[j]) >= 0)
	    {
	      /* check if we have a sibling on the same level */
	      for(jj = j + 1; jj < 8; jj++)
		if((pp = suns[jj]) >= 0)
		  break;

	      if(jj < 8)	/* yes, we do */
		nextsib = pp;
	      else
		nextsib = sib;

	      force_update_node_recursive(p, nextsib, no, LocSoft);

	      if(p >= All.MaxPart)	/* an internal node or pseudo particle */
		{
		  if(p >= All.MaxPart + MaxNodes)	/* a pseudo particle */
		    {
		      /* nothing to be done here because the mass of the
		       * pseudo-particle is still zero. This will be changed
		       * later.
		       */
		    }
		  else
		    {
		      mass += (Nodes[p].u.d.mass);
		      s[0] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[0]);
		      s[1] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[1]);
		      s[2] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[2]);
		      vs[0] += (Nodes[p].u.d.mass * Extnodes[p].vs[0]);
		      vs[1] += (Nodes[p].u.d.mass * Extnodes[p].vs[1]);
		      vs[2] += (Nodes[p].u.d.mass * Extnodes[p].vs[2]);
		      if(Nodes[p].u.d.mass > 0)
			{
			  if(Nodes[p].u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES))
			    count_particles += 2;
			  else
			    count_particles++;
			}

		      if(Extnodes[p].hmax > hmax)
			hmax = Extnodes[p].hmax;

		      if(Extnodes[p].vmax > vmax)
			vmax = Extnodes[p].vmax;

		      if(Extnodes[p].divVmax > divVmax)
			divVmax = Extnodes[p].divVmax;


#ifdef UNEQUALSOFTENINGS
		      diffsoftflag |= maskout_different_softening_flag(Nodes[p].u.d.bitflags);

		      if(maxsofttype == 7)
			maxsofttype = extract_max_softening_type(Nodes[p].u.d.bitflags);
		      else
			{
			  current_maxsofttype = extract_max_softening_type(Nodes[p].u.d.bitflags);
			  if(current_maxsofttype != 7)
			    {
			      if(All.ForceSoftening[current_maxsofttype] > All.ForceSoftening[maxsofttype])
				{
				  maxsofttype = current_maxsofttype;
				  diffsoftflag = (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE);
				}
			      else
				{
				  if(All.ForceSoftening[current_maxsofttype] <
				     All.ForceSoftening[maxsofttype])
				    diffsoftflag = (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE);
				}
			    }
			}

		    }
		}
	      else		/* a particle */
		{
		  count_particles++;

		  pa = &P[p];

		  mass += (pa->Mass);
		  s[0] += (pa->Mass * pa->Pos[0]);
		  s[1] += (pa->Mass * pa->Pos[1]);
		  s[2] += (pa->Mass * pa->Pos[2]);

		  vs[0] += (pa->Mass * pa->Vel[0]);
		  vs[1] += (pa->Mass * pa->Vel[1]);
		  vs[2] += (pa->Mass * pa->Vel[2]);

#endif

		  if(pa->Type == 0)
		    {
		      if(PPP[p].Hsml > hmax)
			hmax = PPP[p].Hsml;

		      if(SphP[p].v.DivVel > divVmax)
			divVmax = SphP[p].v.DivVel;
		    }


		  for(k = 0; k < 3; k++)
		    if((v = fabs(pa->Vel[k])) > vmax)
		      vmax = v;


#ifdef UNEQUALSOFTENINGS
		  if(maxsofttype == 7)
		    maxsofttype = pa->Type;
		  else
		    {
		      if(All.ForceSoftening[pa->Type] > All.ForceSoftening[maxsofttype])
			{
			  maxsofttype = pa->Type;
			  diffsoftflag = (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE);
			}
		      else
			{
			  if(All.ForceSoftening[pa->Type] < All.ForceSoftening[maxsofttype])
			    diffsoftflag = (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE);
			}
		    }
#endif
		}
	    }
	}


      if(mass)
	{
	  s[0] /= mass;
	  s[1] /= mass;
	  s[2] /= mass;
	  vs[0] /= mass;
	  vs[1] /= mass;
	  vs[2] /= mass;
	}
      else
	{
	  s[0] = Nodes[no].center[0];
	  s[1] = Nodes[no].center[1];
	  s[2] = Nodes[no].center[2];
	  vs[0] = 0;
	  vs[1] = 0;
	  vs[2] = 0;
	}



      Nodes[no].Ti_current = All.Ti_Current;
      Nodes[no].u.d.mass = mass;
      Nodes[no].u.d.s[0] = s[0];
      Nodes[no].u.d.s[1] = s[1];
      Nodes[no].u.d.s[2] = s[2];
      Nodes[no].GravCost = 0;

      Extnodes[no].Ti_lastkicked = All.Ti_Current;
      Extnodes[no].Flag = GlobFlag;
      Extnodes[no].vs[0] = vs[0];
      Extnodes[no].vs[1] = vs[1];
      Extnodes[no].vs[2] = vs[2];
      Extnodes[no].hmax = hmax;
      Extnodes[no].vmax = vmax;
      Extnodes[no].dp[0] = 0;
      Extnodes[no].dp[1] = 0;
      Extnodes[no].dp[2] = 0;


      if(count_particles > 1)	/* this flags that the node represents more than one particle */
	multiple_flag = (1 << BITFLAG_MULTIPLEPARTICLES);
      else
	multiple_flag = 0;

      Nodes[no].u.d.bitflags = multiple_flag;

#ifdef UNEQUALSOFTENINGS
      Nodes[no].u.d.bitflags |= diffsoftflag + (maxsofttype << BITFLAG_MAX_SOFTENING_TYPE);
#endif
      Nodes[no].u.d.sibling = sib;
      Nodes[no].u.d.father = father;
    }
  else				/* single particle or pseudo particle */
    {
      if(last >= 0)
	{
	  if(last >= All.MaxPart)
	    {
	      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
		Nextnode[last - MaxNodes] = no;
	      else
		Nodes[last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[last] = no;
	}

      last = no;

      if(no < All.MaxPart)	/* only set it for single particles */
	Father[no] = father;

    }
}





/*! This function communicates the values of the multipole moments of the
 *  top-level tree-nodes of the domain grid.  This data can then be used to
 *  update the pseudo-particles on each CPU accordingly.
 */
void force_exchange_pseudodata(void)
{
  int i, no, m, ta, recvTask;
  int *recvcounts, *recvoffset;
  struct DomainNODE
  {
    MyFloat s[3];
    MyFloat vs[3];
    MyFloat mass;
    MyFloat hmax;
    MyFloat vmax;

    unsigned int bitflags;
  }
   *DomainMoment;


  DomainMoment = (struct DomainNODE *) mymalloc("DomainMoment", NTopleaves * sizeof(struct DomainNODE));

  for(m = 0; m < MULTIPLEDOMAINS; m++)
    for(i = DomainStartList[ThisTask * MULTIPLEDOMAINS + m];
	i <= DomainEndList[ThisTask * MULTIPLEDOMAINS + m]; i++)
      {
	no = DomainNodeIndex[i];

	/* read out the multipole moments from the local base cells */
	DomainMoment[i].s[0] = Nodes[no].u.d.s[0];
	DomainMoment[i].s[1] = Nodes[no].u.d.s[1];
	DomainMoment[i].s[2] = Nodes[no].u.d.s[2];
	DomainMoment[i].vs[0] = Extnodes[no].vs[0];
	DomainMoment[i].vs[1] = Extnodes[no].vs[1];
	DomainMoment[i].vs[2] = Extnodes[no].vs[2];
	DomainMoment[i].mass = Nodes[no].u.d.mass;
	DomainMoment[i].hmax = Extnodes[no].hmax;
	DomainMoment[i].vmax = Extnodes[no].vmax;
	DomainMoment[i].bitflags = Nodes[no].u.d.bitflags;


      }

  /* share the pseudo-particle data accross CPUs */

  recvcounts = (int *) mymalloc("recvcounts", sizeof(int) * NTask);
  recvoffset = (int *) mymalloc("recvoffset", sizeof(int) * NTask);

  for(m = 0; m < MULTIPLEDOMAINS; m++)
    {
      for(recvTask = 0; recvTask < NTask; recvTask++)
	{
	  recvcounts[recvTask] =
	    (DomainEndList[recvTask * MULTIPLEDOMAINS + m] - DomainStartList[recvTask * MULTIPLEDOMAINS + m] +
	     1) * sizeof(struct DomainNODE);
	  recvoffset[recvTask] = DomainStartList[recvTask * MULTIPLEDOMAINS + m] * sizeof(struct DomainNODE);
	}

      MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &DomainMoment[0], recvcounts, recvoffset, MPI_BYTE,
		     MYMPI_COMM_WORLD);
    }

  myfree(recvoffset);
  myfree(recvcounts);


  for(ta = 0; ta < NTask; ta++)
    if(ta != ThisTask)
      for(m = 0; m < MULTIPLEDOMAINS; m++)
	for(i = DomainStartList[ta * MULTIPLEDOMAINS + m]; i <= DomainEndList[ta * MULTIPLEDOMAINS + m]; i++)
	  {
	    no = DomainNodeIndex[i];

	    Nodes[no].u.d.s[0] = DomainMoment[i].s[0];
	    Nodes[no].u.d.s[1] = DomainMoment[i].s[1];
	    Nodes[no].u.d.s[2] = DomainMoment[i].s[2];
	    Extnodes[no].vs[0] = DomainMoment[i].vs[0];
	    Extnodes[no].vs[1] = DomainMoment[i].vs[1];
	    Extnodes[no].vs[2] = DomainMoment[i].vs[2];
	    Nodes[no].u.d.mass = DomainMoment[i].mass;
	    Extnodes[no].hmax = DomainMoment[i].hmax;
	    Extnodes[no].vmax = DomainMoment[i].vmax;

	    Nodes[no].u.d.bitflags =
	      (Nodes[no].u.d.bitflags & (~BITFLAG_MASK)) | (DomainMoment[i].bitflags & BITFLAG_MASK);

	  }

  myfree(DomainMoment);
}



/*! This function updates the top-level tree after the multipole moments of
 *  the pseudo-particles have been updated.
 */
void force_treeupdate_pseudos(int no)
{
  int j, p, count_particles, multiple_flag;
  MyFloat hmax, vmax;
  MyFloat s[3], vs[3], mass;



#ifdef UNEQUALSOFTENINGS
  int maxsofttype, diffsoftflag, current_maxsofttype;
#endif



  mass = 0;
  s[0] = 0;
  s[1] = 0;
  s[2] = 0;
  vs[0] = 0;
  vs[1] = 0;
  vs[2] = 0;
  hmax = 0;
  vmax = 0;

  count_particles = 0;
#ifdef UNEQUALSOFTENINGS
  maxsofttype = 7;
  diffsoftflag = 0;
#endif


  p = Nodes[no].u.d.nextnode;

  for(j = 0; j < 8; j++)	/* since we are dealing with top-level nodes, we now that there are 8 consecutive daughter nodes */
    {
      if(p >= All.MaxPart && p < All.MaxPart + MaxNodes)	/* internal node */
	{
	  if(Nodes[p].u.d.bitflags & (1 << BITFLAG_INTERNAL_TOPLEVEL))
	    force_treeupdate_pseudos(p);

	  mass += (Nodes[p].u.d.mass);
	  s[0] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[0]);
	  s[1] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[1]);
	  s[2] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[2]);
	  vs[0] += (Nodes[p].u.d.mass * Extnodes[p].vs[0]);
	  vs[1] += (Nodes[p].u.d.mass * Extnodes[p].vs[1]);
	  vs[2] += (Nodes[p].u.d.mass * Extnodes[p].vs[2]);

	  if(Extnodes[p].hmax > hmax)
	    hmax = Extnodes[p].hmax;
	  if(Extnodes[p].vmax > vmax)
	    vmax = Extnodes[p].vmax;

	  if(Nodes[p].u.d.mass > 0)
	    {
	      if(Nodes[p].u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES))
		count_particles += 2;
	      else
		count_particles++;
	    }

#ifdef UNEQUALSOFTENINGS
	  diffsoftflag |= maskout_different_softening_flag(Nodes[p].u.d.bitflags);

	  if(maxsofttype == 7)
	    maxsofttype = extract_max_softening_type(Nodes[p].u.d.bitflags);
	  else
	    {
	      current_maxsofttype = extract_max_softening_type(Nodes[p].u.d.bitflags);
	      if(current_maxsofttype != 7)
		{
		  if(All.ForceSoftening[current_maxsofttype] > All.ForceSoftening[maxsofttype])
		    {
		      maxsofttype = current_maxsofttype;
		      diffsoftflag = (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE);
		    }
		  else
		    {
		      if(All.ForceSoftening[current_maxsofttype] < All.ForceSoftening[maxsofttype])
			diffsoftflag = (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE);
		    }
		}
	    }
#endif

	}
      else
	endrun(6767);		/* may not happen */

      p = Nodes[p].u.d.sibling;
    }

  if(mass)
    {
      s[0] /= mass;
      s[1] /= mass;
      s[2] /= mass;
      vs[0] /= mass;
      vs[1] /= mass;
      vs[2] /= mass;
    }
  else
    {
      s[0] = Nodes[no].center[0];
      s[1] = Nodes[no].center[1];
      s[2] = Nodes[no].center[2];
      vs[0] = 0;
      vs[1] = 0;
      vs[2] = 0;
    }



  Nodes[no].u.d.s[0] = s[0];
  Nodes[no].u.d.s[1] = s[1];
  Nodes[no].u.d.s[2] = s[2];
  Extnodes[no].vs[0] = vs[0];
  Extnodes[no].vs[1] = vs[1];
  Extnodes[no].vs[2] = vs[2];
  Nodes[no].u.d.mass = mass;

  Extnodes[no].hmax = hmax;
  Extnodes[no].vmax = vmax;

  Extnodes[no].Flag = GlobFlag;


  if(count_particles > 1)
    multiple_flag = (1 << BITFLAG_MULTIPLEPARTICLES);
  else
    multiple_flag = 0;

  Nodes[no].u.d.bitflags &= (~BITFLAG_MASK);	/* this clears the bits */

  Nodes[no].u.d.bitflags |= multiple_flag;

#ifdef UNEQUALSOFTENINGS
  Nodes[no].u.d.bitflags |= diffsoftflag + (maxsofttype << BITFLAG_MAX_SOFTENING_TYPE);
#endif

}



/*! This function flags nodes in the top-level tree that are dependent on
 *  local particle data.
 */
void force_flag_localnodes(void)
{
  int no, i, m;

  /* mark all top-level nodes */

  for(i = 0; i < NTopleaves; i++)
    {
      no = DomainNodeIndex[i];

      while(no >= 0)
	{
	  if(Nodes[no].u.d.bitflags & (1 << BITFLAG_TOPLEVEL))
	    break;

	  Nodes[no].u.d.bitflags |= (1 << BITFLAG_TOPLEVEL);

	  no = Nodes[no].u.d.father;
	}

      /* mark also internal top level nodes */

      no = DomainNodeIndex[i];
      no = Nodes[no].u.d.father;

      while(no >= 0)
	{
	  if(Nodes[no].u.d.bitflags & (1 << BITFLAG_INTERNAL_TOPLEVEL))
	    break;

	  Nodes[no].u.d.bitflags |= (1 << BITFLAG_INTERNAL_TOPLEVEL);

	  no = Nodes[no].u.d.father;
	}
    }

  /* mark top-level nodes that contain local particles */

  for(m = 0; m < MULTIPLEDOMAINS; m++)
    for(i = DomainStartList[ThisTask * MULTIPLEDOMAINS + m];
	i <= DomainEndList[ThisTask * MULTIPLEDOMAINS + m]; i++)
      {
	no = DomainNodeIndex[i];

	if(DomainTask[i] != ThisTask)
	  endrun(131231231);

	while(no >= 0)
	  {
	    if(Nodes[no].u.d.bitflags & (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS))
	      break;

	    Nodes[no].u.d.bitflags |= (1 << BITFLAG_DEPENDS_ON_LOCAL_MASS);

	    no = Nodes[no].u.d.father;
	  }
      }
}


/*! When a new additional star particle is created, we can put it into the
 *  tree at the position of the spawning gas particle. This is possible
 *  because the Nextnode[] array essentially describes the full tree walk as a
 *  link list. Multipole moments of tree nodes need not be changed.
 */
/* GM: not needed now
void force_add_star_to_tree(int igas, int istar)
{
  int no;

  no = Nextnode[igas];
  Nextnode[igas] = istar;
  Nextnode[istar] = no;
  Father[istar] = Father[igas];
}
*/



/*! This function allocates the memory used for storage of the tree and of
 *  auxiliary arrays needed for tree-walk and link-lists.  Usually,
 *  maxnodes approximately equal to 0.7*maxpart is sufficient to store the
 *  tree for up to maxpart particles.
 */
void force_treeallocate(int maxnodes, int maxpart)
{
  int i;
  size_t bytes;
  double allbytes = 0, allbytes_topleaves = 0;
  double u;

  tree_allocated_flag = 1;
  DomainNodeIndex = (int *) mymalloc("DomainNodeIndex", bytes = NTopleaves * sizeof(int));
  allbytes_topleaves += bytes;
  MaxNodes = maxnodes;
  if(!(Nodes_base = (struct NODE *) mymalloc("Nodes_base", bytes = (MaxNodes + 1) * sizeof(struct NODE))))
    {
      printf("failed to allocate memory for %d tree-nodes (%g MB).\n", MaxNodes, bytes / (1024.0 * 1024.0));
      endrun(3);
    }
  allbytes += bytes;
  if(!
     (Extnodes_base =
      (struct extNODE *) mymalloc("Extnodes_base", bytes = (MaxNodes + 1) * sizeof(struct extNODE))))
    {
      printf("failed to allocate memory for %d tree-extnodes (%g MB).\n",
	     MaxNodes, bytes / (1024.0 * 1024.0));
      endrun(3);
    }
  allbytes += bytes;
  Nodes = Nodes_base - All.MaxPart;
  Extnodes = Extnodes_base - All.MaxPart;
  if(!(Nextnode = (int *) mymalloc("Nextnode", bytes = (maxpart + NTopnodes) * sizeof(int))))
    {
      printf("Failed to allocate %d spaces for 'Nextnode' array (%g MB)\n",
	     maxpart + NTopnodes, bytes / (1024.0 * 1024.0));
      exit(0);
    }
  allbytes += bytes;
  if(!(Father = (int *) mymalloc("Father", bytes = (maxpart) * sizeof(int))))
    {
      printf("Failed to allocate %d spaces for 'Father' array (%g MB)\n", maxpart, bytes / (1024.0 * 1024.0));
      exit(0);
    }
  allbytes += bytes;
  if(first_flag == 0)
    {
      first_flag = 1;
      if(ThisTask == 0)
	printf
	  ("\nAllocated %g MByte for BH-tree, and %g Mbyte for top-leaves.  (presently allocted %g MB)\n\n",
	   allbytes / (1024.0 * 1024.0), allbytes_topleaves / (1024.0 * 1024.0),
	   AllocatedBytes / (1024.0 * 1024.0));
      for(i = 0; i < NTAB; i++)
	{
	  u = 3.0 / NTAB * (i + 0.5);
	  shortrange_table[i] = erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u);
	  shortrange_table_potential[i] = erfc(u);
#ifdef DISTORTIONTENSORPS
	  shortrange_table_tidal[i] = 4.0 * u * u * u / sqrt(M_PI) * exp(-u * u);
#endif
	}
    }
}


/*! This function frees the memory allocated for the tree, i.e. it frees
 *  the space allocated by the function force_treeallocate().
 */
void force_treefree(void)
{
  if(tree_allocated_flag)
    {
      myfree(Father);
      myfree(Nextnode);
      myfree(Extnodes_base);
      myfree(Nodes_base);
      myfree(DomainNodeIndex);
      tree_allocated_flag = 0;
    }
}





/*! This function dumps some of the basic particle data to a file. In case
 *  the tree construction fails, it is called just before the run
 *  terminates with an error message. Examination of the generated file may
 *  then give clues to what caused the problem.
 */
/* GM: I leave this one here!! */
void dump_particles(void)
{
  FILE *fd;
  char buffer[200];
  int i;

  sprintf(buffer, "particles%d.dat", ThisTask);
  fd = fopen(buffer, "w");
  my_fwrite(&NumPart, 1, sizeof(int), fd);
  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Pos[0], 3, sizeof(MyFloat), fd);
  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Vel[0], 3, sizeof(MyFloat), fd);
  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].ID, 1, sizeof(int), fd);
  fclose(fd);
}


/* GM: note, this was in gravtree.c */
/*! This function is used as a comparison kernel in a sort routine. It is
 *  used to group particles in the communication buffer that are going to
 *  be sent to the same CPU.
    */
int data_index_compare(const void *a, const void *b)
{
  if(((struct data_index *) a)->Task < (((struct data_index *) b)->Task))
    return -1;

  if(((struct data_index *) a)->Task > (((struct data_index *) b)->Task))
    return +1;

  if(((struct data_index *) a)->Index < (((struct data_index *) b)->Index))
    return -1;

  if(((struct data_index *) a)->Index > (((struct data_index *) b)->Index))
    return +1;

  if(((struct data_index *) a)->IndexGet < (((struct data_index *) b)->IndexGet))
    return -1;

  if(((struct data_index *) a)->IndexGet > (((struct data_index *) b)->IndexGet))
    return +1;

  return 0;
}



/*! This routine computes the gravitational force for a given local
 *  particle, or for a particle in the communication buffer. Depending on
 *  the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
*/
       /* GM: here we had these functions, now, not needed
int force_treeevaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex)
       */

#ifdef PMGRID
/*! In the TreePM algorithm, the tree is walked only locally around the
 *  target coordinate.  Tree nodes that fall outside a box of half
 *  side-length Rcut= RCUT*ASMTH*MeshSize can be discarded. The short-range
 *  potential is modified by a complementary error function, multiplied
 *  with the Newtonian form. The resulting short-range suppression compared
 *  to the Newtonian force is tabulated, because looking up from this table
 *  is faster than recomputing the corresponding factor, despite the
 *  memory-access panelty (which reduces cache performance) incurred by the
 *  table.
 */
       /* GM: now not needed
int force_treeevaluate_shortrange(int target, int mode, int *exportflag, int *exportnodecount,
				  int *exportindex)
       */
#endif






#ifdef PERIODIC
/*! This function computes the Ewald correction, and is needed if periodic
 *  boundary conditions together with a pure tree algorithm are used. Note
 *  that the ordinary tree walk does not carry out this correction directly
 *  as it was done in Gadget-1.1. Instead, the tree is walked a second
 *  time. This is actually faster because the "Ewald-Treewalk" can use a
 *  different opening criterion than the normal tree walk. In particular,
 *  the Ewald correction is negligible for particles that are very close,
 *  but it is large for particles that are far away (this is quite
 *  different for the normal direct force). So we can here use a different
 *  opening criterion. Sufficient accuracy is usually obtained if the node
 *  length has dropped to a certain fraction ~< 0.25 of the
 *  BoxLength. However, we may only short-cut the interaction list of the
 *  normal full Ewald tree walk if we are sure that the whole node and all
 *  daughter nodes "lie on the same side" of the periodic boundary,
 *  i.e. that the real tree walk would not find a daughter node or particle
 *  that was mapped to a different nearest neighbour position when the tree
 *  walk would be further refined.
 */
       /* GM: now not needed
int force_treeevaluate_ewald_correction(int target, int mode, int *exportflag, int *exportnodecount,
					int *exportindex)
       */
#endif






/*! This routine computes the gravitational potential by walking the
 *  tree. The same opening criteria is used as for the gravitational force
 *  walk.
 */
       /* GM: now not needed
int force_treeevaluate_potential(int target, int mode, int *nexport, int *nsend_local)
       */


/* GM: here we had a SubFind part, now not needed
#ifdef SUBFIND
int subfind_force_treeevaluate_potential(int target, int mode, int *nexport, int *nsend_local)
       
#endif
*/


#ifdef PMGRID
/*! This function computes the short-range potential when the TreePM
 *  algorithm is used. This potential is the Newtonian potential, modified
 *  by a complementary error function.
 */
       /* GM: now not needed
int force_treeevaluate_potential_shortrange(int target, int mode, int *nexport, int *nsend_local)
       */
#endif




/*! This function does the force computation with direct summation for the
 *  specified particle in the communication buffer. This can be useful for
 *  debugging purposes, in particular for explicit checks of the force
 *  accuracy.
 */
/* GM: now not needed
#ifdef FORCETEST
int force_treeevaluate_direct(int target, int mode)
#endif
*/

/* GM: now not needed
#ifdef PHIDOT

int phidot_treeevaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex)

int phidot_treeevaluate_ewald_correction(int target, int mode, int *exportflag, int *exportnodecount,
					 int *exportindex)

double phidot_ewald_corr(double dx, double dy, double dz, double vx, double vy, double vz)
#endif // end PHIDOT
*/




#ifdef PERIODIC

/*! This function initializes tables with the correction force and the
 *  correction potential due to the periodic images of a point mass located
 *  at the origin. These corrections are obtained by Ewald summation. (See
 *  Hernquist, Bouchet, Suto, ApJS, 1991, 75, 231) The correction fields
 *  are used to obtain the full periodic force if periodic boundaries
 *  combined with the pure tree algorithm are used. For the TreePM
 *  algorithm, the Ewald correction is not used.
 *
 *  The correction fields are stored on disk once they are computed. If a
 *  corresponding file is found, they are loaded from disk to speed up the
 *  initialization.  The Ewald summation is done in parallel, i.e. the
 *  processors share the work to compute the tables if needed.
 */
/* GM: now not needed
void ewald_init(void)
*/


/*! This function looks up the correction force due to the infinite number
 *  of periodic particle/node images. We here use trilinear interpolation
 *  to get it from the precomputed tables, which contain one octant
 *  around the target particle at the origin. The other octants are
 *  obtained from it by exploiting the symmetry properties.
 */
/* GM: now not needed
#ifdef FORCETEST
void ewald_corr(double dx, double dy, double dz, double *fper)
#endif
*/


/*! This function looks up the correction potential due to the infinite
 *  number of periodic particle/node images. We here use tri-linear
 *  interpolation to get it from the precomputed table, which contains
 *  one octant around the target particle at the origin. The other
 *  octants are obtained from it by exploiting symmetry properties.
 */
/* GM: now not needed
double ewald_pot_corr(double dx, double dy, double dz)
*/



/*! This function computes the potential correction term by means of Ewald
 *  summation.
 */
/* GM: now not needed
double ewald_psi(double x[3])
*/


/*! This function computes the force correction term (difference between full
 *  force of infinite lattice and nearest image) by Ewald summation.
 */
/* GM: now not needed
void ewald_force(int iii, int jjj, int kkk, double x[3], double force[3])
*/
#endif

