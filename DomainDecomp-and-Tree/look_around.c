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
#include <gsl/gsl_math.h>


#include "allvars.h"
#include "proto.h"


#ifndef DEBUG
#define NDEBUG
#endif
#include <assert.h>

#define NO_OUTPUT


extern int NextParticle;
extern int Nexport, Nimport;
extern int BufferFullFlag;
extern int NextJ;
extern int TimerFlag;


/*! \file look_around.c
*  \brief 
*
*  This file looks for neighbout of one single particle Iparticle
*  within a cube of size cSize
*  writes the positions of neighbours on a file named part%0d_task%03d.dat
*/

int LocalNgb;

struct la_in
{
  MyDoublePos Pos[3];
  MyIDType Id;
  int NodeList[NODELISTLENGTH];
  int task;
  MyFloat Hsml;
}
 *laDataIn, *laDataGet;


struct la_out
{
  int ExternNgb;
}
 *laDataResult, *laDataOut;


/***************************************/
/*      look around exchange functions */
/***************************************/
void look_around(int Iparticle, double cSize)
{
  int j, k, ngrp, ndone, ndone_flag;
  int sendTask, recvTask, place;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait, tstart, tend, t0, t1;
  int ExternNgb, TotNgb;
  char bfr[100];
  FILE *ff;
  long long n_exported = 0;

  /* allocate buffers to arrange communication */
  ExternNgb=TotNgb=0;

  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct la_in) +
					     sizeof(struct la_out) +
					     sizemax(sizeof(struct la_in),
						     sizeof(struct la_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  t0 = second();


  //  do
  //    {

      BufferFullFlag = 0;
      Nexport = 0;

      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      tstart = second();

      la_evaluate_primary(Iparticle, cSize);	/* do local particle and prepare export list */

      tend = second();
      timecomp1 += timediff(tstart, tend);

      if(BufferFullFlag)
	{
	  int new_export = 0;

	  for(j = 0, k = 0; j < Nexport; j++)
	    if(ProcessedFlag[DataIndexTable[j].Index] != 2)
	      {
		if(k < j + 1)
		  k = j + 1;

		for(; k < Nexport; k++)
		  if(ProcessedFlag[DataIndexTable[k].Index] == 2)
		    {
		      int old_index = DataIndexTable[j].Index;

		      DataIndexTable[j] = DataIndexTable[k];
		      DataNodeList[j] = DataNodeList[k];
		      DataIndexTable[j].IndexGet = j;
		      new_export++;

		      DataIndexTable[k].Index = old_index;
		      k++;
		      break;
		    }
	      }
	    else
	      new_export++;

	  Nexport = new_export;
	}

      n_exported += Nexport;

      for(j = 0; j < NTask; j++)
	Send_count[j] = 0;
      for(j = 0; j < Nexport; j++)
	Send_count[DataIndexTable[j].Task]++;

      qsort(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);

      tstart = second();
      
      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MYMPI_COMM_WORLD);
      /* apparentemente l'ultima particella non riceve giusto */

      tend = second();
      timewait1 += timediff(tstart, tend);

      for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      laDataGet = (struct la_in *) mymalloc("LookArountDataGet", Nimport * sizeof(struct la_in));
      laDataIn = (struct la_in *) mymalloc("LookArountDataIn", Nexport * sizeof(struct la_in));

      /* prepare particle data for export */
      for(j = 0; j < Nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  for(k = 0; k < 3; k++)
	    {
	      laDataIn[j].Pos[k] = P[place].Pos[k];
	    }
	  laDataIn[j].Hsml = cSize;
	  laDataIn[j].Id = P[place].ID;
	  laDataIn[j].task = ThisTask;

	  memcpy(laDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));

	}


      /* exchange particle data */
      tstart = second();
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&laDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct la_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A,
			       &laDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct la_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();
      timecommsumm1 += timediff(tstart, tend);


      myfree(laDataIn);
      laDataResult =
	(struct la_out *) mymalloc("LookAroundDataResult", Nimport * sizeof(struct la_out));
      laDataOut =
	(struct la_out *) mymalloc("LookAroundOut", Nexport * sizeof(struct la_out));

      /* now do the particle that was sent to us, if any */

      tstart = second();

      NextJ = 0;
      la_evaluate_secondary();


      tend = second();
      timecomp2 += timediff(tstart, tend);


      tstart = second();
      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MYMPI_COMM_WORLD);
      tend = second();
      timewait2 += timediff(tstart, tend);

      /* get the result */
      tstart = second();
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&laDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct la_out),
			       MPI_BYTE, recvTask, TAG_HYDRO_B,
			       &laDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct la_out),
			       MPI_BYTE, recvTask, TAG_HYDRO_B, MYMPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();
      timecommsumm2 += timediff(tstart, tend);


      /* add the result to the local particles */
      for(j = 0; j < Nexport; j++)
        {
          place = DataIndexTable[j].Index;
          ExternNgb += laDataOut[j].ExternNgb;
        }


      myfree(laDataOut);
      myfree(laDataResult);
      myfree(laDataGet);
      //    }
//  while(ndone < NTask);

  TotNgb = LocalNgb + ExternNgb;
  sprintf(bfr,"part%06u_task%02d.dat",P[Iparticle].ID, ThisTask);
#ifndef NO_OUTPUT
  ff=fopen(bfr,"a");
  fprintf(ff,"#  local %d extern %d tot %d\n",LocalNgb, ExternNgb, TotNgb); fflush(ff);
  fclose(ff);
#endif

  myfree(DataNodeList);
  myfree(DataIndexTable);

  myfree(Ngblist);



  /* collect some timing information */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1); /* total */

  timecomp = timecomp1 + timecomp2; /* computation */
  timewait = timewait1 + timewait2; /* waiting other procs */
  timecomm = timecommsumm1 + timecommsumm2; /* communications */


  CPU_Step[CPU_COMP] += timecomp; 
  CPU_Step[CPU_COMM] += timecomm;
  CPU_Step[CPU_TIMEWAIT] += timewait;
}



/*! 
*  
*/
int la_evaluate(int target, double cSize, int mode, 
		int *nexport, int *exportnodecount, int *exportindex,
		int *ngblist)
{
  int startnode, numngb, listindex = 0;
  int j, n, dummy;
  MyDoublePos *pos;
  MyFloat h;
  char bfr[100];
  FILE *f;




  if(mode == 0)
    {
      pos = P[target].Pos;
      h = cSize;
      sprintf(bfr,"part%06u_task%02d.dat",P[target].ID, ThisTask);
    } 
  else
    {
      pos = laDataGet[target].Pos;
      h = laDataGet[target].Hsml;
      sprintf(bfr,"part%06u_task%02d.dat",laDataGet[target].Id, laDataGet[target].task);
    }
#ifndef NO_OUTPUT
      f=fopen(bfr,"a");
#endif


  /* Now start the actual neighbour search for this particle */

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = laDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{

	  numngb =
	    ngb_treefind_variable(pos, h, target, &startnode, mode, nexport, &dummy);


	  if(numngb < 0)
	    {
	      if(mode>0)
		laDataResult[target].ExternNgb = 0;
	      return -1;
	    }

	  for(n = 0; n < numngb; n++)
	    {
	      j = ngblist[n];
#ifndef NO_OUTPUT
	      fprintf(f,"%f %f %f   %d %u  %d\n",P[j].Pos[0],P[j].Pos[1],P[j].Pos[2],
		      P[j].Type, P[j].ID, ThisTask); 
#endif
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = laDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }


  if(mode==0)
    {
      LocalNgb = numngb;
    } 
  else 
    {
      laDataResult[target].ExternNgb = numngb;
    }

#ifndef NO_OUTPUT
  fclose(f);
#endif
  return 0;
}




/* This function processes ONE single particle i */
void *la_evaluate_primary(int i, double cS)
{

  int j;

  int *exportflag, *exportnodecount, *exportindex, *ngblist, *nexport;


  ngblist = Ngblist;
  exportflag = Exportflag;
  exportnodecount = Exportnodecount;
  exportindex = Exportindex;
  nexport = &Nexport;

  for(j = 0; j < NTask; j++)
    exportflag[j] = -1;


  la_evaluate(i, cS, 0, nexport, exportnodecount, exportindex, ngblist);


  return NULL;

}



void *la_evaluate_secondary(void)
{

  int j, dummy, *ngblist;

  ngblist = Ngblist;

  while(1)
    {
      j = NextJ;
      NextJ++;

      if(j >= Nimport)
	break;

      la_evaluate(j, dummy, 1, &dummy, &dummy, &dummy, ngblist);

    }

  return NULL;

}

