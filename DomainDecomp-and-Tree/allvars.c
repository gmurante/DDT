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

/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define statements
 *     - add #include "allvars.h"
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#include "allvars.h"




#ifdef PERIODIC
MyDouble boxSize, boxHalf, inverse_boxSize;
#ifdef LONG_X
MyDouble boxSize_X, boxHalf_X, inverse_boxSize_X;
#endif
#ifdef LONG_Y
MyDouble boxSize_Y, boxHalf_Y, inverse_boxSize_Y;
#endif
#ifdef LONG_Z
MyDouble boxSize_Z, boxHalf_Z, inverse_boxSize_Z;
#endif
#endif


/*********************************************************/
/*  Global variables                                     */
/*********************************************************/



int FirstActiveParticle;
int *NextActiveParticle;
unsigned char *ProcessedFlag;

int TimeBinCount[TIMEBINS];
int TimeBinCountSph[TIMEBINS];
int TimeBinActive[TIMEBINS];

int FirstInTimeBin[TIMEBINS];
int LastInTimeBin[TIMEBINS];
int *NextInTimeBin;
int *PrevInTimeBin;

int ThisTask;		/*!< the number of the local processor  */
int NTask;		/*!< number of processors */
int PTask;		/*!< note: NTask = 2^PTask */
MPI_Comm MYMPI_COMM_WORLD;



double CPUThisRun;	/*!< Sums CPU time of current process */


int NumForceUpdate;	/*!< number of active particles on local processor in current timestep  */
long long GlobNumForceUpdate;

int NumSphUpdate;	/*!< number of active SPH particles on local processor in current timestep  */

int MaxTopNodes;	        /*!< Maximum number of nodes in the top-level tree used for domain decomposition */

int RestartFlag;		/*!< taken from command line used to start code. 0 is normal start-up from
				   initial conditions, 1 is resuming a run from a set of restart files, while 2
				   marks a restart from a snapshot file. */
int RestartSnapNum;
int SelRnd;

int TakeLevel;

int *Exportflag;	        /*!< Buffer used for flagging whether a particle needs to be exported to another process */
int *Exportnodecount;
int *Exportindex;

int *Send_offset, *Send_count, *Recv_count, *Recv_offset;


size_t AllocatedBytes;
size_t HighMarkBytes;
size_t FreeBytes;

double CPU_Step[CPU_PARTS];
char CPU_Symbol[CPU_PARTS];
char CPU_SymbolImbalance[CPU_PARTS];
char CPU_String[CPU_STRING_LEN + 1];

double WallclockTime;    /*!< This holds the last wallclock time measurement for timings measurements */

int Flag_FullStep;	/*!< Flag used to signal that the current step involves all particles */

size_t HighMark_run,  HighMark_domain, HighMark_gravtree, HighMark_pmperiodic,
  HighMark_pmnonperiodic,  HighMark_sphdensity, HighMark_sphhydro, HighMark_addSPH;



int TreeReconstructFlag;
int GlobFlag;
char DumpFlag;


int NumPart;		/*!< number of particles on the LOCAL processor */
int N_gas;		/*!< number of gas particles on the LOCAL processor  */


long long Ntype[6];	/*!< total number of particles of each type */
int NtypeLocal[6];	/*!< local number of particles of each type */

gsl_rng *random_generator;	/*!< the random number generator used */


struct topnode_data *TopNodes;

double TimeOfLastTreeConstruction;	/*!< holds what it says */

int *Ngblist;		/*!< Buffer to hold indices of neighbours retrieved by the neighbour search */


int NTopnodes, NTopleaves;

double *R2ngblist;

double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
int *DomainStartList, *DomainEndList;

double *DomainWork;
int *DomainCount;
int *DomainCountSph;
int *DomainTask;
int *DomainNodeIndex;
int *DomainList, DomainNumChanged;

peanokey *Key, *KeySorted;

double RndTable[RNDTABLE];




/* variables for input/output , usually only used on process 0 */


char ParameterFile[100];	/*!< file name of parameterfile used for starting the simulation */

FILE *FdInfo,		/*!< file handle for info.txt log-file. */
 *FdEnergy,			/*!< file handle for energy.txt log-file. */
 *FdTimings,			/*!< file handle for timings.txt log-file. */
 *FdBalance,			/*!< file handle for balance.txt log-file. */
 *FdCPU,			/*!< file handle for cpu.txt log-file. */
 *FdTimebin;


/*! table for the cosmological drift factors */
double DriftTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for gravitational forces */
double GravKickTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for hydrodynmical forces */
double HydroKickTable[DRIFT_TABLE_LENGTH];



void *CommBuffer;	/*!< points to communication buffer, which is used at a few places */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
/* GM: I leave here a number of features that are not needed by DDT, for future use. */
struct global_data_all_processes All;



/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data *P,				/*!< holds particle data on local processor */
 *DomainPartBuf;		/*!< buffer for particle data used in domain decomposition */




/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
struct sph_particle_data  *SphP,				/*!< holds SPH particle data on local processor */
  *DomainSphBuf;			/*!< buffer for SPH particle data in domain decomposition */


peanokey *DomainKeyBuf;



/* global state of system
*/
struct state_of_system SysState, SysStateAtStart, SysStateAtEnd;


/* Various structures for communication during the gravity computation.
 */

struct data_index  *DataIndexTable;		/*!< the particles to be exported are grouped
                                  by task-number. This table allows the
                                  results to be disentangled again and to be
                                  assigned to the correct particle */

struct data_nodelist *DataNodeList;

struct gravdata_in *GravDataIn,			/*!< holds particle data to be exported to other processors */
 *GravDataGet;			/*!< holds particle data imported from other processors */


struct gravdata_out *GravDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *GravDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


struct potdata_out *PotDataResult,		/*!< holds the partial results computed for imported particles. Note: We use GravDataResult = GravDataGet, such that the result replaces the imported data */
 *PotDataOut;			/*!< holds partial results received from other processors. This will overwrite the GravDataIn array */


struct info_block *InfoBlock;


/*! Header for the standard file format.
 */
struct io_header header;				/*!< holds header for snapshot files */



/*
 * Variables for Tree
 * ------------------
 */

int Nexport, Nimport;
int BufferFullFlag;
int NextParticle;
int NextGroup;
int NextJ;
int TimerFlag;

struct NODE *Nodes_base,			/*!< points to the actual memory allocted for the nodes */
 *Nodes;			/*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart]
				   gives the first allocated node */


struct extNODE *Extnodes, *Extnodes_base;


int MaxNodes;		/*!< maximum allowed number of internal nodes */
int Numnodestree;	/*!< number of (internal) nodes in each tree */


int *Nextnode;		/*!< gives next node in tree walk  (nodes array) */
int *Father;		/*!< gives parent node in tree (Prenodes array) */

int maxThreads;




                                                                                                                                                                                                                                                                                                 



