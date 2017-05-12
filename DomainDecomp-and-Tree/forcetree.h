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
#ifndef FORCETREE_H
#define FORCETREE_H

#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif


#define BITFLAG_TOPLEVEL                   0
#define BITFLAG_DEPENDS_ON_LOCAL_MASS      1
#define BITFLAG_MAX_SOFTENING_TYPE         2	/* bits 2-4 */
#define BITFLAG_MIXED_SOFTENINGS_IN_NODE   5
#define BITFLAG_INTERNAL_TOPLEVEL          6
#define BITFLAG_MULTIPLEPARTICLES          7
#define BITFLAG_NODEHASBEENKICKED          8
#define BITFLAG_INSIDE_LINKINGLENGTH       9

#define BITFLAG_MASK  ((1 << BITFLAG_MULTIPLEPARTICLES) + (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE) + (7 << BITFLAG_MAX_SOFTENING_TYPE))
#define maskout_different_softening_flag(x) (x & (1 << BITFLAG_MIXED_SOFTENINGS_IN_NODE))
#define extract_max_softening_type(x) ((x >> BITFLAG_MAX_SOFTENING_TYPE) & 7)


/* forcetree.c */
int    force_treebuild(int, struct unbind_data *);
int    force_treebuild_single(int, struct unbind_data *);
void force_create_empty_nodes(int, int, int, int, int, int, int *, int *);
void force_insert_pseudo_particles(void);
void force_update_node_recursive(int, int, int, float);
void force_exchange_pseudodata(void);
void force_treeupdate_pseudos(int);
void force_flag_localnodes(void);
void   force_treeallocate(int, int);  
void   force_treefree(void);
void   dump_particles(void);


/* these are commented out in forcetree.c 
int force_treeevaluate(int, int, int *, int *, int *);
int force_treeevaluate_shortrange(int, int, int *, int *, int *);
int force_treeevaluate_ewald_correction(int, int, int *, int *, int *);
int force_treeevaluate_potential(int target, int type, int *nexport, int *nsend_local);
int force_treeevaluate_potential_shortrange(int, int, int *, int *);
int    force_treeevaluate_direct(int, int);
*/


/* forcetree_update.c */
void force_update_tree(void);
void force_kick_node(int, MyFloat *);
void force_finish_kick_nodes(void);
void force_drift_node(int no, int time1);
void force_update_hmax(void);
int data_index_compare(const void *, const void *);


/* ngb.c */
MyFloat  INLINE_FUNC ngb_periodic(MyFloat x);
MyFloat  INLINE_FUNC ngb_periodic_longbox(MyFloat x);
int ngb_treefind_pairs(MyDoublePos searchcenter[3], MyFloat, int, int *, int, int *, int *);
int ngb_treefind_variable(MyDoublePos searchcenter[3], MyFloat, int, int *, int, int *, int *);
void ngb_init(void);
void   ngb_treebuild(void);


/* these functions are in modules not currently included 
void *gravity_primary_loop(void *p);
void *gravity_secondary_loop(void *p);
*/


#endif

