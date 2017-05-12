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
#ifndef ALLVARS_H
#include "allvars.h"
#endif
#include "forcetree.h"
#include "domain.h"

/* GM: PLEASE put prototipes under the name of the file that
  contains the function, and IN THE SAME ORDER in which they
  can be found in the file! */


/* general */
static inline double DMAX(double a, double b)
{
  return (a > b) ? a : b;
}

static inline double DMIN(double a, double b)
{
  return (a < b) ? a : b;
}

static inline int IMAX(int a, int b)
{
  return (a > b) ? a : b;
}

static inline int IMIN(int a, int b)
{
  return (a < b) ? a : b;
}


/* main.c */
/* nothing to declare */


/* domain.c */
/* protos in domain.h */


/* begrun.c */
void begrun(void);
void set_units(void);
void open_outputfiles(void);
void read_parameter_file(char *);


/* run.c */
void run(void);
void output_log_messages(void);
void write_cpu_log(void);
void put_symbol(double, double, char);
void check_particles_info(const char *, const char *, int);


/* io.c */ 
void savepositions(int);
void fill_write_buffer(enum iofields, int *, int, int);
int get_bytes_per_blockelement(enum iofields, int);
int get_datatype_in_block(enum iofields);
int get_values_per_blockelement(enum iofields);
int get_particles_in_block(enum iofields, int *);
int blockpresent(enum iofields);
void get_Tab_IO_Label(enum iofields, char *);
void get_dataset_name(enum iofields, char *);
void write_file(char *, int, int);
size_t my_fwrite(void *, size_t, size_t, FILE *);
size_t my_fread(void *, size_t, size_t, FILE *);
void mpi_printf(const char *, ...);


/* init.c */
void init(void);
void test_id_uniqueness(void);
int compare_IDs(const void *, const void *);
void set_softenings(void);


/* forcetree.c */
/* protos in forcetree.h */


/* forcetree_update.c */
/* protos in forcetree.h */


/* read_ic.c */
void read_ic(char *);
void empty_read_buffer(enum iofields, int, int, int);
void read_file(char *, int, int);
void distribute_file(int, int, int, int, int *, int *, int *);
int find_files(char *);
#ifdef AUTO_SWAP_ENDIAN_READIC
void swap_Nbyte(char *, int, int);
void swap_header(void);
#endif
void find_block(char *, FILE *);


/* ngb.c */
int ngb_treefind_pairs(MyDoublePos searchcenter[3], MyFloat, int, int *, int, int *, int *);
int ngb_treefind_variable(MyDoublePos searchcenter[3], MyFloat, int, int *, int, int *, int *);
void ngb_init(void);
void ngb_treebuild(void);


/* system.c */
void enable_core_dumps_and_fpu_exceptions(void);
void catch_abort(int);
void catch_fatal(int);
void terminate_processes(void);
void write_pid_file(void);
double get_random_number(MyIDType);
void set_random_numbers(void);
double second(void);
double measure_time(void);
double report_time(void);
double timediff(double, double);
void minimum_large_ints(int, long long *, long long *);
void sumup_large_ints(int, int *, long long *);
void sumup_longs(int, long long *, long long *);
size_t sizemax(size_t, size_t);
void report_VmRSS(void);
long long report_comittable_memory(long long *, long long *, long long *, long long *);
void mpi_report_comittable_memory(void);
void report_task_list(void);
void write_ps_files(void);


/* allocate.c */
void allocate_memory(void);


/* endrun.c */
void endrun(int ierr);


/* parallel_sort.c */
void parallel_sort(void *, size_t, size_t, int (*compar) (const void *, const void *));
void parallel_sort_comm(void *, size_t, size_t, int (*compar) (const void *, const void *), MPI_Comm);


/* peano.c */
void peano_hilbert_order(void);
int peano_compare_key(const void *, const void *);
void reorder_gas(void);
void reorder_particles(void);
peanokey peano_hilbert_key(int, int, int, int);
peanokey morton_key(int, int, int, int);
peanokey peano_and_morton_key(int, int, int, int, peanokey *);
peanokey peano_hilbert_key_old(int, int, int, int);
peanokey peano_and_morton_key_old(int, int, int, int, peanokey *);
peanokey morton_key_old(int, int, int, int);
void mysort_peano(void *, size_t, size_t, int (*cmp) (const void *, const void *));


/* mpi_util.c */
int mpi_calculate_offsets(int *, int *, int *, int *, int);
int intpointer_compare(const void *, const void *);
void sort_based_on_field(void *, int, int, int, void **);
void mpi_distribute_items_to_tasks(void *, int, int *, int *, int);


/* mymalloc.c */
void mymalloc_init(void);
void report_detailed_memory_usage_of_largest_task(size_t *, const char *, const char *, const char *, int);
void dump_memory_table(void);
void *mymalloc_fullinfo(const char *, size_t, const char *, const char *, int);
void *mymalloc_movable_fullinfo(void *, const char *, size_t, const char *, const char *, int);
void myfree_fullinfo(void *, const char *, const char *, int);
void myfree_movable_fullinfo(void *, const char *, const char *, int);
void *myrealloc_fullinfo(void *, size_t, const char *, const char *, int);
void *myrealloc_movable_fullinfo(void *, size_t, const char *, const char *, int);


/* compile_time_info.c */
void output_compile_time_options(void);


/* look_around.c*/
void look_around(int, double);
int la_evaluate(int target, double, int, int *, int *, int *, int *);
void *la_evaluate_primary(int, double);
void *la_evaluate_secondary(void);
