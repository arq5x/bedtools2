/*  thread_pool_internal.h -- Internal API for the thread pool.

    Copyright (c) 2013-2016 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

/*
 * This file implements a thread pool for multi-threading applications.
 * It consists of two distinct interfaces: thread pools an thread job queues.
 *
 * The pool of threads is given a function pointer and void* data to pass in.
 * This means the pool can run jobs of multiple types, albeit first come
 * first served with no job scheduling except to pick tasks from
 * queues that have room to store the result.
 *
 * Upon completion, the return value from the function pointer is
 * added to back to the queue if the result is required.  We may have
 * multiple queues in use for the one pool.
 *
 * To see example usage, please look at the #ifdef TEST_MAIN code in
 * thread_pool.c.
 */

#ifndef THREAD_POOL_INTERNAL_H
#define THREAD_POOL_INTERNAL_H

#include <pthread.h>
#include <stdint.h>
#include "htslib/thread_pool.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * An input job, before execution.
 */
typedef struct hts_tpool_job {
    void *(*func)(void *arg);
    void *arg;
    struct hts_tpool_job *next;

    struct hts_tpool *p;
    struct hts_tpool_process *q;
    uint64_t serial;
} hts_tpool_job;

/*
 * An output, after job has executed.
 */
struct hts_tpool_result {
    struct hts_tpool_result *next;
    uint64_t serial; // sequential number for ordering
    void *data;      // result itself
};

/*
 * A per-thread worker struct.
 */
typedef struct {
    struct hts_tpool *p;
    int idx;
    pthread_t tid;
    pthread_cond_t  pending_c; // when waiting for a job
} hts_tpool_worker;

/*
 * An IO queue consists of a queue of jobs to execute
 * (the "input" side) and a queue of job results post-
 * execution (the "output" side).
 *
 * We have size limits to prevent either queue from
 * growing too large and serial numbers to ensure
 * sequential consumption of the output.
 *
 * The thread pool may have many hetergeneous tasks, each
 * using its own io_queue mixed into the same thread pool.
 */
struct hts_tpool_process {
    struct hts_tpool *p;             // thread pool
    hts_tpool_job    *input_head;    // input list
    hts_tpool_job    *input_tail;
    hts_tpool_result *output_head;   // output list
    hts_tpool_result *output_tail;
    int qsize;                       // max size of i/o queues
    uint64_t next_serial;            // next serial for output
    uint64_t curr_serial;            // current serial (next input)

    int n_input;                     // no. items in input queue; was njobs
    int n_output;                    // no. items in output queue
    int n_processing;                // no. items being processed (executing)

    int shutdown;                    // true if pool is being destroyed
    int in_only;                     // if true, don't queue result up.
    int wake_dispatch;               // unblocks waiting dispatchers

    int ref_count;                   // used to track safe destruction

    pthread_cond_t output_avail_c;   // Signalled on each new output
    pthread_cond_t input_not_full_c; // Input queue is no longer full
    pthread_cond_t input_empty_c;    // Input queue has become empty
    pthread_cond_t none_processing_c;// n_processing has hit zero

    struct hts_tpool_process *next, *prev;// to form circular linked list.
};

/*
 * The single pool structure itself.
 *
 * This knows nothing about the nature of the jobs or where their
 * output is going, but it maintains a list of queues associated with
 * this pool from which the jobs are taken.
 */
struct hts_tpool {
    int nwaiting; // how many workers waiting for new jobs
    int njobs;    // how many total jobs are waiting in all queues
    int shutdown; // true if pool is being destroyed

    // I/O queues to check for jobs in and to put results.
    // Forms a circular linked list.  (q_head may be amended
    // to point to the most recently updated.)
    hts_tpool_process *q_head;

    // threads
    int tsize;    // maximum number of jobs
    hts_tpool_worker *t;
    // array of worker IDs free
    int *t_stack, t_stack_top;

    // A single mutex used when updating this and any associated structure.
    pthread_mutex_t pool_m;

    // Tracking of average number of running jobs.
    // This can be used to dampen any hysteresis caused by bursty
    // input availability.
    int n_count, n_running;

    // Debugging to check wait time.
    // FIXME: should we just delete these and cull the associated code?
    long long total_time, wait_time;
};

#ifdef __cplusplus
}
#endif

#endif
