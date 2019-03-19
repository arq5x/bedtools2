/*  thread_pool.c -- A pool of generic worker threads

    Copyright (c) 2013-2017 Genome Research Ltd.

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

#ifndef TEST_MAIN
#include <config.h>
#endif

#include <stdlib.h>
#include <inttypes.h>
#include <signal.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include <stdarg.h>
#include <unistd.h>
#include <limits.h>

#include "thread_pool_internal.h"

//#define DEBUG

#ifdef DEBUG
static int worker_id(hts_tpool *p) {
    int i;
    pthread_t s = pthread_self();
    for (i = 0; i < p->tsize; i++) {
        if (pthread_equal(s, p->t[i].tid))
            return i;
    }
    return -1;
}

int DBG_OUT(FILE *fp, char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    return vfprintf(fp, fmt, args);
}
#else
#define DBG_OUT(...) do{}while(0)
#endif

/* ----------------------------------------------------------------------------
 * A process-queue to hold results from the thread pool.
 *
 * Each thread pool may have jobs of multiple types being queued up and
 * interleaved, so we attach several job process-queues to a single pool.
 *
 * The jobs themselves are expected to push their results onto their
 * appropriate results queue.
 */

/*
 * Adds a result to the end of the process result queue.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
static int hts_tpool_add_result(hts_tpool_job *j, void *data) {
    hts_tpool_process *q = j->q;
    hts_tpool_result *r;

    pthread_mutex_lock(&q->p->pool_m);

    DBG_OUT(stderr, "%d: Adding result to queue %p, serial %"PRId64", %d of %d\n",
            worker_id(j->p), q, j->serial, q->n_output+1, q->qsize);

    if (--q->n_processing == 0)
        pthread_cond_signal(&q->none_processing_c);

    /* No results queue is fine if we don't want any results back */
    if (q->in_only) {
        pthread_mutex_unlock(&q->p->pool_m);
        return 0;
    }

    if (!(r = malloc(sizeof(*r))))
        return -1;

    r->next = NULL;
    r->data = data;
    r->serial = j->serial;

    q->n_output++;
    if (q->output_tail) {
        q->output_tail->next = r;
        q->output_tail = r;
    } else {
        q->output_head = q->output_tail = r;
    }

    assert(r->serial >= q->next_serial    // Or it will never be dequeued ...
           || q->next_serial == INT_MAX); // ... unless flush in progress.
    if (r->serial == q->next_serial) {
        DBG_OUT(stderr, "%d: Broadcasting result_avail (id %"PRId64")\n",
                worker_id(j->p), r->serial);
        pthread_cond_broadcast(&q->output_avail_c);
        DBG_OUT(stderr, "%d: Broadcast complete\n", worker_id(j->p));
    }

    pthread_mutex_unlock(&q->p->pool_m);

    return 0;
}

static void wake_next_worker(hts_tpool_process *q, int locked);

/* Core of hts_tpool_next_result() */
static hts_tpool_result *hts_tpool_next_result_locked(hts_tpool_process *q) {
    hts_tpool_result *r, *last;

    if (q->shutdown)
        return NULL;

    for (last = NULL, r = q->output_head; r; last = r, r = r->next) {
        if (r->serial == q->next_serial)
            break;
    }

    if (r) {
        // Remove r from out linked list
        if (q->output_head == r)
            q->output_head = r->next;
        else
            last->next = r->next;

        if (q->output_tail == r)
            q->output_tail = last;

        if (!q->output_head)
            q->output_tail = NULL;

        q->next_serial++;
        q->n_output--;

        if (q->qsize && q->n_output < q->qsize) {
            // Not technically input full, but can guarantee there is
            // room for the input to go somewhere so we still signal.
            // The waiting code will then check the condition again.
            if (q->n_input < q->qsize)
                pthread_cond_signal(&q->input_not_full_c);
            if (!q->shutdown)
                wake_next_worker(q, 1);
        }
    }

    return r;
}

/*
 * Pulls the next item off the process result queue.  The caller should free
 * it (and any internals as appropriate) after use.  This doesn't wait for a
 * result to be present.
 *
 * Results will be returned in strict order.
 *
 * Returns hts_tpool_result pointer if a result is ready.
 *         NULL if not.
 */
hts_tpool_result *hts_tpool_next_result(hts_tpool_process *q) {
    hts_tpool_result *r;

    DBG_OUT(stderr, "Requesting next result on queue %p\n", q);

    pthread_mutex_lock(&q->p->pool_m);
    r = hts_tpool_next_result_locked(q);
    pthread_mutex_unlock(&q->p->pool_m);

    DBG_OUT(stderr, "(q=%p) Found %p\n", q, r);

    return r;
}

/*
 * Pulls the next item off the process result queue.  The caller should free
 * it (and any internals as appropriate) after use.  This will wait for
 * a result to be present if none are currently available.
 *
 * Results will be returned in strict order.
 *
 * Returns hts_tpool_result pointer if a result is ready.
 *         NULL on error or during shutdown.
 */
hts_tpool_result *hts_tpool_next_result_wait(hts_tpool_process *q) {
    hts_tpool_result *r;

    pthread_mutex_lock(&q->p->pool_m);
    while (!(r = hts_tpool_next_result_locked(q))) {
        /* Possible race here now avoided via _locked() call, but incase... */
        struct timeval now;
        struct timespec timeout;

        gettimeofday(&now, NULL);
        timeout.tv_sec = now.tv_sec + 10;
        timeout.tv_nsec = now.tv_usec * 1000;

        q->ref_count++;
        if (q->shutdown) {
            int rc = --q->ref_count;
            pthread_mutex_unlock(&q->p->pool_m);
            if (rc == 0)
                hts_tpool_process_destroy(q);
            return NULL;
        }
        pthread_cond_timedwait(&q->output_avail_c, &q->p->pool_m, &timeout);

        q->ref_count--;
    }
    pthread_mutex_unlock(&q->p->pool_m);

    return r;
}

/*
 * Returns true if there are no items in the process results queue and
 * also none still pending.
 */
int hts_tpool_process_empty(hts_tpool_process *q) {
    int empty;

    pthread_mutex_lock(&q->p->pool_m);
    empty = q->n_input == 0 && q->n_processing == 0 && q->n_output == 0;
    pthread_mutex_unlock(&q->p->pool_m);

    return empty;
}

void hts_tpool_process_ref_incr(hts_tpool_process *q) {
    pthread_mutex_lock(&q->p->pool_m);
    q->ref_count++;
    pthread_mutex_unlock(&q->p->pool_m);
}

void hts_tpool_process_ref_decr(hts_tpool_process *q) {
    pthread_mutex_lock(&q->p->pool_m);
    if (--q->ref_count <= 0) {
        pthread_mutex_unlock(&q->p->pool_m);
        hts_tpool_process_destroy(q);
        return;
    }

    // maybe also call destroy here if needed?
    pthread_mutex_unlock(&q->p->pool_m);
}

/*
 * Returns the number of completed jobs in the process results queue.
 */
int hts_tpool_process_len(hts_tpool_process *q) {
    int len;

    pthread_mutex_lock(&q->p->pool_m);
    len = q->n_output;
    pthread_mutex_unlock(&q->p->pool_m);

    return len;
}

/*
 * Returns the number of completed jobs in the process results queue plus the
 * number running and queued up to run.
 */
int hts_tpool_process_sz(hts_tpool_process *q) {
    int len;

    pthread_mutex_lock(&q->p->pool_m);
    len = q->n_output + q->n_input + q->n_processing;
    pthread_mutex_unlock(&q->p->pool_m);

    return len;
}

/*
 * Shutdown a process.
 *
 * This sets the shutdown flag and wakes any threads waiting on process
 * condition variables.
 */
void hts_tpool_process_shutdown(hts_tpool_process *q) {
    pthread_mutex_lock(&q->p->pool_m);
    q->shutdown = 1;
    pthread_cond_broadcast(&q->output_avail_c);
    pthread_cond_broadcast(&q->input_not_full_c);
    pthread_cond_broadcast(&q->input_empty_c);
    pthread_cond_broadcast(&q->none_processing_c);
    pthread_mutex_unlock(&q->p->pool_m);
}

/*
 * Frees a result 'r' and if free_data is true also frees
 * the internal r->data result too.
 */
void hts_tpool_delete_result(hts_tpool_result *r, int free_data) {
    if (!r)
        return;

    if (free_data && r->data)
        free(r->data);

    free(r);
}

/*
 * Returns the data portion of a hts_tpool_result, corresponding
 * to the actual "result" itself.
 */
void *hts_tpool_result_data(hts_tpool_result *r) {
    return r->data;
}

/*
 * Initialises a thread process-queue.
 *
 * In_only, if true, indicates that the process generates does not need to
 * hold any output.  Otherwise an output queue is used to store the results
 * of processing each input job.
 *
 * Results hts_tpool_process pointer on success;
 *         NULL on failure
 */
hts_tpool_process *hts_tpool_process_init(hts_tpool *p, int qsize, int in_only) {
    hts_tpool_process *q = malloc(sizeof(*q));

    pthread_cond_init(&q->output_avail_c,   NULL);
    pthread_cond_init(&q->input_not_full_c, NULL);
    pthread_cond_init(&q->input_empty_c,    NULL);
    pthread_cond_init(&q->none_processing_c,NULL);

    q->p           = p;
    q->input_head  = NULL;
    q->input_tail  = NULL;
    q->output_head = NULL;
    q->output_tail = NULL;
    q->next_serial = 0;
    q->curr_serial = 0;
    q->n_input     = 0;
    q->n_output    = 0;
    q->n_processing= 0;
    q->qsize       = qsize;
    q->in_only     = in_only;
    q->shutdown    = 0;
    q->wake_dispatch = 0;
    q->ref_count   = 1;

    q->next        = NULL;
    q->prev        = NULL;

    hts_tpool_process_attach(p, q);

    return q;
}

/* Deallocates memory for a thread process-queue.
 * Must be called before the thread pool is destroyed.
 */
void hts_tpool_process_destroy(hts_tpool_process *q) {
    DBG_OUT(stderr, "Destroying results queue %p\n", q);

    if (!q)
        return;

    // Ensure it's fully drained before destroying the queue
    hts_tpool_process_reset(q, 0);
    pthread_mutex_lock(&q->p->pool_m);
    hts_tpool_process_detach(q->p, q);
    hts_tpool_process_shutdown(q);

    // Maybe a worker is scanning this queue, so delay destruction
    if (--q->ref_count > 0) {
        pthread_mutex_unlock(&q->p->pool_m);
        return;
    }

    pthread_cond_destroy(&q->output_avail_c);
    pthread_cond_destroy(&q->input_not_full_c);
    pthread_cond_destroy(&q->input_empty_c);
    pthread_cond_destroy(&q->none_processing_c);
    pthread_mutex_unlock(&q->p->pool_m);

    free(q);

    DBG_OUT(stderr, "Destroyed results queue %p\n", q);
}


/*
 * Attach and detach a thread process-queue with / from the thread pool
 * scheduler.
 *
 * We need to do attach after making a thread process, but may also wish
 * to temporarily detach if we wish to stop running jobs on a specific
 * process while permitting other process to continue.
 */
void hts_tpool_process_attach(hts_tpool *p, hts_tpool_process *q) {
    pthread_mutex_lock(&p->pool_m);
    if (p->q_head) {
        q->next = p->q_head;
        q->prev = p->q_head->prev;
        p->q_head->prev->next = q;
        p->q_head->prev = q;
    } else {
        q->next = q;
        q->prev = q;
    }
    p->q_head = q;
    assert(p->q_head && p->q_head->prev && p->q_head->next);
    pthread_mutex_unlock(&p->pool_m);
}

void hts_tpool_process_detach(hts_tpool *p, hts_tpool_process *q) {
    pthread_mutex_lock(&p->pool_m);
    if (!p->q_head || !q->prev || !q->next)
        goto done;

    hts_tpool_process *curr = p->q_head, *first = curr;
    do {
        if (curr == q) {
            q->next->prev = q->prev;
            q->prev->next = q->next;
            p->q_head = q->next;
            q->next = q->prev = NULL;

            // Last one
            if (p->q_head == q)
                p->q_head = NULL;
            break;
        }

        curr = curr->next;
    } while (curr != first);

 done:
    pthread_mutex_unlock(&p->pool_m);
}


/* ----------------------------------------------------------------------------
 * The thread pool.
 */

#define TDIFF(t2,t1) ((t2.tv_sec-t1.tv_sec)*1000000 + t2.tv_usec-t1.tv_usec)

/*
 * A worker thread.
 *
 * Once woken, each thread checks each process-queue in the pool in turn,
 * looking for input jobs that also have room for the output (if it requires
 * storing).  If found, we execute it and repeat.
 *
 * If we checked all input queues and find no such job, then we wait until we
 * are signalled to check again.
 */
static void *tpool_worker(void *arg) {
    hts_tpool_worker *w = (hts_tpool_worker *)arg;
    hts_tpool *p = w->p;
    hts_tpool_job *j;

    for (;;) {
        // Pop an item off the pool queue
        pthread_mutex_lock(&p->pool_m);

        assert(p->q_head == 0 || (p->q_head->prev && p->q_head->next));

        int work_to_do = 0;
        hts_tpool_process *first = p->q_head, *q = first;
        do {
            if (p->shutdown)
                break;

            // Iterate over queues, finding one with jobs and also
            // room to put the result.
            //if (q && q->input_head && !hts_tpool_process_output_full(q)) {
            if (q && q->input_head && q->qsize - q->n_output > p->tsize - p->nwaiting) {
                //printf("Work\n");
                work_to_do = 1;
                break;
            }

            if (q) q = q->next;
        } while (q && q != first);

        if (p->shutdown) {
        shutdown:
#ifdef DEBUG
            fprintf(stderr, "%d: Shutting down\n", worker_id(p));
#endif
            pthread_mutex_unlock(&p->pool_m);
            return NULL;
        }

        if (!work_to_do) {
            // We scanned all queues and cannot process any, so we wait.
            p->nwaiting++;

            // Push this thread to the top of the waiting stack
            if (p->t_stack_top == -1 || p->t_stack_top > w->idx)
                p->t_stack_top = w->idx;

            p->t_stack[w->idx] = 1;
//            printf("%2d: no work.  In=%d Proc=%d Out=%d  full=%d\n",
//                   w->idx, p->q_head->n_input, p->q_head->n_processing, p->q_head->n_output,
//                   hts_tpool_process_output_full(p->q_head));
            pthread_cond_wait(&w->pending_c, &p->pool_m);
            p->t_stack[w->idx] = 0;

            /* Find new t_stack_top */
            int i;
            p->t_stack_top = -1;
            for (i = 0; i < p->tsize; i++) {
                if (p->t_stack[i]) {
                    p->t_stack_top = i;
                    break;
                }
            }

            p->nwaiting--;
            pthread_mutex_unlock(&p->pool_m);
            continue; // To outer for(;;) loop.
        }

        // Otherwise work_to_do, so process as many items in this queue as
        // possible before switching to another queue.  This means threads
        // often end up being dedicated to one type of work.
        q->ref_count++;
        while (q->input_head && q->qsize - q->n_output > q->n_processing) {
            if (p->shutdown)
                goto shutdown;

            j = q->input_head;
            assert(j->p == p);

            if (!(q->input_head = j->next))
                q->input_tail = NULL;

            // Transitioning from full queue to not-full means we can wake up
            // any blocked dispatch threads.  We broadcast this as it's only
            // happening once (on the transition) rather than every time we
            // are below qsize.
            // (I wish I could remember why io_lib rev 3660 changed this from
            //  == to >=, but keeping it just incase!)
            q->n_processing++;
            if (q->n_input-- >= q->qsize)
                pthread_cond_broadcast(&q->input_not_full_c);

            if (q->n_input == 0)
                pthread_cond_signal(&q->input_empty_c);

            p->njobs--; // Total number of jobs; used to adjust to CPU scaling

            pthread_mutex_unlock(&p->pool_m);

            DBG_OUT(stderr, "%d: Processing queue %p, serial %"PRId64"\n",
                    worker_id(j->p), q, j->serial);

            hts_tpool_add_result(j, j->func(j->arg));
            //memset(j, 0xbb, sizeof(*j));
            free(j);

            pthread_mutex_lock(&p->pool_m);
        }
        if (--q->ref_count == 0) // we were the last user
            hts_tpool_process_destroy(q);
        else
            // Out of jobs on this queue, so restart search from next one.
            // This is equivalent to "work-stealing".
            p->q_head = q->next;

        pthread_mutex_unlock(&p->pool_m);
    }
}

static void wake_next_worker(hts_tpool_process *q, int locked) {
    hts_tpool *p = q->p;
    if (!locked)
        pthread_mutex_lock(&p->pool_m);

    // Update the q_head to be this queue so we'll start processing
    // the queue we know to have results.
    assert(q->prev && q->next); // attached
    p->q_head = q;

    // Wake up if we have more jobs waiting than CPUs. This partially combats
    // CPU frequency scaling effects.  Starting too many threads and then
    // running out of jobs can cause each thread to have lots of start/stop
    // cycles, which then translates often to CPU frequency scaling
    // adjustments.  Instead it is better to only start as many threads as we
    // need to keep the throughput up, meaning some threads run flat out and
    // others are idle.
    //
    // This isn't perfect as we need to know how many can actually start,
    // rather than how many are waiting.  A limit on output queue size makes
    // these two figures different.
    assert(p->njobs >= q->n_input);

    int running = p->tsize - p->nwaiting;
    int sig = p->t_stack_top >= 0 && p->njobs > p->tsize - p->nwaiting
        && (!q || q->n_processing < q->qsize - q->n_output);

//#define AVG_USAGE
#ifdef AVG_USAGE
    // Track average number of running threads and try to keep close.
    // We permit this to change, but slowly.  This avoids "boom and bust" cycles
    // where we read a lot of data, start a lot of jobs, then become idle again.
    // This way some threads run steadily and others dormant, which is better
    // for throughput.
    //
    // It's 50:50 if this is a good thing.  It helps some tasks quite significantly
    // while slightly hindering other (perhaps more usual) jobs.

    if (++p->n_count == 256) {
        p->n_count >>= 1;
        p->n_running >>= 1;
    }
    p->n_running += running;
    // Built in lag to avoid see-sawing.  Is this safe in all cases?
    if (sig && p->n_count >= 128 && running*p->n_count > p->n_running+1) sig=0;
#endif

    if (0) {
        printf("%d waiting, %d running, %d output, %d, arun %d => %d\t", p->njobs,
               running, q->n_output, q->qsize - q->n_output,
               p->n_running/p->n_count, sig);
        int i;
        for (i = 0; i < p->tsize; i++)
            putchar("x "[p->t_stack[i]]);
        putchar('\n');
    }

    if (sig)
        pthread_cond_signal(&p->t[p->t_stack_top].pending_c);

    if (!locked)
        pthread_mutex_unlock(&p->pool_m);
}

/*
 * Creates a worker pool with n worker threads.
 *
 * Returns pool pointer on success;
 *         NULL on failure
 */
hts_tpool *hts_tpool_init(int n) {
    int i;
    hts_tpool *p = malloc(sizeof(*p));
    p->tsize = n;
    p->njobs = 0;
    p->nwaiting = 0;
    p->shutdown = 0;
    p->q_head = NULL;
    p->t_stack = NULL;
    p->n_count = 0;
    p->n_running = 0;
    p->t = malloc(n * sizeof(p->t[0]));

    pthread_mutexattr_t attr;
    pthread_mutexattr_init(&attr);
    pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
    pthread_mutex_init(&p->pool_m, &attr);
    pthread_mutexattr_destroy(&attr);

    if (!(p->t_stack = malloc(n * sizeof(*p->t_stack))))
        return NULL;
    p->t_stack_top = -1;

    pthread_mutex_lock(&p->pool_m);

    for (i = 0; i < n; i++) {
        hts_tpool_worker *w = &p->t[i];
        p->t_stack[i] = 0;
        w->p = p;
        w->idx = i;
        pthread_cond_init(&w->pending_c, NULL);
        if (0 != pthread_create(&w->tid, NULL, tpool_worker, w)) {
            pthread_mutex_unlock(&p->pool_m);
            return NULL;
        }
    }

    pthread_mutex_unlock(&p->pool_m);

    return p;
}

/*
 * Returns the number of requested threads for a pool.
 */
int hts_tpool_size(hts_tpool *p) {
    return p->tsize;
}

/*
 * Adds an item to the work pool.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int hts_tpool_dispatch(hts_tpool *p, hts_tpool_process *q,
                    void *(*func)(void *arg), void *arg) {
    return hts_tpool_dispatch2(p, q, func, arg, 0);
}

/*
 * As above but optional non-block flag.
 *
 * nonblock  0 => block if input queue is full
 * nonblock +1 => don't block if input queue is full, but do not add task
 * nonblock -1 => add task regardless of whether queue is full (over-size)
 */
int hts_tpool_dispatch2(hts_tpool *p, hts_tpool_process *q,
                     void *(*func)(void *arg), void *arg, int nonblock) {
    hts_tpool_job *j;

    pthread_mutex_lock(&p->pool_m);

    DBG_OUT(stderr, "Dispatching job for queue %p, serial %"PRId64"\n",
            q, q->curr_serial);

    if (q->n_input >= q->qsize && nonblock == 1) {
        pthread_mutex_unlock(&p->pool_m);
        errno = EAGAIN;
        return -1;
    }

    if (!(j = malloc(sizeof(*j)))) {
        pthread_mutex_unlock(&p->pool_m);
        return -1;
    }
    j->func = func;
    j->arg = arg;
    j->next = NULL;
    j->p = p;
    j->q = q;
    j->serial = q->curr_serial++;

    if (nonblock == 0) {
        while (q->n_input >= q->qsize && !q->shutdown && !q->wake_dispatch)
            pthread_cond_wait(&q->input_not_full_c, &q->p->pool_m);
        if (q->shutdown) {
            free(j);
            pthread_mutex_unlock(&p->pool_m);
            return -1;
        }
        if (q->wake_dispatch) {
            //fprintf(stderr, "Wake => non-block for this operation\n");
            q->wake_dispatch = 0;
        }
    }

    p->njobs++;    // total across all queues
    q->n_input++;  // queue specific

    if (q->input_tail) {
        q->input_tail->next = j;
        q->input_tail = j;
    } else {
        q->input_head = q->input_tail = j;
    }

    DBG_OUT(stderr, "Dispatched (serial %"PRId64")\n", j->serial);

    // Let a worker know we have data.
    // Keep incoming queue at 1 per running thread, so there is always
    // something waiting when they end their current task.  If we go above
    // this signal to start more threads (if available). This has the effect
    // of concentrating jobs to fewer cores when we are I/O bound, which in
    // turn benefits systems with auto CPU frequency scaling.
    if (!q->shutdown)
        wake_next_worker(q, 1);

    pthread_mutex_unlock(&p->pool_m);

    return 0;
}

/*
 * Wakes up a single thread stuck in dispatch and make it return with
 * errno EAGAIN.
 */
void hts_tpool_wake_dispatch(hts_tpool_process *q) {
    pthread_mutex_lock(&q->p->pool_m);
    q->wake_dispatch = 1;
    pthread_cond_signal(&q->input_not_full_c);
    pthread_mutex_unlock(&q->p->pool_m);
}

/*
 * Flushes the process-queue, but doesn't exit. This simply drains the queue
 * and ensures all worker threads have finished their current tasks
 * associated with this process.
 *
 * NOT: This does not mean the worker threads are not executing jobs in
 * another process-queue.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
int hts_tpool_process_flush(hts_tpool_process *q) {
     int i;
     hts_tpool *p = q->p;

     DBG_OUT(stderr, "Flushing pool %p\n", p);

     // Drains the queue
     pthread_mutex_lock(&p->pool_m);

     // Wake up everything for the final sprint!
     for (i = 0; i < p->tsize; i++)
         if (p->t_stack[i])
             pthread_cond_signal(&p->t[i].pending_c);

     // Ensure there is room for the final sprint.
     // Shouldn't be possible to get here, but just incase.
     if (q->qsize < q->n_output + q->n_input + q->n_processing)
         q->qsize = q->n_output + q->n_input + q->n_processing;

     // Wait for n_input and n_processing to hit zero.
     while (q->n_input || q->n_processing) {
         while (q->n_input)
             pthread_cond_wait(&q->input_empty_c, &p->pool_m);
         if (q->shutdown) break;
         while (q->n_processing)
             pthread_cond_wait(&q->none_processing_c, &p->pool_m);
         if (q->shutdown) break;
    }

     pthread_mutex_unlock(&p->pool_m);

     DBG_OUT(stderr, "Flushed complete for pool %p, queue %p\n", p, q);

     return 0;
}

/*
 * Resets a process to the initial state.
 *
 * This removes any queued up input jobs, disables any notification of
 * new results/output, flushes what is left and then discards any
 * queued output.  Anything consumer stuck in a wait on results to
 * appear should stay stuck and will only wake up when new data is
 * pushed through the queue.
 *
 * Returns 0 on success;
 *        -1 on failure
 */
int hts_tpool_process_reset(hts_tpool_process *q, int free_results) {
    pthread_mutex_lock(&q->p->pool_m);
    // prevent next_result from returning data during our flush
    q->next_serial = INT_MAX;

    // Purge any queued input not yet being acted upon
    hts_tpool_job *j, *jn;
    for (j = q->input_head; j; j = jn) {
        //fprintf(stderr, "Discard input %d\n", j->serial);
        jn = j->next;
        free(j);
    }
    q->input_head = q->input_tail = NULL;
    q->n_input = 0;

    // Purge any queued output, thus ensuring we have room to flush.
    hts_tpool_result *r, *rn;
    for (r = q->output_head; r; r = rn) {
        //fprintf(stderr, "Discard output %d\n", r->serial);
        rn = r->next;
        hts_tpool_delete_result(r, free_results);
    }
    q->output_head = q->output_tail = NULL;
    q->n_output = 0;
    pthread_mutex_unlock(&q->p->pool_m);

    // Wait for any jobs being processed to complete.
    // (TODO: consider how to cancel any currently processing jobs.
    // Probably this is too hard.)
    if (hts_tpool_process_flush(q) != 0)
        return -1;

    // Discard any new output.
    pthread_mutex_lock(&q->p->pool_m);
    for (r = q->output_head; r; r = rn) {
        //fprintf(stderr, "Discard output %d\n", r->serial);
        rn = r->next;
        hts_tpool_delete_result(r, free_results);
    }
    q->output_head = q->output_tail = NULL;
    q->n_output = 0;

    // Finally reset the serial back to the starting point.
    q->next_serial = q->curr_serial = 0;
    pthread_cond_signal(&q->input_not_full_c);
    pthread_mutex_unlock(&q->p->pool_m);

    return 0;
}

/* Returns the process queue size */
int hts_tpool_process_qsize(hts_tpool_process *q) {
    return q->qsize;
}

/*
 * Destroys a thread pool.  The threads are joined into the main
 * thread so they will finish their current work load.
 */
void hts_tpool_destroy(hts_tpool *p) {
    int i;

    DBG_OUT(stderr, "Destroying pool %p\n", p);

    /* Send shutdown message to worker threads */
    pthread_mutex_lock(&p->pool_m);
    p->shutdown = 1;

    DBG_OUT(stderr, "Sending shutdown request\n");

    for (i = 0; i < p->tsize; i++)
        pthread_cond_signal(&p->t[i].pending_c);

    pthread_mutex_unlock(&p->pool_m);

    DBG_OUT(stderr, "Shutdown complete\n");

    for (i = 0; i < p->tsize; i++)
        pthread_join(p->t[i].tid, NULL);

    pthread_mutex_destroy(&p->pool_m);
    for (i = 0; i < p->tsize; i++)
        pthread_cond_destroy(&p->t[i].pending_c);

    if (p->t_stack)
        free(p->t_stack);

    free(p->t);
    free(p);

    DBG_OUT(stderr, "Destroyed pool %p\n", p);
}


/*
 * Destroys a thread pool without waiting on jobs to complete.
 * Use hts_tpool_kill(p) to quickly exit after a fatal error.
 */
void hts_tpool_kill(hts_tpool *p) {
    int i;

    DBG_OUT(stderr, "Destroying pool %p, kill=%d\n", p, kill);

    for (i = 0; i < p->tsize; i++)
        pthread_kill(p->t[i].tid, SIGINT);

    pthread_mutex_destroy(&p->pool_m);
    for (i = 0; i < p->tsize; i++)
        pthread_cond_destroy(&p->t[i].pending_c);

    if (p->t_stack)
        free(p->t_stack);

    free(p->t);
    free(p);

    DBG_OUT(stderr, "Destroyed pool %p\n", p);
}


/*=============================================================================
 * Test app.
 *
 * This can be considered both as a basic test and as a worked example for
 * various usage patterns.
 *=============================================================================
 */

#ifdef TEST_MAIN

#include <stdio.h>

#ifndef TASK_SIZE
#define TASK_SIZE 1000
#endif

/*-----------------------------------------------------------------------------
 * Unordered x -> x*x test.
 * Results arrive in order of completion.
 */
void *doit_square_u(void *arg) {
    int job = *(int *)arg;

    usleep(random() % 100000); // to coerce job completion out of order

    printf("RESULT: %d\n", job*job);

    free(arg);
    return NULL;
}

int test_square_u(int n) {
    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q = hts_tpool_process_init(p, n*2, 1);
    int i;

    // Dispatch jobs
    for (i = 0; i < TASK_SIZE; i++) {
        int *ip = malloc(sizeof(*ip));
        *ip = i;
        hts_tpool_dispatch(p, q, doit_square_u, ip);
    }

    hts_tpool_process_flush(q);
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    return 0;
}


/*-----------------------------------------------------------------------------
 * Ordered x -> x*x test.
 * Results arrive in numerical order.
 *
 * This implementation uses a non-blocking dispatch to avoid dead-locks
 * where one job takes too long to complete.
 */
void *doit_square(void *arg) {
    int job = *(int *)arg;
    int *res;

    // One excessively slow, to stress test output queue filling and
    // excessive out of order scenarios.
    usleep(500000 * ((job&31)==31) + random() % 10000);

    res = malloc(sizeof(*res));
    *res = (job<0) ? -job*job : job*job;

    free(arg);
    return res;
}

int test_square(int n) {
    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q = hts_tpool_process_init(p, n*2, 0);
    int i;
    hts_tpool_result *r;

    // Dispatch jobs
    for (i = 0; i < TASK_SIZE; i++) {
        int *ip = malloc(sizeof(*ip));
        *ip = i;
        int blk;

        do {
            // In the situation where some jobs take much longer than
            // others, we could end up blocking here as we haven't got
            // any room in the output queue to place it. (We don't launch a
            // job if the output queue is full.)

            // This happens when the next serial number to fetch is, eg, 50
            // but jobs 51-100 have all executed really fast and appeared in
            // the output queue before 50.  A dispatch & check-results
            // alternating loop can fail to find job 50 many times over until
            // eventually the dispatch blocks before it arrives.

            // Our solution is to dispatch in non-blocking mode so we are
            // always to either dispatch or consume a result.
            blk = hts_tpool_dispatch2(p, q, doit_square, ip, 1);

            // Check for results.
            if ((r = hts_tpool_next_result(q))) {
                printf("RESULT: %d\n", *(int *)hts_tpool_result_data(r));
                hts_tpool_delete_result(r, 1);
            }
            if (blk == -1) {
                // The alternative is a separate thread for dispatching and/or
                // consumption of results. See test_squareB.
                putchar('.'); fflush(stdout);
                usleep(10000);
            }
        } while (blk == -1);
    }

    // Wait for any input-queued up jobs or in-progress jobs to complete.
    hts_tpool_process_flush(q);

    while ((r = hts_tpool_next_result(q))) {
        printf("RESULT: %d\n", *(int *)hts_tpool_result_data(r));
        hts_tpool_delete_result(r, 1);
    }

    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    return 0;
}

/*-----------------------------------------------------------------------------
 * Ordered x -> x*x test.
 * Results arrive in numerical order.
 *
 * This implementation uses separate dispatching threads and job consumption
 * threads (main thread).  This means it can use a blocking calls for
 * simplicity elsewhere.
 */
struct squareB_opt {
    hts_tpool *p;
    hts_tpool_process *q;
    int n;
};
static void *test_squareB_dispatcher(void *arg) {
    struct squareB_opt *o = (struct squareB_opt *)arg;
    int i, *ip;

    for (i = 0; i < o->n; i++) {
        ip = malloc(sizeof(*ip));
        *ip = i;

        hts_tpool_dispatch(o->p, o->q, doit_square, ip);
    }

    // Dispatch an sentinel job to mark the end
    *(ip = malloc(sizeof(*ip))) = -1;
    hts_tpool_dispatch(o->p, o->q, doit_square, ip);
    pthread_exit(NULL);
}

int test_squareB(int n) {
    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q = hts_tpool_process_init(p, n*2, 0);
    struct squareB_opt o = {p, q, TASK_SIZE};
    pthread_t tid;

    // Launch our job creation thread.
    pthread_create(&tid, NULL, test_squareB_dispatcher, &o);

    // Consume all results until we find the end-of-job marker.
    for(;;) {
        hts_tpool_result *r = hts_tpool_next_result_wait(q);
        int x = *(int *)hts_tpool_result_data(r);
        hts_tpool_delete_result(r, 1);
        if (x == -1)
            break;
        printf("RESULT: %d\n", x);
    }

    // Wait for any input-queued up jobs or in-progress jobs to complete.
    // This should do nothing as we've been executing until the termination
    // marker of -1.
    hts_tpool_process_flush(q);
    assert(hts_tpool_next_result(q) == NULL);

    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
    pthread_join(tid, NULL);

    return 0;
}


/*-----------------------------------------------------------------------------
 * A simple pipeline test.
 * We use a dedicated input thread that does the initial generation of job
 * and dispatch, several execution steps running in a shared pool, and a
 * dedicated output thread that prints up the final result.  It's key that our
 * pipeline execution stages can run independently and don't themselves have
 * any waits.  To achieve this we therefore also use some dedicated threads
 * that take the output from one queue and resubmits the job as the input to
 * the next queue.
 *
 * More generally this could perhaps be a single pipeline thread that
 * marshalls multiple queues and their interactions, but this is simply a
 * demonstration of a single pipeline.
 *
 * Our process fills out the bottom byte of a 32-bit int and then shifts it
 * left one byte at a time.  Only the final stage needs to be ordered.  Each
 * stage uses its own queue.
 *
 * Possible improvement: we only need the last stage to be ordered.  By
 * allocating our own serial numbers for the first job and manually setting
 * these serials in the last job, perhaps we can permit out of order execution
 * of all the inbetween stages.  (I doubt it'll affect speed much though.)
 */

static void *pipe_input_thread(void *arg);
static void *pipe_stage1(void *arg);
static void *pipe_stage2(void *arg);
static void *pipe_stage3(void *arg);
static void *pipe_output_thread(void *arg);

typedef struct {
    hts_tpool *p;
    hts_tpool_process *q1;
    hts_tpool_process *q2;
    hts_tpool_process *q3;
    int n;
} pipe_opt;

typedef struct {
    pipe_opt *o;
    unsigned int x;
    int eof; // set with last job.
} pipe_job;

static void *pipe_input_thread(void *arg) {
    pipe_opt *o = (pipe_opt *)arg;

    int i;
    for (i = 1; i <= o->n; i++) {
        pipe_job *j = malloc(sizeof(*j));
        j->o = o;
        j->x = i;
        j->eof = (i == o->n);

        printf("I  %08x\n", j->x);

        if (hts_tpool_dispatch(o->p, o->q1, pipe_stage1, j) != 0) {
            free(j);
            pthread_exit((void *)1);
        }
    }

    pthread_exit(NULL);
}

static void *pipe_stage1(void *arg) {
    pipe_job *j = (pipe_job *)arg;

    j->x <<= 8;
    usleep(random() % 10000); // fast job
    printf("1  %08x\n", j->x);

    return j;
}

static void *pipe_stage1to2(void *arg) {
    pipe_opt *o = (pipe_opt *)arg;
    hts_tpool_result *r;

    while ((r = hts_tpool_next_result_wait(o->q1))) {
        pipe_job *j = (pipe_job *)hts_tpool_result_data(r);
        hts_tpool_delete_result(r, 0);
        if (hts_tpool_dispatch(j->o->p, j->o->q2, pipe_stage2, j) != 0)
            pthread_exit((void *)1);
        if (j->eof)
            break;
    }

    pthread_exit(NULL);
}

static void *pipe_stage2(void *arg) {
    pipe_job *j = (pipe_job *)arg;

    j->x <<= 8;
    usleep(random() % 100000); // slow job
    printf("2  %08x\n", j->x);

    return j;
}

static void *pipe_stage2to3(void *arg) {
    pipe_opt *o = (pipe_opt *)arg;
    hts_tpool_result *r;

    while ((r = hts_tpool_next_result_wait(o->q2))) {
        pipe_job *j = (pipe_job *)hts_tpool_result_data(r);
        hts_tpool_delete_result(r, 0);
        if (hts_tpool_dispatch(j->o->p, j->o->q3, pipe_stage3, j) != 0)
            pthread_exit((void *)1);
        if (j->eof)
            break;
    }

    pthread_exit(NULL);
}

static void *pipe_stage3(void *arg) {
    pipe_job *j = (pipe_job *)arg;

    usleep(random() % 10000); // fast job
    j->x <<= 8;
    return j;
}

static void *pipe_output_thread(void *arg) {
    pipe_opt *o = (pipe_opt *)arg;
    hts_tpool_result *r;

    while ((r = hts_tpool_next_result_wait(o->q3))) {
        pipe_job *j = (pipe_job *)hts_tpool_result_data(r);
        int eof = j->eof;
        printf("O  %08x\n", j->x);
        hts_tpool_delete_result(r, 1);
        if (eof)
            break;
    }

    pthread_exit(NULL);
}

int test_pipe(int n) {
    hts_tpool *p = hts_tpool_init(n);
    hts_tpool_process *q1 = hts_tpool_process_init(p, n*2, 0);
    hts_tpool_process *q2 = hts_tpool_process_init(p, n*2, 0);
    hts_tpool_process *q3 = hts_tpool_process_init(p, n*2, 0);
    pipe_opt o = {p, q1, q2, q3, TASK_SIZE};
    pthread_t tidIto1, tid1to2, tid2to3, tid3toO;
    void *retv;
    int ret;

    // Launch our data source and sink threads.
    pthread_create(&tidIto1, NULL, pipe_input_thread,  &o);
    pthread_create(&tid1to2, NULL, pipe_stage1to2,     &o);
    pthread_create(&tid2to3, NULL, pipe_stage2to3,     &o);
    pthread_create(&tid3toO, NULL, pipe_output_thread, &o);

    // Wait for tasks to finish.
    ret = 0;
    pthread_join(tidIto1, &retv); ret |= (retv != NULL);
    pthread_join(tid1to2, &retv); ret |= (retv != NULL);
    pthread_join(tid2to3, &retv); ret |= (retv != NULL);
    pthread_join(tid3toO, &retv); ret |= (retv != NULL);
    printf("Return value %d\n", ret);

    hts_tpool_process_destroy(q1);
    hts_tpool_process_destroy(q2);
    hts_tpool_process_destroy(q3);
    hts_tpool_destroy(p);

    return 0;
}

/*-----------------------------------------------------------------------------*/
int main(int argc, char **argv) {
    int n;
    srandom(0);

    if (argc < 3) {
        fprintf(stderr, "Usage: %s command n_threads\n", argv[0]);
        fprintf(stderr, "Where commands are:\n\n");
        fprintf(stderr, "unordered       # Unordered output\n");
        fprintf(stderr, "ordered1        # Main thread with non-block API\n");
        fprintf(stderr, "ordered2        # Dispatch thread, blocking API\n");
        fprintf(stderr, "pipe            # Multi-stage pipeline, several queues\n");
        exit(1);
    }

    n = atoi(argv[2]);
    if (strcmp(argv[1], "unordered") == 0) return test_square_u(n);
    if (strcmp(argv[1], "ordered1") == 0)  return test_square(n);
    if (strcmp(argv[1], "ordered2") == 0)  return test_squareB(n);
    if (strcmp(argv[1], "pipe") == 0)      return test_pipe(n);

    fprintf(stderr, "Unknown sub-command\n");
    exit(1);
}
#endif
