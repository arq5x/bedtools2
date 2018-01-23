Thread pool tests
=================
The thread_pool.c file has a built-in test program which is enabled when compiling with TEST_MAIN defined. The test program can be run in four different modes by giving a command-line parameter: unordered, ordered1, ordered2, and pipe. The modes and their expected outputs are described below.

unordered
---------
Dispatches TASK_SIZE (=1000) jobs to the thread pool and waits for them to finish. The job index (0..TASK_SIZE-1) is passed as a parameter. The job function is doit_square_u, which sleeps for a while and then prints the square of its input parameter to stdout.

Expected output when n = 1:
```
RESULT: 0
...
RESULT: 998001
```

Expected output when n > 1: same, but in jumbled up order.

ordered1
--------
Dispatches TASK_SIZE (=1000) jobs to the thread pool in non-blocking mode. Results are returned on the result queue and are pulled in order. The job index (0..TASK_SIZE-1) is passed as a parameter. The job function is doit_square, which sleeps for a while and then returns the square of its input parameter as a result. Some of the jobs take way longer than the others to finish.

The expected output is the results printed in order, regardless of n.

ordered2
--------
Starts a dispatcher thread which dispatches jobs to the thread pool. After all regular jobs have been dispatched, a sentinel job follows where the input parameter is set to -1, which receives special handling in doit_square to return the -1 as the result.

Results are consumed on the main thread using hts_tpool_next_result_wait, until the end-of-job marker is found.

The expected output is the results printed in order, regardless of n.

pipe
----
This program uses one thread pool (hts_tpool) and three queues (hts_tpool_process) shared across threads using a pipe_opt struct. There are four threads: input, stage1to2, stage2to3, and output.

The input thread (pipe_input_thread procedure) dispatches jobs to the thread pool with the job number (1..TASK_SIZE) and an end-of-job flag as parameters. The jobs are executed by the pipe_stage1 procedure, which multiplies by 256 and sleeps for a short while.

The stage1to2 thread (pipe_stage1to2 procedure) pulls results from the first queue (q1) and passes them to new jobs in the thread pool. These jobs are executed by the pipe_stage2 procedure, which does the same as pipe_stage1, only slower.

The stage2to3 thread is similar to the stage1to2 thread. It pulls from the second queue and dispatches new jobs to be executed by the pipe_stage3 procedure. pipe_stage3 is similar to pipe_stage1.

The output thread pulls from the third queue.

Expected output:
```
I 00000001
1 00000100
2 00010000
O 01000000
...
I 000003e8
1 0003e800
2 03e80000
O e8000000
```
...but not in order, because the input queues might be served in any order.

However, if only the lines from the output thread are printed, they should be in order regardless of the number of threads.
