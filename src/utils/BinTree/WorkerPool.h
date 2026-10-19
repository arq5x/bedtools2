/*
 * WorkerPool.h
 *
 * The threads of intersect's -threads N: the thread that makes the pool is
 * thread 0, and N - 1 workers wait for jobs.
 */

#ifndef WORKERPOOL_H_
#define WORKERPOOL_H_

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <system_error>

class WorkerPool {
public:
	// Starts as many of the numThreads - 1 workers as the system allows.
	explicit WorkerPool(int numThreads)
	: _generation(0), _running(0), _stop(false) {
		_threads.reserve(numThreads > 1 ? numThreads - 1 : 0);
		for (int t = 1; t < numThreads; t++) {
			try {
				_threads.push_back(std::thread(&WorkerPool::loop, this, t));
			} catch (const std::system_error &) {
				break;
			}
		}
	}
	~WorkerPool() {
		{
			std::lock_guard<std::mutex> lock(_mutex);
			_stop = true;
		}
		_wake.notify_all();
		for (size_t i = 0; i < _threads.size(); i++) {
			_threads[i].join();
		}
	}
	// The threads, this one included.
	int size() const { return (int)_threads.size() + 1; }

	// Runs job(t) on each worker t (1 to size() - 1) and returns at once;
	// wait() returns once every worker has returned from it.
	void start(const std::function<void(int)> &job) {
		{
			std::lock_guard<std::mutex> lock(_mutex);
			_job = job;
			_running = (int)_threads.size();
			_generation++;
		}
		_wake.notify_all();
	}
	void wait() {
		std::unique_lock<std::mutex> lock(_mutex);
		while (_running > 0) {
			_done.wait(lock);
		}
	}
	// job(t) on every thread, this one as thread 0, returning when all are done.
	void run(const std::function<void(int)> &job) {
		start(job);
		job(0);
		wait();
	}

	WorkerPool(const WorkerPool &) = delete;
	WorkerPool &operator=(const WorkerPool &) = delete;

private:
	std::vector<std::thread> _threads;
	std::mutex _mutex;
	std::condition_variable _wake;
	std::condition_variable _done;
	std::function<void(int)> _job;
	unsigned long _generation;
	int _running;
	bool _stop;

	void loop(int t) {
		unsigned long seen = 0;
		for (;;) {
			std::function<void(int)> job;
			{
				std::unique_lock<std::mutex> lock(_mutex);
				while (!_stop && _generation == seen) {
					_wake.wait(lock);
				}
				if (_stop) {
					return;
				}
				seen = _generation;
				job = _job;
			}
			job(t);
			std::lock_guard<std::mutex> lock(_mutex);
			if (--_running == 0) {
				_done.notify_all();
			}
		}
	}
};

#endif /* WORKERPOOL_H_ */
