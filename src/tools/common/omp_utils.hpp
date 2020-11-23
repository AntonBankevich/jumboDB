//
// Created by anton on 7/27/20.
//
#pragma once
#include "logging.hpp"
#include <parallel/algorithm>
#include <omp.h>
#include <utility>
#include <numeric>
#include <wait.h>


class ParallelCounter {
    std::vector<size_t> cnt;
public:
    explicit ParallelCounter(size_t thread_num) : cnt(thread_num){
    }

    void operator++() {
        cnt[omp_get_thread_num()]++;
    }

    void operator+=(const size_t val) {
        cnt[omp_get_thread_num()] += val;
    }

    size_t get() const {
        return std::accumulate(cnt.begin(), cnt.end(), size_t(0));
    }
};

template<class T>
class ParallelRecordCollector {
    std::vector<std::vector<T>> recs;
public:
    friend class Iterator;
    class Iterator : public std::iterator<std::forward_iterator_tag, T, size_t,  T*, T&>{
    private:
        ParallelRecordCollector<T> &data;
        size_t row;
        size_t col;
    public:
        explicit Iterator(ParallelRecordCollector<T> &_data, size_t _row = 0, size_t _col = 0) : data(_data), row(_row), col(_col) {
            while(row < data.recs.size() && col == data.recs[row].size()) {
                row += 1;
                col = 0;
            }
        }

        void operator++() {
            col += 1;
            while(row < data.recs.size() && col == data.recs[row].size()) {
                row += 1;
                col = 0;
            }
        }

        T &operator *() {
            return data.recs[row][col];
        }

        bool operator==(const Iterator &other) {
            return row == other.row && col == other.col;
        }
        bool operator!=(const Iterator &other) {
            return row != other.row || col != other.col;
        }

    };
    explicit ParallelRecordCollector(size_t thread_num) : recs(thread_num){
    }

    void add(const T &rec) {
        recs[omp_get_thread_num()].emplace_back(rec);
    }

    template<class I>
    void addAll(I begin, I end) {
        recs[omp_get_thread_num()].insert(recs[omp_get_thread_num()].end(), begin, end);
    }

    template< class... Args >
    void emplace_back( Args&&... args ) {
        recs[omp_get_thread_num()].emplace_back(args...);
    }

    Iterator begin() {
        return Iterator(*this, 0, 0);
    }

    Iterator end() {
        return Iterator(*this, recs.size(), 0);
    }

    size_t size() const {
        size_t res = 0;
        for (const std::vector<T> & row : recs) {
            res += row.size();
        }
        return res;
    }

    bool empty() const {
        return size() == 0;
    }

    std::vector<T> collect() {
        std::vector<T> res;
        for(std::vector<T> &row : recs) {
            res.insert(res.end(), row.begin(), row.end());
            row.clear();
        }
        return std::move(res);
    }

    void clear() {
        for(std::vector<T> &row : recs) {
            row.clear();
        }
    }

    std::vector<T> collectUnique() {
        std::vector<T> res = collect();
        __gnu_parallel::sort(res.begin(), res.end());
        res.erase(std::unique(res.begin(), res.end()), res.end());
        return std::move(res);
    }
};

template<class T>
std::ostream& operator<<(std::ostream& out, const ParallelRecordCollector<T>& tree) {
    if(tree.size() == 0) {
        return out << "[]" << std::endl;
    }
    out << "[";
    for(const T & item: tree) {
        out << item << ", ";
    }
    return out << "]";
}


template<class V>
class ParallelProcessor {
public:
    std::function<void(V &)> task = [] (V &) {};
    std::function<void ()> doBefore = [] () {};
    std::function<void ()> doAfter = [] () {};
    std::function<void ()> doInParallel = [] () {};
    std::function<void (V&)> doInOneThread = [] (V &) {};
    std::function<void ()> doInTheEnd = [] () {};
    logging::Logger &logger;
    size_t threads;

    ParallelProcessor(std::function<void(V &)> _task, logging::Logger & _logger, size_t _threads) :
                    task(_task), logger(_logger), threads(_threads) {
    }

//This method expects iterator to be a generator, i.e. it returns temporary objects. Thus we have to store them in a buffer and
//keep track of total size of stored objects.
    template<class I>
    void processRecords(I begin, I end, size_t bucket_length = 1024 * 1024) {
        logger.info() << "Starting parallel calculation" << std::endl;
        omp_set_num_threads(threads);
//        size_t bucket_length = 1024 * 1024;
        size_t buffer_size = 1024 * 1024;
        size_t max_length = 1024 * 1024 * 1024;
        ParallelProcessor<V> &self = *this;
        size_t total = 0;
        size_t total_len = 0;
        while(begin != end) {
            size_t clen = 0;
            std::vector<V> items;
            items.reserve(buffer_size);
            doBefore();
#pragma omp parallel default(none) shared(begin, end, items, buffer_size, clen, max_length, bucket_length, self, std::cout)
            {
#pragma omp single
                {
#pragma omp task default(none) shared(self, std::cout)
                    {
                        self.doInParallel();
                    }
                    while (begin != end && items.size() < buffer_size && clen < max_length) {
                        size_t left = items.size();
                        size_t right = items.size();
                        size_t cur_length = 0;
                        while (begin != end && items.size() < buffer_size && clen < max_length && cur_length < bucket_length) {
                            items.emplace_back(*begin);
                            ++begin;
                            right += 1;
                            V &item = items.back();
                            clen += item.size();
                            cur_length += item.size();
                            self.doInOneThread(item);
                        }
#pragma omp task default(none) shared(items, self, std::cout) firstprivate(left, right)
                        {
                            for(size_t i = left; i < right; i++)
                                self.task(items[i]);
                        }
                    }
                }
            }
            doAfter();
            logger.info() << items.size() << " items of total length "<< clen << " processed " << std::endl;
            total += items.size();
            items.clear();
            total_len += clen;
        }
        doInTheEnd();
        logger.info() << "Finished parallel processing. Processed " << total <<
               " items with total length " << total_len << std::endl;
    }


    //This method expects that iterators return references to objects instead of temporary objects.
    template<class I>
    void processObjects(I begin, I end, size_t bucket_size = 1024) {
        logger.info() << "Starting parallel calculation" << std::endl;
        omp_set_num_threads(threads);
        ParallelProcessor<V> &self = *this;
        size_t buffer_size = 1024 * 1024;
        size_t total = 0;
        while(begin != end) {
            std::vector<V*> items;
            items.reserve(buffer_size);
            doBefore();
#pragma omp parallel default(none) shared(begin, end, items, buffer_size, bucket_size, self)
            {
#pragma omp single
                {
#pragma omp task default(none) shared(self)
                    {
                        self.doInParallel();
                    }
                    while (begin != end && items.size() < buffer_size) {
                        size_t left = items.size();
                        size_t right = items.size();
                        while (begin != end && items.size() < buffer_size && right - left < bucket_size) {
                            items.push_back(&(*begin));
                            self.doInOneThread(*items.back());
                            ++begin;
                            right += 1;
                        }
#pragma omp task default(none) shared(items, self) firstprivate(left, right)
                        {
                            for(size_t i = left; i < right; i++)
                                self.task(*items[i]);
                        }
                    }
                }
            }
            doAfter();
            logger.info() << "Processed " << items.size() << " items" << std::endl;
            total += items.size();
            items.clear();
        }
        doInTheEnd();
        logger.info() << "Finished parallel processing. Processed " << total << " items " << std::endl;
    }

};

//This method expects that iterators return references to objects instead of temporary objects.
template<class I>
void processObjects(I begin, I end, logging::Logger &logger, size_t threads, std::function<void(typename I::value_type &)> task,
                    size_t bucket_size = 1024) {
    typedef typename I::value_type V;
    ParallelProcessor<V>(task, logger, threads).processObjects(begin, end, bucket_size);
}

//This method expects iterator to be a generator, i.e. it returns temporary objects. Thus we have to store them in a buffer and
//keep track of total size of stored objects.
template<class I>
void processRecords(I begin, I end, logging::Logger &logger, size_t threads, std::function<void(typename I::value_type &)> task,
                    size_t bucket_length = 1024 * 1024) {
    typedef typename I::value_type V;
    ParallelProcessor<V>(task, logger, threads).processRecords(begin, end, bucket_length);
}

inline void runInFork(const std::function<void()>& f) {
    pid_t p = fork();
    if (p < 0) {
        std::cout << "Fork failed" << std::endl;
        exit(1);
    }
    if(p == 0) {
        f();
        exit(0);
    } else {
        int status = 0;
        waitpid(p, &status, 0);
        if (WEXITSTATUS(status) || WIFSIGNALED(status)) {
            std::cout << "Child process crashed" << std::endl;
            exit(1);
        }
    }
}