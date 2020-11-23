//
// Created by anton on 8/17/20.
//

#pragma once
#include "sparse_dbg.hpp"
#include <sequences/sequence.hpp>
#include <common/logging.hpp>
#include <vector>


template<typename htype>
size_t maxOutgoingCov(const Vertex<htype> &rec) {
    size_t res = 0;
    for(const Edge<htype> &edge : rec.getOutgoing())
        res = std::max(res, edge.end()->coverage());
    return res;
}

template<typename htype>
void processVertex(Vertex<htype> &rec, size_t threshold, size_t k) {
    rec.lock();
    rec.rc().lock();
    if (rec.outDeg() > 1 && rec.coverage() > threshold && maxOutgoingCov(rec) < threshold) {
        for (const Edge<htype> &edge : rec.getOutgoing()) {
            VERIFY(edge.end() != nullptr);
            std::vector<Edge<htype>> path = rec.walkForward(edge);
            Vertex<htype> &end = path.back().end()->rc();
            if (end.hash() <= rec.hash())
                continue;
            end.lock();
            end.rc().lock();
            if (end.inDeg() > 1)
                continue;
            if (end.outDeg() > 1 && maxOutgoingCov(end) <= threshold) {
                continue;
            }
            bool ok = true;
            size_t len = 0;
            for (const Edge<htype> &path_edge : path) {
                len += path_edge.size();
                if (path_edge.end()->coverage() > threshold) {
                    ok = false;
                    break;
                }
            }
            const Edge<htype> & rcStart = end.removeEdgesTo();
            end.unlock();
            end.rc().unlock();
        }
    }
    rec.unlock();
    rec.rc().unlock();
}

template<typename htype>
void simplifySparseGraph(logging::Logger & logger, SparseDBG<htype> &sdbg, size_t threshold, size_t threads) {
#pragma omp parallel default(none) shared(sdbg, threshold)
    {
#pragma omp single
        {
            for(auto it = sdbg.begin(); it != sdbg.end(); ++it) {
                Vertex<htype> &rec = it->second;
                htype hash = it->first;
                if(!rec.isJunction())
                    continue;
#pragma omp task default(none) shared(sdbg, rec, hash)
                {
                }
            }
        }
    }
}
