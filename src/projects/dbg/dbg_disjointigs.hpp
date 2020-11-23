//
// Created by anton on 8/3/20.
//
#pragma once

#include "sparse_dbg.hpp"
#include <common/bloom_filter.hpp>
#include <common/cl_parser.hpp>
#include <common/simple_computation.hpp>
#include "common/logging.hpp"

template<typename htype>
Sequence buildDisjointig(const Vertex<htype> &rec,
                         const std::vector<Edge<htype>*> &path) {
    Sequence disjointig = rec.pathSeq(path);
    const Vertex<htype> &last = path.back()->end()->rc();
    const Edge<htype> &lastEdge = path.size() > 1 ? path[path.size() - 2]->end()->sparseRcEdge(*path.back()) :
            rec.sparseRcEdge(*path.back());
//    Sequence seq = rec.pathSeq(path);
//    std::string kmer = "TCTGAGCACAGTGTGACACTCACTAGAGTGAGAGCAGATGATGCTCAGTGCTCACAGAGATGACACAGCACTCAGTCACTCACAGCATGAGTGTGCAGCTGAGAGAGCATCAGAGCTGACAGTCTCTCAGCTCAGCAGCTGAGCTGACAGCATCACTCAGTATAGATGCTCACTCTGTCAGTGAGCTCTGCAGTCAGTCAGAGACTGACAGCATGACTGTAGTGCTGTCTCAGACACAGACACATCTGACATCATCACATACACAGAGCTCTACAGTCACGAGACAGACACATGTGCAGAGAGACAGAGTGCAGAGCTGATGAGTCATATCACAGACAGATCTGCACTGCTGATCACACTCATCATCACTCAGTGAGAGACTGATCACAGATGAGCATCTGCACAGAGATACTATCAGTCACTATGATAGTGTGATGACTAGCTATCGAGTGTGACTCAGTCACGTACTCTACAGCAGATGATACAGTATCAGAGCTGTGTCAGCTATCGT";
//    std::string kmer1 = "ACGATAGCTGACACAGCTCTGATACTGTATCATCTGCTGTAGAGTACGTGACTGAGTCACACTCGATAGCTAGTCATCACACTATCATAGTGACTGATAGTATCTCTGTGCAGATGCTCATCTGTGATCAGTCTCTCACTGAGTGATGATGAGTGTGATCAGCAGTGCAGATCTGTCTGTGATATGACTCATCAGCTCTGCACTCTGTCTCTCTGCACATGTGTCTGTCTCGTGACTGTAGAGCTCTGTGTATGTGATGATGTCAGATGTGTCTGTGTCTGAGACAGCACTACAGTCATGCTGTCAGTCTCTGACTGACTGCAGAGCTCACTGACAGAGTGAGCATCTATACTGAGTGATGCTGTCAGCTCAGCTGCTGAGCTGAGAGACTGTCAGCTCTGATGCTCTCTCAGCTGCACACTCATGCTGTGAGTGACTGAGTGCTGTGTCATCTCTGTGAGCACTGAGCATCATCTGCTCTCACTCTAGTGAGTGTCACACTGTGCTCAGA";
//    if(seq.str().find(kmer) != size_t(-1) || seq.str().find(kmer1) != size_t(-1)) {
//        std::cout << "Found " << &rec << " " << rec.outDeg() << " " << rec.inDeg() << " " << &last << " " << last.outDeg() << " " << last.inDeg() <<
//                        " " << seq.str().find(kmer) << " " << seq.str().find(kmer1) << std::endl <<
//                        path[0].intCov() << " " << lastEdge.intCov() << " " << rec.seq.size() << " " << disjointig.size() << std::endl
//                        << rec.seq << std::endl;
//        std::cout << seq.str() << std::endl;
//        for(size_t i = 0; i < rec.outDeg(); i++) {
//            if(rec.getOutgoing()[i].seq == path[0].seq) {
//                std::cout << "* " << seq.str().find(kmer) << " " << seq.str().find(kmer1) << std::endl;
//            }
//            Sequence tmp = rec.seq + rec.getOutgoing()[i].seq;
//            std::cout <<"out " << rec.getOutgoing()[i].intCov() << " " << tmp.str().find(kmer) << " " << tmp.str().find(kmer1) <<
//                        std::endl << rec.getOutgoing()[i].seq << std::endl;
//        }
//        for(size_t i = 0; i < last.outDeg(); i++) {
//            Sequence tmp = last.seq + last.getOutgoing()[i].seq;
//            std::cout <<"inc " << last.getOutgoing()[i].intCov() << " " << tmp.str().find(kmer) << " " << tmp.str().find(kmer1)
//                        <<std::endl << last.getOutgoing()[i].seq << std::endl;
//        }
//    }
    if(path[0]->intCov() + lastEdge.intCov() + rec.seq.size() > disjointig.size())
        return Sequence{};
    disjointig = disjointig.Subseq(path[0]->intCov(), disjointig.size() - lastEdge.intCov());
    if (rec.inDeg() > 1 && rec.outDeg() == 1) {
        VERIFY(path[0]->intCov() == 0);
        const Edge<htype>& extra = rec.rc().getOutgoing()[0];
        disjointig = !(extra.seq.Subseq(0, extra.intCov())) + disjointig;
    }
    if(last.inDeg() > 1 && last.outDeg() == 1) {
        VERIFY(lastEdge.intCov() == 0);
        const Edge<htype>& extra = last.rc().getOutgoing()[0];
        disjointig = disjointig + extra.seq.Subseq(0, extra.intCov());
    }
    return disjointig;
}

template<typename htype>
void processVertex(Vertex<htype> &rec, ParallelRecordCollector<Sequence> &res) {
    for(Edge<htype> & edge : rec.getOutgoing()) {
        VERIFY(edge.end() != nullptr);
        VERIFY(!rec.seq.empty());
        std::vector<Edge<htype>*> path = rec.walkForward(edge);
//        Vertex<htype> &rec1 = path.back().end()->rc();
//        const Edge<htype> &edge1 = path.size() > 1 ? path[path.size() - 2].end()->sparseRcEdge(path.back()) :
//                                   rec.sparseRcEdge(path.back());
//        std::vector<Edge<htype>> path1 = rec1.walkForward(edge1);
//        bool oppa1 = (rec < path.back().end()->rc() || (rec == path.back().end()->rc() && rec.pathSeq(path) <= !rec.pathSeq(path)));
//        bool oppa2 = (rec1 < path1.back().end()->rc() || (rec1 == path1.back().end()->rc() && rec1.pathSeq(path1) <= !rec1.pathSeq(path1)));
//        htype oppa = htype(63100013936723) * 1000000000000 * 1000000000000 + htype(716503108335) * 1000000000000  + 449034034478;
//        if(rec.hash() == oppa || path.back().end()->hash() == oppa ||
//        if(     (oppa1 == oppa2 && rec.pathSeq(path) != !rec.pathSeq(path)) ||
//                (!oppa1 && !oppa2 && rec.pathSeq(path) == !rec.pathSeq(path))) {
//            std::cout << "Incorrect comparison of vertices. Antisymmetry is broken " << oppa1 << " " << oppa2 <<
//                        " " << (rec.pathSeq(path) == !rec.pathSeq(path)) << std::endl;
//            std::cout << rec.hash()<< " "<< rec.isCanonical() << rec1.hash() << " " << rec1.isCanonical() << std::endl;
//            std::cout << rec.pathSeq(path) << std::endl << !rec.pathSeq(path)<< std::endl << rec1.pathSeq(path1) << std::endl;
//            VERIFY(false)
//        }
        if(rec < path.back()->end()->rc() || (rec == path.back()->end()->rc() && rec.pathSeq(path) <= !rec.pathSeq(path))) {
            Sequence disjointig = buildDisjointig(rec, path);
            for(size_t i = 0; i + 1 < path.size(); i++) {
                path[i]->end()->clearSequence();
            }
            if (!disjointig.empty())
                res.add(disjointig.copy());
        }
    }
}

template<typename htype>
void prepareVertex(Vertex<htype> &vertex) {
    vertex.sortOutgoing();
    for(size_t i = 1; i < vertex.outDeg(); i++) {
        vertex.getOutgoing()[i].incCov(vertex.getOutgoing()[i].seq.commonPrefix(vertex.getOutgoing()[i - 1].seq));
    }
    if(vertex.outDeg() > 1) {
        vertex.getOutgoing()[0].incCov(vertex.getOutgoing()[1].intCov());
    }
}


template<typename htype>
void extractLinearDisjointigs(SparseDBG<htype> &sdbg, ParallelRecordCollector<Sequence> &res, logging::Logger & logger, size_t threads) {
//    TODO support sorted edge list at all times since we compare them during construction anyway
    std::function<void(std::pair<const htype, Vertex<htype>> &)> prepare_task =
            [&sdbg, &res](std::pair<const htype, Vertex<htype>> & pair) {
                if(pair.second.isJunction()) {
                    prepareVertex(pair.second);
                    prepareVertex(pair.second.rc());
                }
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, prepare_task);
    std::function<void(std::pair<const htype, Vertex<htype>> &)> task =
            [&sdbg, &res](std::pair<const htype, Vertex<htype>> & pair) {
                htype hash = pair.first;
                Vertex<htype> &rec = pair.second;
                if(rec.isJunction()) {
                    processVertex(rec, res);
                    processVertex(rec.rc(), res);
                    if (rec.inDeg() != 1 && rec.outDeg() != 1 &&  (rec.inDeg() != 0 || rec.outDeg() != 0)) {
                        Sequence disjointig  = rec.seq;
                        if (rec.inDeg() > 0) {
                            Edge<htype> &e1 = rec.rc().getOutgoing()[0];
                            disjointig = !(e1.seq.Subseq(0, e1.intCov())) + disjointig;
                        }
                        if (rec.outDeg() > 0) {
                            Edge<htype> &e2 = rec.getOutgoing()[0];
                            disjointig = disjointig + e2.seq.Subseq(0, e2.intCov());
                        }
                        res.add(disjointig.copy());
                    }
                }
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);
}

template<typename htype>
void extractCircularDisjointigs(SparseDBG<htype> &sdbg, ParallelRecordCollector<Sequence> &res, logging::Logger & logger, size_t threads) {
    std::function<void(std::pair<const htype, Vertex<htype>> &)> task =
            [&sdbg, &res](std::pair<const htype, Vertex<htype>> & pair) {
                Vertex<htype> &rec = pair.second;
                htype hash = pair.first;
                if(rec.isJunction() || rec.seq.empty())
                    return;
                Edge<htype> &edge = rec.getOutgoing()[0];
                VERIFY(edge.end() != nullptr);
                std::vector<Edge<htype> *> path = rec.walkForward(edge);
                if(path.back()->end() != &rec) {
                    std::cout << path.back()->end() << std::endl;
                    std::cout << path.back()->end()->hash() << std::endl;
                    std::cout << rec.seq << std::endl;
                }
                VERIFY(path.back()->end() == &rec);
                for(size_t i = 0; i + 1 < path.size(); i++) {
                    if(*(path[i]->end()) < rec) {
                        return;
                    }
                }
                Sequence tmp = rec.seq;
                rec.clearSequence();
                Sequence disjointig = rec.pathSeq(path);
                res.add(tmp + disjointig + disjointig);
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);
}

template<typename htype>
std::vector<Sequence> extractDisjointigs(logging::Logger & logger, SparseDBG<htype> &sdbg, size_t threads) {
    logger.info() << "Starting to extract disjointigs." << std::endl;
    ParallelRecordCollector<Sequence> res(threads);
    logger.info() << "Extracting linear disjointigs." << std::endl;
    extractLinearDisjointigs(sdbg, res, logger, threads);
    logger.info() << "Finished extracting linear disjointigs." << std::endl;
    logger.info() << "Extracting circular disjointigs." << std::endl;
    extractCircularDisjointigs(sdbg, res, logger, threads);
    logger.info() << "Finished extracting circular disjointigs." << std::endl;
    std::vector<Sequence> rres = res.collect();
    std::sort(rres.begin(), rres.end(), [] (const Sequence& lhs, const Sequence& rhs) {
        return lhs.size() > rhs.size();
    });
    logger.info() << "Finished extracting " << rres.size() << " disjointigs of total size " << total_size(rres) << std::endl;
    return rres;
}

template<typename htype>
std::vector<Sequence> constructDisjointigs(const RollingHash<htype> &hasher, size_t w, const io::Library &reads_file,
                                           const std::vector<htype128> & hash_list, size_t threads,
                                           logging::Logger & logger) {
    std::vector<Sequence> disjointigs;
    SparseDBG<htype> sdbg = constructSparseDBGFromReads(logger, reads_file, threads, hasher, hash_list, w);
//    sdbg.printStats(logger);
    sdbg.checkSeqFilled(threads, logger);

    tieTips(logger, sdbg, w, threads);
    sdbg.checkSeqFilled(threads, logger);
    logger.info() << "Statistics for sparse de Bruijn graph:" << std::endl;
    sdbg.printStats(logger);
//    std::ofstream os;
//    os.open("sdbg.fasta");
//    sdbg.printFasta(os);
//    os.close();

    disjointigs = extractDisjointigs(logger, sdbg, threads);
    return disjointigs;
}

