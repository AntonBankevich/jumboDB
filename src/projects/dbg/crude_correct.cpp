#include "crude_correct.hpp"
#include "dbg_construction.hpp"
#include "rolling_hash.hpp"
#include "sequences/seqio.hpp"
#include <iostream>
#include <unordered_set>

double avgCoverage(const SparseDBG<htype128> &dbg) {
    size_t cov = 0;
    size_t len = 0;
    for(const std::pair<const htype128, Vertex<htype128>> &it : dbg) {
        const Vertex<htype128> &v = it.second;
        for (const Edge<htype128> &edge : v.getOutgoing()) {
            if (edge.intCov() >= 2 * edge.size()) {
                cov += edge.intCov();
                len += edge.size();
            }
        }
        for (const Edge<htype128> &edge : v.rc().getOutgoing()) {
            if (edge.intCov() >= 2 * edge.size()) {
                cov += edge.intCov();
                len += edge.size();
            }
        }
    }
    double avg_cov = double(cov) / len;
    return avg_cov;
}

std::unordered_set<const Edge<htype128> *> filterEdges(const SparseDBG<htype128> &dbg, size_t threshold, double avg_cov) {
    std::unordered_set<const Edge<htype128> *> to_skip;
    for(auto & it : dbg) {
        for (const Vertex<htype128> *pv : {&it.second, &it.second.rc()}) {
            for (const Edge<htype128> &edge : pv->getOutgoing())
                if (edge.getCoverage() < threshold) {
                    to_skip.emplace(&edge);
                    to_skip.emplace(&pv->rcEdge(edge));
                    std::cout << "Skip " << edge.end()->hash() << edge.end()->isCanonical() << " " << edge.size() << " "
                              << edge.getCoverage() << std::endl;
                }
        }
    }
    for(auto & it : dbg) {
        for(const Vertex<htype128> *pv : {&it.second, &it.second.rc()}) {
            size_t indeg = 0;
            for(const Edge<htype128> &edge : pv->rc().getOutgoing())
                if (to_skip.find(&edge) == to_skip.end()){
                    indeg += 1;
                }
            if(indeg != 0)
                continue;
            size_t cov = 0;
            size_t len = 0;
            std::vector<const Edge<htype128> *> tip;
            std::vector<const Edge<htype128> *> rc_tip;
            while(true) {
                size_t outdeg = 0;
                indeg = 0;
                const Edge<htype128> *next = nullptr;
                for(const Edge<htype128> &edge : pv->getOutgoing())
                    if (to_skip.find(&edge) == to_skip.end()){
                        outdeg += 1;
                        next = &edge;
                    }
                for(const Edge<htype128> &edge : pv->rc().getOutgoing())
                    if (to_skip.find(&edge) == to_skip.end()){
                        indeg += 1;
                    }
                if(outdeg == 1 && indeg <= 1) {
                    tip.push_back(next);
                    rc_tip.push_back(&pv->rcEdge(*next));
                    cov += next->intCov();
                    len += next->size();
                    pv = next->end();
                } else {
                    if(indeg > 1 && double(cov) / len < avg_cov / 2 && len < pv->seq.size() * 2) {
                        for(size_t i = 0; i < tip.size(); i++) {
                            to_skip.emplace(tip[i]);
                            to_skip.emplace(rc_tip[i]);
                            std::cout << "Skip " << tip[i]->end()->hash() << tip[i]->end()->isCanonical() << " " <<tip[i]->size() << " " << tip[i]->getCoverage() << std::endl;
                        }
                    }
                    break;
                }
            }
        }
    }
    return std::move(to_skip);
}

SparseDBG<htype128> simplifyGraph(logging::Logger &logger, SparseDBG<htype128> &dbg,
                                  std::unordered_set<const Edge<htype128> *> &to_skip, size_t threads) {
    std::vector<Sequence> edges;
    std::vector<htype128> vertices_again;
    for(auto & it : dbg) {
        Vertex<htype128> &vert = it.second;
        bool add = false;
        for(Vertex<htype128> *pv : {&it.second, &it.second.rc()}) {
            for(Edge<htype128> & edge : pv->getOutgoing()) {
                if (to_skip.find(&edge) == to_skip.end()) {
                    edges.push_back(pv->seq + edge.seq);
                    add = true;
                }
            }
        }
        if (add)
            vertices_again.push_back(vert.hash());
    }
    SparseDBG<htype128> simp_dbg(vertices_again.begin(), vertices_again.end(), dbg.hasher());
    simp_dbg.fillSparseDBGEdges(edges.begin(), edges.end(), logger, threads, 0);
    for(auto & it : simp_dbg) {
        Vertex<htype128> &vert = it.second;
        Vertex<htype128> &other = dbg.getVertex(vert.hash());
        bool add = false;
        for(Edge<htype128> & edge : vert.getOutgoing()) {
            edge.incCov(other.getOutgoing(edge.seq[0]).intCov());
        }
        for(Edge<htype128> & edge : vert.rc().getOutgoing()) {
            edge.incCov(other.rc().getOutgoing(edge.seq[0]).intCov());
        }
    }
    simp_dbg.checkConsistency(threads, logger);
    simp_dbg.checkSeqFilled(threads, logger);
    mergeAll(logger, simp_dbg, threads);
    return std::move(simp_dbg);
}

void handleSimpleBulge(Vertex<htype128> &v, std::unordered_map<const Edge<htype128> *, Sequence> &edge_map, double avg_cov,
                  logging::Logger &logger) {
    Edge<htype128> &edge1 = v.getOutgoing()[0].getCoverage() < v.getOutgoing()[1].getCoverage() ? v.getOutgoing()[0] : v.getOutgoing()[1];
    Edge<htype128> &edge2 = v.getOutgoing()[0].getCoverage() < v.getOutgoing()[1].getCoverage() ? v.getOutgoing()[1] : v.getOutgoing()[0];
    logger.info() << "Handling simple bulge " << v.hash() << v.isCanonical() << " " << edge1.end()->hash() << edge1.end()->isCanonical() << std::endl;
    Edge<htype128> &edge0 = v.rc().rcEdge(v.rc().getOutgoing()[0]);
    bool c1 = edge0.getCoverage() < avg_cov * 3 / 2 && edge1.getCoverage() + edge2.getCoverage() < avg_cov * 3 / 2;
    bool c2 = edge0.size() > 20000 && edge0.getCoverage() < avg_cov * 3 / 2;
    if (c1 || c2) {
        Sequence new_seq(v.seq + edge2.seq);
        edge_map[&edge1] = new_seq;
        edge_map[&v.rcEdge(edge1)] = !new_seq;
        Sequence seq1 = v.seq + edge1.seq;
        Sequence seq2 = v.seq + edge2.seq;
        size_t l = std::max<size_t>(0, std::min(seq1.size(), seq2.size()) / 2 - 50);
        logger  << "New bulge mapping " << edge1.getCoverage() << " " << edge2.getCoverage()
                << "\n" << seq1.Subseq(l, seq1.size() - l) << "\n"<< seq2.Subseq(l, seq2.size() - l) << std::endl;
    }
}

std::vector<Edge<htype128> *> checkPath(Vertex<htype128> &start, const Sequence& seq) {
    size_t cur = start.seq.size();
    Vertex<htype128> *vp = &start;
    std::vector<Edge<htype128> *> res;
    while(cur < seq.size()) {
        if(!vp->hasOutgoing(seq[cur]))
            return {};
        Edge<htype128> &next = vp->getOutgoing(seq[cur]);
        if (!seq.Subseq(cur, seq.size()).startsWith(next.seq))
            return {};
        cur += next.size();
        res.push_back(&next);
        vp = next.end();
    }
    return res;
}

std::vector<Edge<htype128> *> choosePath(Vertex<htype128> &start, Edge<htype128> &bulge, Edge<htype128> &alternative) {
    std::vector<Edge<htype128> *> res;
    Sequence alt_seq =start.seq + alternative.seq;
    std::cout << bulge.size() << " " << alternative.size() << " " << alt_seq.size() << std::endl;
    for(size_t sz = bulge.size() + start.seq.size() - 20; sz <= bulge.size() + start.seq.size() + 20; sz++) {
        if(alt_seq.size() + start.seq.size() < sz)
            continue;
        size_t over_size = alt_seq.size() + start.seq.size() - sz;
        if(alt_seq.Subseq(alt_seq.size() - over_size, alt_seq.size()) == bulge.end()->seq.Subseq(0, over_size)) {
            Sequence seq = alt_seq.Subseq(0, sz - start.seq.size()) + bulge.end()->seq;
            std::vector<Edge<htype128> *> path = checkPath(start, seq);
            if(!path.empty()) {
                if (res.empty()) {
                    res = std::move(path);
                } else
                    return {};
            }
        }
    }
    return res;
}

void handleComplexBulge(Vertex<htype128> &v, std::unordered_set<const Edge<htype128> *> &to_skip, std::unordered_map<const Edge<htype128> *, Sequence> &edge_map, double avg_cov,
                       logging::Logger &logger) {
    logger.info() << "Handling complex bulge " << v.hash() << v.isCanonical() << std::endl;
    Edge<htype128> &edge1 = v.getOutgoing()[0].getCoverage() < v.getOutgoing()[1].getCoverage() ? v.getOutgoing()[0] : v.getOutgoing()[1];
    Edge<htype128> &edge2 = v.getOutgoing()[0].getCoverage() < v.getOutgoing()[1].getCoverage() ? v.getOutgoing()[1] : v.getOutgoing()[0];
    Edge<htype128> &edge0 = v.rc().rcEdge(v.rc().getOutgoing()[0]);
    bool c1 = edge0.getCoverage() < avg_cov * 3 / 2 && edge1.getCoverage() + edge2.getCoverage() < avg_cov * 3 / 2;
    bool c2 = edge0.size() > 20000 && edge0.getCoverage() < avg_cov * 3 / 2;
    if (c1 || c2) {
        Edge<htype128> *bulge = nullptr;
        Edge<htype128> *alt = nullptr;
        if (edge1.size() > v.seq.size() - 500 && edge1.size() < v.seq.size() + 500) {
            bulge = &edge1;
            alt = &edge2;
        }
        if (edge2.size() > v.seq.size() - 500 && edge2.size() < v.seq.size() + 500) {
            bulge = &edge2;
            alt = &edge1;
        }
        if (bulge == nullptr) {
            logger.info() << "Finished handling complex bulge " << std::endl;
            return;
        }
        std::vector<Edge<htype128> *> path = choosePath(v, *bulge, *alt);
        if(path.empty()) {
            logger.info() << "Finished handling complex bulge " << std::endl;
            return;
        }
        if(bulge->getCoverage() < avg_cov / 2) {
            SequenceBuilder sb;
            sb.append(v.seq);
            for(Edge<htype128> * edge : path) {
                sb.append(edge->seq);
            }
            Sequence bseq = sb.BuildSequence();
            edge_map[bulge] = bseq;
            edge_map[&v.rcEdge(*bulge)] = !bseq;
            logger  << "New complex bulge mapped "  << v.hash() << v.isCanonical() << " " << bulge->getCoverage();
            for(Edge<htype128> * edge : path) {
                logger << " " << edge->end()->hash() << edge->end()->isCanonical() << "(" << edge->getCoverage() << ")";
            }
            logger << std::endl;
        } else {
            logger  << "New complex bulge removed "  << v.hash() << v.isCanonical() << " " << bulge->getCoverage();
            for(size_t i = 0; i < path.size(); i++) {
                Edge<htype128> * edge = path[i];
                if(edge->getCoverage() < avg_cov / 2) {
                    logger << " " << edge->end()->hash() << edge->end()->isCanonical() << "(" << edge->getCoverage() << ")";
                    to_skip.emplace(edge);
                    if (i == 0)
                        to_skip.emplace(&v.rcEdge(*edge));
                    else
                        to_skip.emplace(&(*path[i - 1]->end()).rcEdge(*edge));
                }
                logger << std::endl;
            }
        }
    }
    logger.info() << "Finished handling complex bulge " << std::endl;
}


bool checkComponent(Vertex<htype128> *start, Vertex<htype128> *end, const std::vector<Edge<htype128> *> &path,
                    std::unordered_set<const Edge<htype128> *> &to_skip) {
    std::vector<Vertex<htype128> *> queue;
    std::unordered_set<Vertex<htype128> *> good;
    std::unordered_set<Vertex<htype128> *> visited;
    queue.push_back(start);
    good.emplace(start);
    while(!queue.empty()) {
        Vertex<htype128> *next = queue.back();
        visited.emplace(next);
        queue.pop_back();
        bool ok = true;
        for(Edge<htype128> & edge : next->rc().getOutgoing()) {
            if (good.find(&edge.end()->rc()) == good.end())
                ok = false;
        }
        if(ok || next == start) {
            good.emplace(next);
            if(next == end) {
                break;
            }
            for(Edge<htype128> & edge : next->getOutgoing()) {
                if(edge.size() > start->seq.size() * 2)
                    return false;
                queue.push_back(edge.end());
            }
        }
    }

    if(good.find(end) == good.end())
        return false;
    for(Vertex<htype128> *v : visited) {
        if(good.find(v) == good.end())
            return false;
    }
    std::unordered_set<Edge<htype128> *> path_set;
    for(Edge<htype128> *ep : path) {
        path_set.emplace(ep);
    }
    for(Vertex<htype128> *v : good) {
        if(v == end)
            continue;
        for(Edge<htype128> & edge : v->getOutgoing()) {
            if(path_set.find(&edge) == path_set.end()) {
                to_skip.emplace(&edge);
                to_skip.emplace(&v->rcEdge(edge));
            }
        }
    }
    return true;
}

bool handleSubcomponent(Vertex<htype128> &v, double avg_cov, std::unordered_map<const Edge<htype128> *, Sequence> &edge_map,
                        std::unordered_set<const Edge<htype128> *> &to_skip, logging::Logger &logger) {
    logger.info() << "Handling component " << v.hash() << v.isCanonical() << std::endl;
    std::vector<Edge<htype128> *> path;
    std::vector<Edge<htype128> *> rc_path;
    size_t plen = 0;
    Vertex<htype128> * cv = &v;
    while (true) {
        Edge<htype128> *next = nullptr;
        for(Edge<htype128> &cand : cv->getOutgoing()) {
            if (next == nullptr || cand.getCoverage() > next->getCoverage()) {
                next = &cand;
            }
        }
        if(next == nullptr) {
            logger.info() << "Finished handling component " << std::endl;
            return false;
        }
        path.push_back(next);
        plen += next->size();
        rc_path.push_back(&cv->rcEdge(*next));
        cv = next->end();
        if(plen >= v.seq.size() * 2) {
            logger.info() << "Finished handling component " << std::endl;
            return false;
        }
        if (checkComponent(&v, path.back()->end(), path, to_skip)) {
            for(size_t i = 0; i < path.size(); i++) {
                Edge<htype128> *edge = path[i];
                Edge<htype128> *rcedge = rc_path[i];
                edge_map[edge] = rc_path[i]->end()->rc().seq + edge->seq;
                edge_map[rcedge] = !(rc_path[i]->end()->rc().seq + edge->seq);
            }
            logger.info() << "Removed new subcomponent. Retained path " << v.hash() << v.isCanonical();
            for(Edge<htype128> * edge : path) {
                logger << " " << edge->end()->hash() << edge->end()->isCanonical() << "(" << edge->getCoverage() << ")";
            }
            logger << std::endl;
            logger.info() << "Finished handling component " << std::endl;
            return true;
        }
    }
    logger.info() << "Finished handling component " << std::endl;
    return  false;
}

std::experimental::filesystem::path CrudeCorrect(logging::Logger &logger, SparseDBG <htype128> &dbg,
                                                 const std::experimental::filesystem::path &dir,
                                                 const size_t w, const io::Library &reads_lib, size_t threads, size_t threshold) {
    const RollingHash<htype128> &hasher = dbg.hasher();
    logger.info() << "Crude error correction based of removing low covered edges and heterozygous bulges" << std::endl;
    logger.info() << "Removing tips and low covered edges" << std::endl;
    double avg_cov = avgCoverage(dbg);
    logger.info() << "Estimated average coverage as " << avg_cov << std::endl;
    std::unordered_set<const Edge<htype128> *> to_skip = filterEdges(dbg, threshold, avg_cov);
    SparseDBG<htype128> simp_dbg(simplifyGraph(logger, dbg, to_skip, threads));
    logger.info() << "Printing initial simplified graph to fasta" << std::endl;
    std::ofstream simp_os;
    simp_os.open(dir / "simp_graph.fasta");
    simp_dbg.printFasta(simp_os);
    simp_os.close();
    logger.info() << "Printing initial simplified graph to dot" << std::endl;
    std::ofstream dot;
    dot.open(dir / "simp_graph.dot");
    simp_dbg.printDot(dot, true);
    dot.close();

    simp_dbg.fillAnchors(w, logger, threads);

    std::unordered_map<const Edge<htype128> *, Sequence> edge_map;
    for(std::pair<const htype128, Vertex<htype128>> &it : simp_dbg) {
        for (Vertex<htype128> *vit : {&it.second, &it.second.rc()}) {
            Vertex<htype128> &v = *vit;
            if (v.inDeg() == 1 && v.outDeg() == 2 && v.rc().getOutgoing()[0].getCoverage() < avg_cov * 3 / 2 &&
                    edge_map.find(&v.getOutgoing()[0]) == edge_map.end() &&
                    edge_map.find(&v.getOutgoing()[1]) == edge_map.end()) {
                if (v.getOutgoing()[0].end() == v.getOutgoing()[1].end()) {
                    handleSimpleBulge(v, edge_map, avg_cov, logger);
                }
                else {
                    bool is_subcomponent = handleSubcomponent(v, avg_cov, edge_map, to_skip, logger);
                    if(!is_subcomponent) {
                        handleComplexBulge(v, to_skip, edge_map, avg_cov, logger);
                    }
                }
            }
        }
    }

    logger.info() << "Finished graph correction. Correcting reads." << std::endl;
    ParallelRecordCollector<Contig> alignment_results(threads);

    std::function<void(StringContig &)> task = [&simp_dbg, &dbg, &to_skip, &alignment_results, &hasher, &threshold, w, &edge_map](StringContig & contig) {
        Contig read = contig.makeContig();
        if(read.size() < w + hasher.k - 1)
            return;
        GraphAlignment<htype128> old_al = dbg.align(read.seq);
        while (old_al.size() > 0 && to_skip.find(&old_al.front().contig()) != to_skip.end()) {
            old_al = old_al.subalignment(1, old_al.size());
        }
        while (old_al.size() > 0 && to_skip.find(&old_al.back().contig()) != to_skip.end()) {
            old_al = old_al.subalignment(0, old_al.size() - 1);
        }
        if (old_al.size() == 0)
            return;
        for(Segment<Edge<htype128>> seg : old_al) {
            if (to_skip.find(&seg.contig()) != to_skip.end()) {
                return;
            }
        }
        Sequence seq = old_al.Seq();
        if(seq.size() < w + hasher.k - 1)
            return;
        GraphAlignment<htype128> gal = simp_dbg.align(seq);
        for(Segment<Edge<htype128>> seg : gal) {
            if (to_skip.find(&seg.contig()) != to_skip.end()) {
                return;
            }
        }
        if (gal.size() > 0) {
            Sequence res = gal.map(edge_map);
            alignment_results.emplace_back(res, read.id);
        }
    };
    io::SeqReader reader(reads_lib);
    processRecords(reader.begin(), reader.end(), logger, threads, task);
    std::experimental::filesystem::path corrected_reads = dir / "corrected.fasta";
    std::ofstream os(corrected_reads);
    for(Contig & rec : alignment_results) {
        os << ">" << rec.id << "\n" << rec.seq << "\n";
    }
    os.close();
    return corrected_reads;
}
