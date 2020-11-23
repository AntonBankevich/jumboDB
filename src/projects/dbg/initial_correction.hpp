#pragma once
#include "compact_path.hpp"

size_t edit_distance(Sequence s1, Sequence s2) {
    size_t left_skip = 0;
    while(left_skip < s1.size() && left_skip < s2.size() && s1[left_skip] == s2[left_skip]) {
        left_skip++;
    }
    s1 = s1.Subseq(left_skip, s1.size());
    s2 = s2.Subseq(left_skip, s2.size());
    size_t right_skip = 0;
    while(right_skip < s1.size() && right_skip < s2.size() && s1[s1.size() - 1 - right_skip] == s2[s2.size() - 1 - right_skip]) {
        right_skip++;
    }
    s1 = s1.Subseq(0, s1.size() - right_skip);
    s2 = s2.Subseq(0, s2.size() - right_skip);
    std::vector<std::vector<size_t>> d(s1.size() + 1, std::vector<size_t>(s2.size() + 1));
    d[0][0] = 0;
    for(unsigned int i = 1; i <= s1.size(); ++i) d[i][0] = i;
    for(unsigned int i = 1; i <= s2.size(); ++i) d[0][i] = i;

    for(unsigned int i = 1; i <= s1.size(); ++i)
        for(unsigned int j = 1; j <= s2.size(); ++j)
            d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
    return d[s1.size()][s2.size()];
}

size_t bestPrefix(const Sequence &s1, const Sequence &s2) {
    std::vector<std::vector<size_t>> d(s1.size() + 1, std::vector<size_t>(s2.size() + 1));
    d[0][0] = 0;
    for(unsigned int i = 1; i <= s1.size(); ++i) d[i][0] = i;
    for(unsigned int i = 1; i <= s2.size(); ++i) d[0][i] = i;

    for(unsigned int i = 1; i <= s1.size(); ++i)
        for(unsigned int j = 1; j <= s2.size(); ++j)
            d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
    size_t res = s2.size();
    for(size_t j = 0; j < s2.size(); j++)
        if(d[s1.size()][j] < d[s1.size()][res])
            res = j;
    return res;
}

size_t tournament(const Sequence &bulge, const std::vector<Sequence> &candidates) {
    size_t winner = 0;
    std::vector<size_t> dists;
    for(size_t i = 0; i < candidates.size(); i++) {
        dists.push_back(edit_distance(bulge, candidates[i]));
        if (dists.back() < dists[winner])
            winner = i;
    }
    size_t max_dist = std::max<size_t>(20, bulge.size() / 100);
    if(dists[winner] > max_dist)
        return -1;
    for(size_t i = 0; i < candidates.size(); i++) {
        if(i != winner) {
            size_t diff = edit_distance(candidates[winner], candidates[i]);
            VERIFY(dists[winner] <= dists[i] + diff);
            VERIFY(dists[i] <= dists[winner] + diff);
            if(dists[i] < max_dist && dists[i] != dists[winner] + diff)
                return -1;
        }
    }
    return winner;
}

template <typename htype>
std::vector<Path<htype>> FindBulgeAlternatives(const Path<htype> &path, size_t max_diff) {
    size_t k = path.start().seq.size();
    std::vector<GraphAlignment<htype>> als = GraphAlignment<htype>(path.start()).allExtensions(max_diff);
    max_diff = std::min(max_diff, path.len());
    std::vector<Path<htype>> res;
    Sequence path_seq = path.truncSeq();
    for(GraphAlignment<htype> &diff_al : als) {
        size_t path_pos = 0;
        size_t edge_pos = size_t (-1);
        for(size_t i = 0; i < max_diff; i++) {
            GraphAlignment<htype> al = diff_al;
            if(i > 0 && al.size() > 0 && al.lastNucl() == path[path_pos].seq[edge_pos])
                continue;
            Sequence seq = path_seq.Subseq(i, path_seq.size());
            al.extend(seq);
            if(al.valid() && al.endClosed() && al.back().contig().end() == &path.finish()){
                res.emplace_back(al.path());
            }
            edge_pos += 1;
            if(edge_pos == path[path_pos].size()) {
                path_pos += 1;
                edge_pos = 0;
            }
        }
    }
    return res;
}

template <typename htype>
std::vector<GraphAlignment<htype>> FilterAlternatives(logging::Logger &logger1, const GraphAlignment<htype> &initial, std::vector<GraphAlignment<htype>> &als,
                                            size_t max_diff, double threshold) {
    size_t len = initial.len();
    std::vector<GraphAlignment<htype>> res;
    size_t k = initial.getVertex(0).seq.size();
    for(GraphAlignment<htype> &al : als) {
        CompactPath<htype> cpath(al);
        bool ok = true;
        for(size_t i = 0; i < al.size(); i++) {
            if(al[i].contig().getCoverage() < threshold) {
//                logger << "Filtered " << cpath << " based on coverage " << i << " " << al[i].contig().getCoverage() << std::endl;
                ok = false;
                break;
            }
        }
        if(!ok) {
            continue;
        }
        size_t al_len = al.len();
        if(len > al_len + max_diff || al_len > len + max_diff) {
//            logger << "Filtered " << cpath << " based on length " << al_len << " " << len << std::endl;
//            logger << al.truncSeq() << std::endl;
//            logger << initial.truncSeq() << std::endl;
            continue;
        }
        res.emplace_back(al);
    }
    return res;
}

template<typename htype>
GraphAlignment<htype> processBulge(logging::Logger &logger1, std::ostream &out, const GraphAlignment<htype> &bulge,
                         const RecordStorage<htype> &reads_storage, double threshold) {
    size_t size = bulge.len();
    out << size << " bulge " << bulge.size() << " " << bulge.minCoverage();
//    logger  << "New bulge from "  << bulge.start().hash() << bulge.start().isCanonical() << " "
//            << bulge.start().outDeg() << " " << bulge.start().inDeg() << std::endl
//            << "To " << bulge.finish().hash() << bulge.finish().isCanonical() << " "
//            << bulge.finish().outDeg() << " " << bulge.finish().inDeg() << std::endl;
//    std::vector<Path<htype>> alts = FindBulgeAlternatives(bulge, 10);
//    logger << "Filtering graph alternatives" << std::endl;
//    std::vector<Path<htype>> filtered = FilterAlternatives(logger, bulge, alts, std::max<size_t>(10, bulge.len() / 200), threshold);
//    out << " " << filtered.size();
    std::vector<GraphAlignment<htype>> read_alternatives = reads_storage.getRecord(bulge.start()).getBulgeAlternatives(bulge.finish(), threshold);
//    logger << "Filtering read alternatives" << std::endl;
    std::vector<GraphAlignment<htype>> read_alternatives_filtered = FilterAlternatives(logger1, bulge, read_alternatives,
                                                                                       std::max<size_t>(100, bulge.len() / 100), threshold);
//    logger << read_alternatives.size() << " " << read_alternatives_filtered.size() << std::endl;
//    logger << reads_storage.getRecord(bulge.start());
//    logger << "Read alternatives" << std::endl;
//    for(GraphAlignment<htype> &candidate : read_alternatives) {
//        logger << CompactPath<htype>(candidate) << std::endl;
//    }
//    logger << "Result " << read_alternatives_filtered.size() << std::endl;
    out << " " << read_alternatives.size() << " " << read_alternatives_filtered.size();
    if(read_alternatives_filtered.size() > 1) {
//        logger << "Multiple choice in bulge " << read_alternatives_filtered.size() << std::endl << bulge.truncSeq() << std::endl;
        Sequence old = bulge.truncSeq();
        std::vector<Sequence> candidates;
        for(GraphAlignment<htype> &cand : read_alternatives_filtered) {
//            logger << cand.truncSeq() << std::endl;
            candidates.push_back(cand.truncSeq());
        }
        size_t winner = tournament(old, candidates);
//        logger << "Winner found " << winner << std::endl;
        if(winner != size_t(-1)) {
            read_alternatives_filtered = {read_alternatives_filtered[winner]};
        }
    }
    out << " " << read_alternatives_filtered.size() << std::endl;
    if(read_alternatives_filtered.size() == 1) {
        return std::move(read_alternatives_filtered[0]);
    } else {
        return std::move(bulge);
    }
}

template<typename htype>
GraphAlignment<htype> processTip(logging::Logger &logger1, std::ostream &out, const GraphAlignment<htype> &tip,
                         const RecordStorage<htype> &reads_storage, double threshold) {
    size_t size = tip.len();
    out << size << " tip " << tip.size() << " " << tip.minCoverage();
//    logger  << "New tip from "  << tip.start().hash() << tip.start().isCanonical() << " "
//            << tip.start().outDeg() << " " << tip.start().inDeg() << std::endl
//            << "To " << tip.finish().hash() << tip.finish().isCanonical() << " "
//            << tip.finish().outDeg() << " " << tip.finish().inDeg() << std::endl;
    std::vector<GraphAlignment<htype>> read_alternatives = reads_storage.getRecord(tip.start()).getTipAlternatives(tip.len(), threshold);
//    logger << "Filtering read alternatives" << std::endl;
    std::vector<GraphAlignment<htype>> read_alternatives_filtered =
            FilterAlternatives(logger1, tip, read_alternatives, size_t(-1) / 2, threshold);
//    logger << read_alternatives.size() << " " << read_alternatives_filtered.size() << std::endl;
//    logger << reads_storage.getRecord(tip.start());
//    logger << "Read alternatives" << std::endl;
//    for(GraphAlignment<htype> &candidate : read_alternatives) {
//        logger << CompactPath<htype>(candidate) << std::endl;
//    }
    out << " " << read_alternatives.size() << " " << read_alternatives_filtered.size();
    if(read_alternatives_filtered.size() > 1) {
//        logger << "Multiple choice for tip " << read_alternatives_filtered.size() << std::endl << tip.truncSeq() << std::endl;
        Sequence old = tip.truncSeq();
        std::vector<Sequence> candidates;
        for(GraphAlignment<htype> &cand : read_alternatives_filtered) {
//            logger << cand.truncSeq() << std::endl;
            Sequence candSeq = cand.truncSeq();
            Sequence prefix = candSeq.Subseq(0, bestPrefix(old, candSeq));
//            logger << prefix << std::endl;
            candidates.push_back(prefix);
            cand.cutBack(cand.len() - prefix.size());
        }
        size_t winner = tournament(old, candidates);
//        logger << "Winner found " << winner << std::endl;
        if(winner != size_t(-1)) {
            read_alternatives_filtered = {read_alternatives_filtered[winner]};
        }
    }
//    logger << "Result " << read_alternatives_filtered.size() << std::endl;
    out << " " << read_alternatives_filtered.size();
    out << std::endl;
    if(read_alternatives_filtered.size() == 1) {
        return std::move(read_alternatives_filtered[0]);
    } else {
        return std::move(tip);
    }
}

template<typename htype>
size_t correctLowCoveredRegions(logging::Logger &logger, RecordStorage<htype> &reads_storage,
                              const std::experimental::filesystem::path &out_file,
                              double threshold, size_t k, size_t threads) {
    ParallelRecordCollector<std::string> results(threads);
    logger.info() << "Correcting low covered regions in reads" << std::endl;
#pragma omp parallel for default(none) shared(reads_storage, results, threshold, k, logger)
    for(size_t read_ind = 0; read_ind < reads_storage.size(); read_ind++) {
        std::stringstream ss;
        AlignedRead<htype> &alignedRead = reads_storage[read_ind];
        CompactPath<htype> &initial_cpath = alignedRead.path;
        GraphAlignment<htype> path = initial_cpath.getAlignment();
        GraphAlignment<htype> corrected_path(path.start());
        bool corrected = false;
        for(size_t path_pos = 0; path_pos < path.size(); path_pos++) {
            VERIFY(corrected_path.finish() == path.getVertex(path_pos));
            Edge<htype> &edge = path[path_pos].contig();
            if (edge.getCoverage() >= threshold || edge.size() > 5 * k) {
//                if(edge.size() > 5 * k) {
//                    logger << "Very long read path segment with low coverage. Skipping." << std::endl;
//                }
                corrected_path.push_back(path[path_pos]);
                continue;
            }
            size_t step_back = 0;
            size_t step_front = 0;
            size_t size = edge.size();
            while(step_back < corrected_path.size() && corrected_path[corrected_path.size() - step_back - 1].contig().getCoverage() < 10) {
                size += corrected_path[corrected_path.size() - step_back - 1].size();
                step_back += 1;
            }
            while(step_front + path_pos + 1 < path.size() && path[step_front + path_pos + 1].contig().getCoverage() < 10) {
                size += path[step_front + path_pos + 1].size();
                step_front += 1;
            }
            Vertex<htype> &start = corrected_path.getVertex(corrected_path.size() - step_back);
            Vertex<htype> &end = path.getVertex(path_pos + 1 + step_front);
            GraphAlignment<htype> badPath = corrected_path.subPath(corrected_path.size() - step_back, corrected_path.size())
                                            + path.subPath(path_pos, path_pos + 1 + step_front);
            corrected_path.pop_back(step_back);
//            logger << "Bad read segment "  << alignedRead.id << " " << path_pos << " " << step_back << " " << step_front << " " << path.size()
//                   << " " << size << " " << edge.getCoverage() << " size " << step_back + step_front + 1 << std::endl;
//            if(step_back < corrected_path.size()) {
//                logger  << "Start stop " << step_back << " " << corrected_path.getVertex(corrected_path.size() - step_back).hash() << " "
//                        << corrected_path.getVertex(corrected_path.size() - step_back).isCanonical()
//                        << " " << corrected_path[corrected_path.size() - step_back - 1].contig().getCoverage()
//                        << corrected_path.getVertex(corrected_path.size() - step_back - 1).rcEdge(corrected_path[corrected_path.size() - step_back - 1].contig()).getCoverage() << std::endl;
//                Vertex<htype> &tmpv = corrected_path.getVertex(corrected_path.size() - step_back);
//                logger << tmpv.outDeg() << " " << tmpv.inDeg() << std::endl;
//                for(Edge<htype> &e : tmpv.getOutgoing())
//                    logger << "Edge out " << e.size() << " " << e.getCoverage() << std::endl;
//                for(Edge<htype> &e : tmpv.rc().getOutgoing())
//                    logger << "Edge in " << e.size() << " " << e.getCoverage() << std::endl;
//            }
//            if(path_pos + 1 + step_front < path.size()) {
//                logger  << "End stop " << step_front << " " << path.getVertex(path_pos + 1 + step_front).hash()
//                        << " " << path.getVertex(path_pos + 1 + step_front).isCanonical()
//                        << " " << path[path_pos + step_front].contig().getCoverage() << std::endl;
//                Vertex<htype> &tmpv = path.getVertex(path_pos + step_front + 1);
//                logger << tmpv.outDeg() << " " << tmpv.inDeg() << std::endl;
//                for(Edge<htype> &e : tmpv.getOutgoing())
//                    logger << "Edge out " << e.size() << " " << e.getCoverage() << std::endl;
//                for(Edge<htype> &e : tmpv.rc().getOutgoing())
//                    logger << "Edge in " << e.size() << " " << e.getCoverage() << std::endl;
//            }
            if(step_back == corrected_path.size() && step_front == path.size() - path_pos - 1) {
//                logger << "Whole read has low coverage. Skipping." << std::endl;
                for(const Segment<Edge<htype>> &seg : badPath) {
                    corrected_path.push_back(seg);
                }
            } else if(size > 5 * k) {
//                logger << "Very long read path segment with low coverage. Skipping." << std::endl;
                for(const Segment<Edge<htype>> &seg : badPath) {
                    corrected_path.push_back(seg);
                }
            } else if(step_back == corrected_path.size()) {
//                logger << "Processing incoming tip" << std::endl;
                GraphAlignment<htype> rcBadPath = badPath.RC();
                GraphAlignment<htype> substitution = processTip(logger, ss, rcBadPath, reads_storage, threshold);
                GraphAlignment<htype> rcSubstitution = substitution.RC();
                corrected_path = std::move(rcSubstitution);
            } else if(step_front == path.size() - path_pos - 1) {
//                logger << "Processing outgoing tip" << std::endl;
                GraphAlignment<htype> substitution = processTip(logger, ss, badPath, reads_storage, threshold);
                for(const Segment<Edge<htype>> &seg : substitution) {
                    corrected_path.push_back(seg);
                }
            } else {
                GraphAlignment<htype> substitution = processBulge(logger, ss, badPath, reads_storage, threshold);
                for (const Segment<Edge<htype>> &seg : substitution) {
                    corrected_path.push_back(seg);
                }
            }
            path_pos = path_pos + step_front;
        }
        if(path != corrected_path) {
//            logger << "Corrected read " << alignedRead.id << std::endl;
            reads_storage.reroute(alignedRead, path, corrected_path);
        }
        results.emplace_back(ss.str());
    }
    std::ofstream out;
    out.open(out_file);
    size_t res = 0;
    for(std::string &s : results) {
        if(!s.empty() && s.find(" 1\n") != size_t(-1))
            res += 1;
        out << s;
    }
    out.close();
    return res;
}

template<typename htype>
size_t correctAT(logging::Logger &logger, RecordStorage<htype> &reads_storage, size_t k, size_t threads) {
    ParallelRecordCollector<std::string> results(threads);
    logger.info() << "Correcting dinucleotide errors in reads" << std::endl;
    ParallelCounter cnt(threads);
//#pragma omp parallel for default(none) shared (reads_storage, results, k, logger, cnt, std::cout)
    for(size_t read_ind = 0; read_ind < reads_storage.size(); read_ind++) {
        AlignedRead<htype> &alignedRead = reads_storage[read_ind];
        CompactPath<htype> &initial_cpath = alignedRead.path;
        GraphAlignment<htype> path = initial_cpath.getAlignment();
        size_t corrected = 0;
        for (size_t path_pos = 0; path_pos < path.size(); path_pos++) {
            if(path[path_pos].left > 0 || path[path_pos].right < path[path_pos].contig().size())
                continue;
//            std::cout << alignedRead.id << " " << path_pos << " " << path.size() << std::endl;
            Sequence seq = path.getVertex(path_pos).seq;
            size_t at_cnt1 = 2;
            while (at_cnt1 < seq.size() && seq[seq.size() - at_cnt1 - 1] == seq[seq.size() - at_cnt1 + 1])
                at_cnt1 += 1;
            if (at_cnt1 < 8 || at_cnt1 == seq.size())
                continue;
            Sequence extension = path.truncSeq(path_pos, k * 2);
            size_t at_cnt2 = 0;
            while (at_cnt2 < extension.size() && extension[at_cnt2] == seq[seq.size() - 2 + at_cnt2 % 2])
                at_cnt2 += 1;
            if(at_cnt2 >= k)
                continue;
            size_t max_variation = std::max<size_t>(3, (at_cnt1 + at_cnt2) / 50);
            if (extension.size() < k + max_variation * 2) {
                continue;
            }
            extension = extension.Subseq(0, k + max_variation * 2);
//            Sequence extension = seq.Subseq(seq.size() - 2 * max_variation, seq.size()) + path.truncSeq(path_pos, k + max_variation * 2);
            SequenceBuilder sb;
            for (size_t i = 0; i < max_variation; i++) {
                sb.append(seq.Subseq(seq.size() - 2, seq.size()));
            }
            sb.append(extension);
            Sequence longest_candidate = sb.BuildSequence();
            Sequence best_seq;
            size_t best_support = 0;
            const VertexRecord<htype> &rec = reads_storage.getRecord(path.getVertex(path_pos));
            size_t initial_support = 0;
            for (size_t i = 0; i <= std::min(2 * max_variation, max_variation + at_cnt2 / 2); i++) {
                Sequence candidate_seq = longest_candidate.Subseq(i * 2, i * 2 + k);
                GraphAlignment<htype> candidate(path.getVertex(path_pos));
                candidate.extend(candidate_seq);
                if (!candidate.valid())
                    continue;
                CompactPath<htype> ccandidate(candidate);
                size_t support = rec.countStartsWith(ccandidate.seq());
                if(i == max_variation) {
                    initial_support = support;
                    VERIFY(support > 0);
                }
                if (support > best_support) {
                    best_seq = longest_candidate.Subseq(i * 2, longest_candidate.size());
                    best_support = support;
                }
                if (candidate_seq[0] != seq[seq.size() - 2] || candidate_seq[1] != seq[seq.size() - 1])
                    break;
            }
            if (extension.startsWith(best_seq))
                continue;
            logger  << "Correcting AT " << best_support << " " << initial_support  << " " << at_cnt1 << " " << at_cnt2 << " " << max_variation << std::endl
                    << seq << std::endl
                    << extension << std::endl << best_seq << std::endl;
            VERIFY(best_support > 0);
            GraphAlignment<htype> rerouting(path.getVertex(path_pos));
            rerouting.extend(best_seq);
            GraphAlignment<htype> old_path(path.getVertex(path_pos));
            old_path.extend(extension);
            VERIFY(rerouting.valid());
            VERIFY(old_path.valid());
            VERIFY(rerouting.back() == old_path.back());
            if (rerouting.back().right != rerouting.back().contig().size()) {
                rerouting.pop_back();
                old_path.pop_back();
            }
            GraphAlignment<htype> prev_path = path;
            path = path.reroute(path_pos, path_pos + old_path.size(), rerouting);
            reads_storage.reroute(alignedRead, prev_path, path);
//            std::cout << "Rerouted " << alignedRead.id << " " << initial_support << " " << best_support << std::endl;
            if(path.size() > 10000) {
                std::cout << extension << "\n" <<best_seq << std::endl;
            }
            corrected += std::max(prev_path.len(), path.len()) - std::min(prev_path.len(), path.len());
        }
        if(corrected > 0) {
#pragma omp critical
            {
                logger << "ATAT " << alignedRead.id << " " << corrected << std::endl;
                if(corrected > 100) {
                    logger << "oppa" << std::endl;
                }
            }
            ++cnt;
        }
    }
    return cnt.get();
}

template<typename htype>
void initialCorrect(SparseDBG<htype> &sdbg, logging::Logger &logger,
                    const std::experimental::filesystem::path &out_file,
                    const std::experimental::filesystem::path &out_reads,
                    const io::Library &reads_lib,
                    const std::experimental::filesystem::path &ref,
                     double threshold, size_t threads, const size_t min_read_size) {
    size_t k = sdbg.hasher().k;
    logger.info() << "Collecting info from reads" << std::endl;
    size_t extension_size = std::max(std::min(min_read_size * 3 / 4, sdbg.hasher().k * 11 / 2), sdbg.hasher().k * 3 / 2);
    size_t min_extension = sdbg.hasher().k * 2 / 3;
    RecordStorage<htype> reads_storage(sdbg, min_extension, extension_size, true);
    io::SeqReader readReader(reads_lib);
    reads_storage.fill(readReader.begin(), readReader.end(), min_read_size, logger, threads);
//    logger.info() << "Collecting info from reference" << std::endl;
//    RecordStorage<htype> ref_storage(sdbg, min_extension, extension_size, false);
//    io::SeqReader refReader(ref);
//    ref_storage.fill(refReader.begin(), refReader.end(), min_read_size, logger, threads);
    {
        size_t correctedAT = correctAT(logger, reads_storage, k, threads);
        logger.info() << "Corrected " << correctedAT << " dinucleotide sequences" << std::endl;
    }
    size_t corrected_low = correctLowCoveredRegions(logger, reads_storage, out_file, threshold, k, threads);
    logger.info() << "Corrected low covered regions in " << corrected_low << " reads" << std::endl;
    {
        size_t correctedAT = correctAT(logger, reads_storage, k, threads);
        logger.info() << "Corrected " << correctedAT << " dinucleotide sequences" << std::endl;
    }
    logger.info() << "Printing reads to disk" << std::endl;
    std::ofstream ors;
    ors.open(out_reads);
    for(auto it = reads_storage.begin(); it != reads_storage.end(); ++it) {
        AlignedRead<htype> &alignedRead = *it;
        ors << ">" << alignedRead.id << "\n" << alignedRead.path.getAlignment().Seq() << "\n";
    }
    ors.close();
}