//
// Created by anton on 8/26/20.
//

#pragma once

#include <queue>
#include "common/simple_computation.hpp"
#include "sparse_dbg.hpp"
#include "common/output_utils.hpp"

namespace error_correction {
    inline std::ostream& operator<<(std::ostream& out, const unsigned __int128& item) {
        std::vector<char> res;
        unsigned __int128 tmp = item;
        while(tmp != 0) {
            res.push_back(char((tmp % 10) + '0'));
            tmp /= 10;
        }
        return out << std::string(res.rbegin(), res.rend());
    }
    template<typename htype>
    struct State {
        State() = default;

        State(size_t lastMatch, size_t diff, Vertex<htype> *graphPosition, bool match) : last_match(lastMatch),
                                                                                         diff(diff),
                                                                                         graph_position(graphPosition),
                                                                                         match(match) {}

        size_t last_match;
        size_t diff;
        Vertex<htype> *graph_position;
        bool match;

        bool operator==(const State<htype> &other) const {
            return last_match == other.last_match && diff == other.diff &&
                   graph_position == other.graph_position && match == other.match;
        }

        size_t hash() const {
            return std::hash<size_t>()(last_match) ^ std::hash<size_t>()(diff) ^
                   std::hash<void *>()(graph_position) ^ std::hash<bool>()(match);
        }
    };

    template<class htype>
    std::ostream& operator<<(std::ostream& out, const State<htype> & item) {
        return out << "(" << item.last_match << " " << item.diff << " " << item.graph_position->hash() << " " << item.match << ")";
    }

}

namespace std {
    template<typename htype>
    struct hash<error_correction::State<htype>> {
        std::size_t operator()(const error_correction::State<htype> &state) const {
            return std::hash<size_t>()(state.last_match) ^ std::hash<size_t>()(state.diff) ^
                   std::hash<void *>()(state.graph_position) ^ std::hash<bool>()(state.match);
        }
    };
}

namespace error_correction {

    template<typename htype>
    struct ResRecord {
        ResRecord(ResRecord<htype> *prev, Edge<htype> *lastEdge, bool goodmatch) :
                    prev(prev), last_edge(lastEdge), good_match(goodmatch) {}
        ResRecord<htype> *prev;
        Edge<htype> *last_edge;
        bool good_match;
    };

    template<typename htype>
    struct ScoredState {
        State<htype> state;
        ResRecord<htype> resRecord;
        size_t score;

        ScoredState(const State<htype> &state, const ResRecord<htype> &resRecord, size_t score) : state(state),
                                                                                                  resRecord(resRecord),
                                                                                                  score(score) {}

        bool operator<(const ScoredState<htype> &other) const {
            return score > other.score;
        }
    };

    template<typename htype>
    struct CorrectionResult {
        CorrectionResult(const Path<htype> &path, size_t score, size_t iterations) : path(path),
                                                                                                  score(score),
                                                                                                  iterations(
                                                                                                          iterations) {}

        Path<htype> path;
        size_t score;
        size_t iterations;
    };

    struct ScoreScheme {
        size_t cov_threshold = 8;
        size_t alternative_penalty = 2;
        size_t diff_penalty = 20;
        size_t max_diff = 10;
        size_t low_coverage_penalty = 10;
        size_t very_bad_coverage = 3;
        size_t very_bad_covereage_penalty = 50;
        size_t high_coverage_penalty = 50;
        size_t match_to_freeze = 100;

        template<typename htype>
        size_t scoreMatchingEdge(const Edge<htype> &edge) const {
            if (edge.getCoverage() <= very_bad_coverage) {
                return very_bad_covereage_penalty * edge.size();
            } else if (edge.getCoverage() < cov_threshold) {
                return low_coverage_penalty * edge.size();
            } else {
                return 0;
            }
        }

        template<typename htype>
        size_t scoreAlternativeEdge(const Edge<htype> &edge) const {
            if (edge.getCoverage() <= very_bad_coverage) {
                return (very_bad_covereage_penalty + alternative_penalty) * edge.size();
            } else if (edge.getCoverage() < cov_threshold) {
                return (low_coverage_penalty + alternative_penalty) * edge.size();
            } else {
                return alternative_penalty * edge.size();
            }
        }


        template<typename htype>
        size_t scoreRemovedEdge(const Edge<htype> &edge) const {
            if (edge.getCoverage() > cov_threshold) {
                return high_coverage_penalty * edge.size();
            } else {
                return 0;
            }
        }
    };

    template<typename htype>
    std::vector<Edge<htype>*> restoreResult(const std::unordered_map<State<htype>, ResRecord<htype>> stateMap,
            ResRecord<htype> *resRecord) {
        std::vector<Edge<htype>*> res;
        while(resRecord->prev != nullptr) {
            if (resRecord->last_edge != nullptr) {
                res.push_back(resRecord->last_edge);
            }
            resRecord = resRecord->prev;
        }
        return {res.rbegin(), res.rend()};
    }

    template<typename htype>
    size_t checkPerfect(const std::unordered_map<State<htype>, ResRecord<htype>> stateMap,
                                           ResRecord<htype> *resRecord, size_t length) {
        size_t len = 0;
        size_t cnt = 0;
        while(resRecord->prev != nullptr) {
            if (!resRecord->good_match) {
                return size_t(-1);
            }
            len += resRecord->last_edge->size();
            cnt += 1;
            if(len >= length)
                return cnt;
            resRecord = resRecord->prev;
        }
        return size_t(-1);
    }

    template<typename htype>
    CorrectionResult<htype> correct(Path<htype> & initial_path, ScoreScheme scores = {}) {
//        std::cout << "New read " << path.size() << std::endl;
//        for (size_t i = 0; i < path.size(); i++) {
//            std::cout << " " << path[i].end()->hash() << " " << path[i].getCoverage() << " " << path[i].size();
//        }
//        std::cout << std::endl;
        size_t from = 0;
        size_t cut_len_start = 0;
        size_t to = initial_path.size();
        while(from < to && initial_path[from].getCoverage() < scores.cov_threshold && cut_len_start < 1000) {
            cut_len_start +=initial_path[from].size();
            from += 1;
        }
        size_t cut_len_end = 0;
        while(from < to && initial_path[to - 1].getCoverage() < scores.cov_threshold && cut_len_end < 1000) {
            cut_len_end += initial_path[to - 1].size();
            to -= 1;
        }
        if (from == to || cut_len_end + cut_len_start > 1500) {
//            std::cout << "Finished fail " << (size_t(-1) >> 1u) << std::endl;
            return {initial_path, size_t(-1) >> 1u, size_t(-1) >> 1u};
        }
        Path<htype> path = initial_path.subPath(from, to);
        std::priority_queue<ScoredState<htype>> queue;
        queue.emplace(State<htype>(0, 0, &path.start(), true), ResRecord<htype>(nullptr, nullptr, false), 0);
        std::unordered_map<State<htype>, ResRecord<htype>> res;
        size_t cnt = 0;
        size_t frozen = 0;
        while(!queue.empty()) {
            ScoredState<htype> next = queue.top();
            queue.pop();
            if(next.state.last_match == path.size()) {
                VERIFY(next.state.match);
//                std::cout << "Finished " << next.score << " " << cnt << std::endl;
                return {Path<htype>(path.start(), restoreResult(res, &next.resRecord)), next.score, cnt};
            }
//            std::cout<< "New state " << next.state << " " << next.score << " " << next.state.graph_position->outDeg();
//            for(size_t i = 0; i < next.state.graph_position->outDeg(); i++) {
//                std::cout << " " << next.state.graph_position->getOutgoing()[i].getCoverage();
//            }
//            if (next.resRecord.last_edge != nullptr)
//                std::cout << " " << next.resRecord.last_edge->getCoverage() << " " << next.resRecord.last_edge->seq;
//            std::cout << std::endl;
            if(res.find(next.state) != res.end() || next.state.last_match < frozen)
                continue;
            ResRecord<htype> &prev = res.emplace(next.state, next.resRecord).first->second;
            State<htype> &state = next.state;
            if(next.state.match && !prev.good_match) {
                size_t perfect_len = 0;
                size_t right = state.last_match;
                while(right < path.size() && path[right].getCoverage() >= scores.cov_threshold) {
                    right += 1;
                }
                size_t left = right;
                while(left > state.last_match && perfect_len + path[left - 1].size() < scores.match_to_freeze) {
                    perfect_len += path[left - 1].size();
                    left -= 1;
                }
                if (left > state.last_match) {
//                    std::cout << "Skip to " << left << " " << &path.getVertex(left) << std::endl;
                    frozen = left;
                    ResRecord<htype> *prev_ptr = &prev;
                    for(size_t i = state.last_match; i + 1 < left; i++) {
                        prev_ptr = &res.emplace(State<htype>(i + 1, 0, &path.getVertex(i + 1), true),
                                ResRecord<htype>(prev_ptr, &path[i], true)).first->second;
                    }
                    queue.emplace(State<htype>(left, 0, &path.getVertex(left), true),
                                                ResRecord<htype>(prev_ptr, &path[left - 1], true), next.score);
                    continue;
                }

            }
            if(next.state.match) {
                queue.emplace(State<htype>(next.state.last_match + 1, 0, &path.getVertex(next.state.last_match + 1), true),
                        ResRecord<htype>(&prev, &path[next.state.last_match], path[next.state.last_match].getCoverage() > scores.cov_threshold),
                        next.score + scores.scoreMatchingEdge(path[next.state.last_match]));
//                std::cout << "Add perfect " << next.score + scores.scoreMatchingEdge(path[next.state.last_match]) << std::endl;
            } else {
                size_t path_dist = 0;
                size_t extra_score = 0;
                for(size_t i = next.state.last_match; i < path.size(); i++) {
                    path_dist += path[i].size();
                    extra_score += scores.scoreRemovedEdge(path[i]);
                    if(path_dist > next.state.diff + scores.max_diff) {
                        break;
                    }
                    if(path_dist > next.state.diff - scores.max_diff && &path.getVertex(i + 1) == next.state.graph_position) {
                        size_t diff = std::min<size_t>(path_dist - next.state.diff, next.state.diff - path_dist);
                        queue.emplace(State<htype>(i + 1, 0, next.state.graph_position, true),
                                      ResRecord<htype>(&prev, nullptr, false),
                                      next.score + extra_score + diff * scores.diff_penalty);
//                        std::cout << "Add back to good " << next.score + extra_score + diff * scores.diff_penalty << std::endl;
                    }
                }
            }
            for(size_t i = 0; i < next.state.graph_position->outDeg(); i++) {
                Edge<htype> &edge = next.state.graph_position->getOutgoing()[i];
                if (!next.state.match || edge.seq != path[next.state.last_match].seq) {
                    queue.emplace(State<htype>(next.state.last_match, next.state.diff + edge.size(), edge.end(), false),
                                  ResRecord<htype>(&prev, &edge, false),
                                  next.score + scores.scoreAlternativeEdge(edge));
//                    std::cout << "Add bad " << next.score + scores.scoreAlternativeEdge(edge) << std::endl;
                }
            }

            cnt += 1;
            if (cnt >= 100000) {
                break;
            }
        }
//        std::cout << "Finished fail " << cnt << std::endl;
        return {initial_path, size_t(-1) >> 1u, cnt};
    }

    template<typename htype, class Iterator>
    void correctSequences(SparseDBG<htype> &sdbg, logging::Logger &logger, Iterator begin, Iterator end,
          const std::experimental::filesystem::path& output_file,
          const std::experimental::filesystem::path& bad_file,
          size_t threads, const size_t min_read_size) {
    //        threads = 1;
    //        omp_set_num_threads(1);
        typedef typename Iterator::value_type ContigType;
        logger.info() << "Starting to correct reads" << std::endl;
        ParallelRecordCollector<size_t> times(threads);
        ParallelRecordCollector<size_t> scores(threads);
        ParallelRecordCollector<Contig> result(threads);
        ParallelRecordCollector<Contig> prev_result(threads);
        ParallelRecordCollector<std::pair<Contig, size_t>> bad_reads(threads);
        std::unordered_map<std::string, size_t> read_order;
        std::unordered_map<std::string, size_t> prev_read_order;
        std::ofstream os;
        os.open(output_file);

        std::function<void(ContigType &)> task = [&sdbg, &times, &scores, min_read_size, &result,  &bad_reads](ContigType & contig) {
            Sequence seq = contig.makeSequence();
            if(seq.size() >= min_read_size) {
                Path<htype> path = sdbg.align(seq).path();
                CorrectionResult<htype> res = correct(path);
                times.emplace_back(res.iterations);
                scores.emplace_back(res.score);
                result.emplace_back(res.path.Seq(), contig.id + " " + std::to_string(res.score));
                if(res.score > 25000)
                    bad_reads.emplace_back(Contig(contig.makeSequence(), contig.id + " " + std::to_string(res.score)), res.score);
            } else {
                result.emplace_back(seq, contig.id + " 0");
            }
        };

        ParallelProcessor<StringContig> processor(task, logger, threads);
        processor.doInOneThread = [&read_order](ContigType & contig) {
            size_t val = read_order.size();
            read_order[split(contig.id)[0]] = val;
        };
        processor.doAfter = [&read_order, &prev_read_order, &result, &prev_result]() {
            std::swap(read_order, prev_read_order);
            std::swap(result, prev_result);
            read_order.clear();
            result.clear();
        };
        processor.doInParallel = [&prev_read_order, &prev_result, &os]() {
            VERIFY(prev_read_order.size() == prev_result.size());
            std::vector<Contig> res(prev_read_order.size());
            for(Contig &contig : prev_result) {
                res[prev_read_order[split(contig.id)[0]]] = std::move(contig);
            }
            for(Contig &contig : res) {
                os << ">" << contig.id << "\n" << contig.seq << "\n";
            }
        };

        processor.doInTheEnd = processor.doInParallel;

        processor.processRecords(begin, end);

        os.close();
        std::vector<size_t> time_hist = histogram(times.begin(), times.end(), 100000, 1000);
        std::vector<size_t> score_hist = histogram(scores.begin(), scores.end(), 100000, 1000);
        logger.info() << "Reader corrected" << std::endl;
        logger.info() << "Iterations" << std::endl;
        for(size_t i = 0; i < time_hist.size(); i++) {
            logger << i * 1000 << " " << time_hist[i] << std::endl;
        }
        logger.info() << "Scores" << std::endl;
        for(size_t i = 0; i < score_hist.size(); i++) {
            logger << i * 1000 << " " << score_hist[i] << std::endl;
        }
        logger.info() << "Printed corrected reads to " << output_file << std::endl;
        logger.info() << "Found " << bad_reads.size() << " bad reads" << std::endl;
        std::ofstream bados;
        bados.open(bad_file);
        for(auto & contig : bad_reads) {
            bados << ">" << contig.first.id << " " << contig.second << "\n" << contig.first.seq << "\n";
        }
        bados.clear();
        logger.info() << "Printed bad reads to " << bad_file << std::endl;
    }
}
