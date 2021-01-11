//
// Created by anton on 7/22/20.
//

#pragma once
#include "sequences/sequence.hpp"
#include "sequences/seqio.hpp"
#include "common/omp_utils.hpp"
#include "common/logging.hpp"
#include "rolling_hash.hpp"
#include "hash_utils.hpp"
#include <vector>
#include <numeric>
#include <unordered_map>
#include <common/oneline_utils.hpp>
#include <unordered_set>

template<class htype>
class Vertex;

template<class htype>
class Edge {
private:
    Vertex<htype> *end_;
    mutable size_t cov;
public:
    mutable size_t extraInfo;
    Sequence seq;
    std::string id;
    friend class Vertex<htype>;

    Edge(Vertex<htype> *_end, const Sequence &_seq) :
            end_(_end), cov(0), extraInfo(-1), seq(_seq), id("") {
    }

    Vertex<htype> *end() const {
        return end_;
    }

    size_t getTipSize() const {
        return extraInfo;
    }

    size_t updateTipSize() const {
        size_t new_val = 0;
        if(extraInfo == size_t(-1) && end_->inDeg() == 1) {
            for (const Edge<htype> & other : end_->getOutgoing()) {
                other.end_->lock();
                new_val = std::max(new_val, other.extraInfo);
                other.end_->unlock();
            }
            if(new_val != size_t(-1))
                new_val += size();
            end_->lock();
            extraInfo = new_val;
            end_->unlock();
        }
        return new_val;
    }

    void bindTip(Vertex<htype> & start, Vertex<htype> & end) {
        VERIFY(end_ == nullptr);
        end_ = &end;
        Sequence rcseq = !(start.seq + seq);
        end.rc().addEdgeLockFree(Edge<htype>(&start.rc(), rcseq.Subseq(start.seq.size())));
    }

    size_t common(Sequence other) const {
        size_t res = 0;
        while(res < seq.size() && res < other.size() && seq[res] == other[res]) {
            res += 1;
        }
        return res;
    }

    size_t size() const {
        return seq.size();
    }

    double getCoverage() const {
        return double(cov) / size();
    }

    double intCov() const {
        return cov;
    }


    void incCov(size_t val) const {
#pragma omp atomic
        cov += val;
    }

    bool operator==(const Edge<htype> &other) const {
        return this == &other;
    }

    bool operator<(const Edge<htype> &other) const {
        return this->seq < other.seq;
    }
};

template<typename htype>
class SparseDBG;

template<class htype>
class Vertex {
private:
    friend class SparseDBG<htype>;
    std::vector<Edge<htype>> outgoing_{};
    Vertex * rc_;
    htype hash_;
    omp_lock_t writelock{};
    size_t coverage_ = 0;
    bool canonical = false;
    bool mark_ = false;


    explicit Vertex(htype hash, Vertex *_rc) : hash_(hash), rc_(_rc), canonical(false) {
        omp_init_lock(&writelock);
    }

public:
    Sequence seq;

    size_t converage() const {
        return coverage_;
    }

    bool isCanonical() const {
        return canonical;
    }

    bool isCanonical(const Edge<htype> &edge) const {
        const Vertex<htype> &other = edge.end()->rc();
        if(hash() != other.hash())
            return hash() < other.hash();
        if (isCanonical() != other.isCanonical())
            return isCanonical();
        const Edge<htype> &rc_edge = rcEdge(edge);
        return edge.seq <= rc_edge.seq;
    }

    std::string edgeId(const Edge<htype> &edge) const {
        std::stringstream ss;
        if(isCanonical(edge)) {
            ss << hash() << isCanonical() << "ACGT"[edge.seq[0]];
            return ss.str();
        } else {
            return edge.end()->rc().edgeId(rcEdge(edge));
        }
    }

    void mark() {
        mark_ = true;
    }

    void unmark() {
        mark_ = false;
    }

    bool marked() const {
        return mark_;
    }

    void clear() {
        outgoing_.clear();
        rc_->outgoing_.clear();
    }

    explicit Vertex(htype hash = 0) : hash_(hash), rc_(new Vertex<htype>(hash, this)), canonical(true) {
        omp_init_lock(&writelock);
    }

    Vertex(Vertex &&other) noexcept : rc_(other.rc_), hash_(other.hash_), canonical(other.canonical) {
        std::swap(outgoing_, other.outgoing_);
        std::swap(seq, other.seq);
        omp_init_lock(&writelock);
        if(other.rc_ != nullptr) {
            other.rc_->rc_ = this;
            other.rc_ = nullptr;
        }
        for(Edge<htype> &edge : rc_->outgoing_) {
            if(edge.end() == &other) {
                edge.end_ = this;
            } else if (edge.end() != nullptr){
                for(Edge<htype> &back : edge.end()->rc().outgoing_) {
                    if(back.end() == &other) {
                        back.end_ = this;
                    }
                }
            }
        }
    }

    void sortOutgoing() {
        std::sort(outgoing_.begin(), outgoing_.end());
    }

    Vertex(const Vertex &) = delete;

    ~Vertex() {
        if(rc_ != nullptr) {
            rc_->rc_ = nullptr;
            delete rc_;
        }
        rc_ = nullptr;
    }

    void checkConsistency() const {
        for(const Edge<htype> & edge : outgoing_) {
            if(edge.end() != nullptr) {
                if(this->rcEdge(edge).end() != &(this->rc())) {
                    std::cout << this << " " << seq << " " << edge.seq << " " << rcEdge(edge).end() << " " << &(this->rc()) <<std::endl;
                }
                VERIFY(this->rcEdge(edge).end() == &(this->rc()));
            }
        }
    }

    htype hash() const {
        return hash_;
    }

    const Vertex<htype> & rc() const {
        return *rc_;
    }

    Vertex<htype> & rc() {
        return *rc_;
    }

    void incCoverage() {
#pragma omp atomic update
        coverage_ += 1;
#pragma omp atomic update
        rc().coverage_ += 1;
    }

    Edge<htype>& rcEdge(const Edge<htype> & edge) {
        Vertex &vend = edge.end()->rc();
        char c;
        if(edge.size() > seq.size()) {
            c = (!edge.seq)[seq.size()];
        } else {
            c = (!seq)[seq.size() - edge.size()];
        }
        return vend.getOutgoing(c);
    }

    Edge<htype>& sparseRcEdge(const Edge<htype> & edge) {
        Vertex &vend = edge.end()->rc();
        VERIFY(seq.size() > 0);
        for(Edge<htype> &candidate : vend.getOutgoing()) {
            if(candidate.end() == rc_ && candidate.size() == edge.size() &&
                    (edge.size() <= seq.size() || candidate.seq.startsWith((!edge.seq).Subseq(seq.size())))) {
                return candidate;
            }
        }
        std::cout << seq + edge.seq << std::endl;
        std::cout << seq << std::endl;
        std::cout << vend.seq << std::endl;
        for(Edge<htype> &candidate : vend.getOutgoing()) {
            if(candidate.end() == rc_ && candidate.size() == edge.size() &&
               (edge.size() <= seq.size() || candidate.seq.startsWith((!edge.seq).Subseq(seq.size())))) {
                std::cout << vend.seq + candidate.seq << std::endl;
            }
        }
        VERIFY(false);
        return vend.getOutgoing()[0];
    }

    const Edge<htype>& sparseRcEdge(const Edge<htype> & edge) const {
        Vertex &vend = edge.end()->rc();
        VERIFY(seq.size() > 0);
        for(Edge<htype> &candidate : vend.getOutgoing()) {
            if(candidate.end() == rc_ && candidate.size() == edge.size() &&
               (edge.size() <= seq.size() || candidate.seq.startsWith((!edge.seq).Subseq(seq.size())))) {
                return candidate;
            }
        }
        std::cout << seq + edge.seq << std::endl;
        std::cout << seq << std::endl;
        std::cout << vend.seq << std::endl;
        for(Edge<htype> &candidate : vend.getOutgoing()) {
            if(candidate.end() == rc_ && candidate.size() == edge.size() &&
               (edge.size() <= seq.size() || candidate.seq.startsWith((!edge.seq).Subseq(seq.size())))) {
                std::cout << vend.seq + candidate.seq << std::endl;
            }
        }
        VERIFY(false);
        return vend.getOutgoing()[0];
    }

    const Edge<htype>& rcEdge(const Edge<htype> & edge) const {
        const Vertex &vend = edge.end()->rc();
        char c;
        if(edge.size() > seq.size()) {
            c = (!edge.seq)[seq.size()];
        } else {
            c = (!seq)[seq.size() - edge.size()];
        }
        return vend.getOutgoing(c);
    }

    Sequence pathSeq(const std::vector<Edge<htype> *> & path) const {
        SequenceBuilder sb;
        sb.append(seq);
        for(const Edge<htype> *e : path) {
            sb.append(e->seq);
        }
        return sb.BuildSequence();
    }

    std::vector<Edge<htype>*> walkForward(Edge<htype> & edge) {
        std::vector<Edge<htype>*> res;
        res.push_back(&edge);
        Vertex<htype> *next = edge.end();
        while(next != nullptr && next != this && !next->isJunction()) {
            res.push_back(&next->getOutgoing()[0]);
            next = next->getOutgoing()[0].end();
        }
        return std::move(res);
    }

    void setSequence(const Sequence &_seq) {
        lock();
        if(seq.empty()) {
            if(seq.empty()) {
                seq = Sequence(_seq.str());
                unlock();
                rc_->lock();
                rc_->seq = !_seq;
                rc_->unlock();
            } else {
                unlock();
            }
        } else {
//            if(seq != _seq) {
//                std::cout << seq << std::endl << _seq << std::endl;
//                VERIFY(false);
//            }
//            VERIFY(_seq == seq);
            unlock();
        }
    }

    void clearSequence() {
        if (!seq.empty()) {
            seq = Sequence();
            rc_->seq = Sequence();
        }
    }

    Edge<htype> & addEdgeLockFree(const Edge<htype> &edge) {
        for(Edge<htype> & e : outgoing_) {
            if (edge.size() <= e.size()) {
                if (edge.seq == e.seq.Subseq(0, edge.size())) {
                    return e;
                }
            } else if (edge.seq.Subseq(0, e.size()) == e.seq) {
                e = edge;
                return e;
            }
        }
        outgoing_.emplace_back(edge);
        return outgoing_.back();
    }

    void addEdge(const Edge<htype> &e) {
        omp_set_lock(&writelock);
        addEdgeLockFree(e);
        omp_unset_lock(&writelock);
    }

    void removeEdgesTo(const Vertex<htype> & other) {
        lock();
        outgoing_ = oneline::filter(other.begin(), other.end(), [&](const Edge<htype> & edge) {return edge.end() != &other;});
        unlock();
    }

    void lock() {
        omp_set_lock(&writelock);
    }

    void unlock() {
        omp_unset_lock(&writelock);
    }

    const std::vector<Edge<htype>> &getOutgoing() const {
        return outgoing_;
    }

    const Edge<htype> & getOutgoing(char c) const {
        for(const Edge<htype> &edge : outgoing_) {
            if(edge.seq[0] == c) {
                return edge;
            }
        }
        std::cout << seq << std::endl;
        std::cout << c << std::endl;
        for(const Edge<htype> &edge : outgoing_) {
            std::cout << edge.seq << std::endl;
        }
        VERIFY(false);
        return getOutgoing()[0];
    }

    Edge<htype> & getOutgoing(char c) {
        for(Edge<htype> &edge : outgoing_) {
            if(edge.seq[0] == c) {
                return edge;
            }
        }
        std::cout << seq << std::endl;
        std::cout << size_t(c) << std::endl;
        for(Edge<htype> &edge : outgoing_) {
            std::cout << edge.seq << std::endl;
        }
        VERIFY(false);
        return getOutgoing()[0];
    }

    bool hasOutgoing(char c) const {
        for(const Edge<htype> &edge : outgoing_) {
            if(edge.seq[0] == c) {
                return true;
            }
        }
        return false;
    }

    std::vector<Edge<htype>> &getOutgoing() {
        return outgoing_;
    }

    size_t outDeg() const {
        return outgoing_.size();
    }

    size_t inDeg() const {
        return rc_->outgoing_.size();
    }

    bool isJunction() const {
        return outDeg() != 1 || inDeg() != 1;
    }

    bool operator==(const Vertex<htype> & other) const {
        return this == &other;
    }

    bool operator!=(const Vertex<htype> & other) const {
        return this != &other;
    }

    bool operator<(const Vertex<htype> & other) const {
        return hash_ < other.hash_ || (hash_ == other.hash_ && canonical && !other.canonical);
    }
};


template<typename htype>
class Path {
private:
    Vertex<htype> *start_;
    std::vector<Edge<htype> *> path;
public:
    Path(Vertex<htype> &_start, std::vector<Edge<htype> *> _path) : start_(&_start), path(std::move(_path)){}
    Path(Vertex<htype> &_start) : start_(&_start) {}

    Path<htype> subPath(size_t from, size_t to) {
        if(from == to)
            return Path(getVertex(from));
        else
            return Path(getVertex(from), std::vector<Edge<htype> *>(path.begin() + from, path.begin() + to));
    }

    Path<htype> RC() {
        std::vector<Edge<htype> *> rcPath;
        for(size_t i = path.size(); i > 0; i--) {
            rcPath.emplace_back(&getVertex(i - 1).rcEdge(*path[i - 1]));
        }
        return Path(back().end()->rc(), rcPath);
    }

    Edge<htype> &operator[](size_t i) const {
        return *path[i];
    }

    Edge<htype> &operator[](size_t i) {
        return *path[i];
    }

    Vertex<htype> &getVertex(size_t i) {
        VERIFY(i <= path.size());
        if(i == 0)
            return *start_;
        else
            return *path[i - 1]->end();
    }

    Vertex<htype> &getVertex(size_t i) const {
        VERIFY(i <= path.size());
        if(i == 0)
            return *start_;
        else
            return *path[i - 1]->end();
    }

    Vertex<htype> &start() const {
        return *start_;
    }

    Vertex<htype> &finish() const {
        if(path.size() == 0)
            return start();
        else
            return *path.back()->end();
    }

    Edge<htype> &back() {
        return *path.back();
    }

    Edge<htype> &front() {
        return *path.front();
    }

    double minCoverage() const {
        double res = 100000;
        for (const Edge<htype> *edge : path) {
            res = std::min(edge->getCoverage(), res);
        }
        return res;
    }

    Sequence Seq() const {
        return start_->pathSeq(path);
    }

    Sequence truncSeq() const {
        SequenceBuilder sb;
        for(const Edge<htype> *e : path) {
            sb.append(e->seq);
        }
        return sb.BuildSequence();
    }

    size_t size() const {
        return path.size();
    }

    typedef typename std::vector<Edge<htype> *>::iterator iterator;
    typedef typename std::vector<Edge<htype> *>::const_iterator const_iterator;
    iterator begin() {
        return path.begin();
    }
    iterator end() {
        return path.end();
    }
    const_iterator begin() const {
        return path.begin();
    }
    const_iterator end() const {
        return path.end();
    }

    size_t len() const {
        size_t res = 0;
        for(Edge<htype> *edge : path)
            res += edge->size();
        return res;
    }

    Path<htype> operator+(const Path<htype> &other) const {
        VERIFY(finish() == *other.start_);
        std::vector<Edge<htype> *> edges = path;
        edges.insert(edges.end(), other.path.begin(), other.path.end());
        return {start(), std::move(edges)};
    }
};


template<typename htype>
class GraphAlignment {
private:
    Vertex<htype> *start_;
    std::vector<Segment<Edge<htype>>> als;
public:
//    TODO change interface
    GraphAlignment(Vertex<htype> *_start, std::vector<Segment<Edge<htype>>> &&_path) : start_(_start), als(std::move(_path)){}

    explicit GraphAlignment(Vertex<htype> &_start) : start_(&_start) {}

    GraphAlignment() : start_(nullptr) {}

    GraphAlignment RC() {
        std::vector<Segment<Edge<htype>>> path;
        for(size_t i = 0; i < als.size(); i++) {
            Edge<htype> &rc_edge = getVertex(i).rcEdge(als[i].contig());
            path.emplace_back(rc_edge, rc_edge.size() - als[i].right, rc_edge.size() - als[i].left);
        }
        GraphAlignment<htype> res= {&finish().rc(), {path.rbegin(), path.rend()}};
        return res;
    }

    void invalidate() {
        start_ = nullptr;
        als.clear();
    }

    bool valid() const {
        return start_ != nullptr;
    }

    void push_back(const Segment<Edge<htype>> &seg) {
        als.push_back(seg);
    }

    void pop_back() {
        als.pop_back();
    }

    void pop_back(size_t len) {
        als.erase(als.end() - len, als.end());
    }

    void cutBack(size_t l) {
        VERIFY(l <= len());
        while(size() > 0 && als.back().size() < l) {
            l -= als.back().size();
            pop_back();
        }
        if(l > 0) {
            als.back().right -= l;
        }
    }


    GraphAlignment<htype> subalignment(size_t from, size_t to) {
        return {&getVertex(from), std::vector<Segment<Edge<htype>>>(als.begin() + from, als.begin() + to)};
    }

    GraphAlignment &addStep() {
        als.back().right += 1;
        return *this;
    }

    GraphAlignment & addStep(Edge<htype> &edge) {
        als.emplace_back(edge, 0, 1);
        return *this;
    }

    GraphAlignment & extend(const Sequence &seq) {
        VERIFY(valid());
        for(size_t cpos = 0; cpos < seq.size(); cpos++) {
            unsigned char c = seq[cpos];
            if(endClosed()) {
                Vertex<htype> &v = finish();
                if(v.hasOutgoing(c)) {
                    Edge<htype> &edge = v.getOutgoing(c);
                    addStep(edge);
                } else {
                    invalidate();
                    return *this;
                }
            } else {
                if(als.back().contig().seq[als.back().right] == c) {
                    addStep();
                } else {
                    invalidate();
                    return *this;
                }
            }
        }
        return *this;
    }

    bool endClosed() const {
        return start_ != nullptr && (als.size() == 0 || als.back().right == als.back().contig().size());
    }

    bool startClosed() const {
        return start_ != nullptr && (als.size() == 0 || als.front().left == 0);
    }

    unsigned char lastNucl() const {
        VERIFY(als.back().size() > 0);
        return als.back().contig().seq[als.back().right - 1];
    }

    size_t leftSkip() const {
        return als.size() == 0 ? 0 : als.front().left;
    }

    size_t rightSkip() const {
        return als.size() == 0 ? 0 : als.back().contig().size() - als.back().right;
    }

    std::vector<GraphAlignment<htype>> allSteps() {
        if(als.size() != 0 && als.back().right < als.back().contig().size()) {
            GraphAlignment<htype> copy = *this;
            return {std::move(copy.addStep())};
        }
        std::vector<GraphAlignment<htype>> res;
        Vertex<htype> &end = als.size() == 0 ? *start_ : *back().contig().end();
        for(Edge<htype> &edge : end.getOutgoing()) {
            GraphAlignment<htype> copy = *this;
            res.emplace_back(std::move(copy.addStep(edge)));
        }
        return res;
    }

    std::vector<GraphAlignment<htype>> allExtensions(size_t len) {
        std::vector<GraphAlignment<htype>> res = {*this};
        size_t left = 0;
        size_t right = 1;
        for(size_t l = 0; l < len; l++) {
            for(size_t i = left; i < right; i++) {
                std::vector<GraphAlignment<htype>> tmp = res[i].allSteps();
                res.insert(res.end(), tmp.begin(), tmp.end());
            }
            left = right;
            right = res.size();
        }
        return std::move(res);
    }

    Sequence map(std::unordered_map<const Edge<htype128> *, Sequence> &edge_map) {
        SequenceBuilder sb;
        bool start = true;
        for(Segment<Edge<htype>> & seg : als) {
            auto it = edge_map.find(&seg.contig());
            if(it == edge_map.end()) {
                if(start) {
                    sb.append((start_->seq + seg.contig().seq).Subseq(seg.left, seg.right + start_->seq.size()));
                    start = false;
                } else {
                    sb.append(seg.seq());
                }
            } else {
                size_t left = start_->seq.size();
                if (start) {
                    left = 0;
                }
                size_t right = start_->seq.size();
                size_t sz = it->second.size() - start_->seq.size();
                if (seg.left == 0 && seg.right == seg.contig().size()) {
                    right += sz;
                } else if (seg.left == 0) {
                    right += std::min(sz, seg.right);
                } else if(seg.right == seg.contig().size()) {
                    left += sz - std::min(sz, seg.size());
                    right += sz;
                } else {
                    size_t l = seg.left * sz / seg.contig().size();
                    left += l;
                    right += std::min(l + seg.size(), sz);
                }
                sb.append(it->second.Subseq(left, right));
                start = false;
            }
        }
        return sb.BuildSequence();
    }

    Sequence Seq() const {
        if(als.size() == 0) {
            return {};
        }
        SequenceBuilder sb;
        size_t k = start_->seq.size();
        if(als[0].left >= k)
            sb.append(als[0].contig().seq.Subseq(als[0].left - k, als[0].right));
        else {
            sb.append(start_->seq.Subseq(als[0].left, k));
            sb.append(als[0].contig().seq.Subseq(0, als[0].right));
        }
        for(size_t i = 1; i < als.size(); i++) {
            sb.append(als[i].seq());
        }
        return sb.BuildSequence();
    }

    Sequence truncSeq()const {
        SequenceBuilder sb;
        for(size_t i = 0; i < als.size(); i++) {
            sb.append(als[i].seq());
        }
        return sb.BuildSequence();
    }

    Sequence truncSeq(size_t start_position, size_t size) const {
        SequenceBuilder sb;
        size_t sz = 0;
        for(size_t i = start_position; i < als.size(); i++) {
//            std::cout << i << " " << sz << " " << size << std::endl;
//            std::cout << als[i].contig().size() << " " << als[i].left << " " << als[i].right << " " << als[i].size() << std::endl;
            if(sz + als[i].size() >= size) {
                sb.append(als[i].seq().Subseq(0, size - sz));
                break;
            } else {
                sb.append(als[i].seq());
                sz += als[i].size();
            }
        }
        return sb.BuildSequence();
    }

    Vertex<htype> &start() const {
        return *start_;
    }

    Vertex<htype> &finish() const {
        return als.size() == 0 ? *start_ : *als.back().contig().end();
    }

    Segment<Edge<htype>> &back() {
        return als.back();
    }

    Segment<Edge<htype>> &front() {
        return als.front();
    }

    const Segment<Edge<htype>> &operator[](size_t i) const {
        return als[i];
    }

    Segment<Edge<htype>> &operator[](size_t i) {
        return als[i];
    }

    Vertex<htype> &getVertex(size_t i) const {
        VERIFY(i <= als.size());
        if(i == 0)
            return *start_;
        else
            return *als[i - 1].contig().end();
    }

    Vertex<htype> &getVertex(size_t i) {
        VERIFY(i <= als.size());
        if(i == 0)
            return *start_;
        else
            return *als[i - 1].contig().end();
    }

    typename std::vector<Segment<Edge<htype>>>::iterator begin() {
        return als.begin();
    }

    typename std::vector<Segment<Edge<htype>>>::iterator end() {
        return als.end();
    }

    typename std::vector<Segment<Edge<htype>>>::const_iterator begin() const {
        return als.begin();
    }

    typename std::vector<Segment<Edge<htype>>>::const_iterator end() const {
        return als.end();
    }

    GraphAlignment<htype> subPath(size_t left, size_t right) const {
        if(left == right)
            return GraphAlignment<htype>(getVertex(left));
        else
            return {&getVertex(left), {als.begin() + left, als.begin() + right}};
    }

    GraphAlignment<htype> reroute(size_t left, size_t right, const GraphAlignment<htype> & rerouting) const {
        VERIFY(getVertex(left) == rerouting.start());
        VERIFY(getVertex(right) == rerouting.finish());
        return subPath(0, left) + rerouting + subPath(right, size());
    }

    GraphAlignment<htype> operator+(const GraphAlignment<htype> &other) const {
        VERIFY(finish() == other.getVertex(0));
        std::vector<Segment<Edge<htype>>> new_als(als.begin(), als.end());
        new_als.insert(new_als.end(), other.begin(), other.end());
        return {&start(), std::move(new_als)};
    }

    double minCoverage() const {
        double res = 100000;
        for (const Segment<Edge<htype>> &seg : als) {
            res = std::min(seg.contig().getCoverage(), res);
        }
        return res;
    }

    Path<htype> path() {
        std::vector<Edge<htype>*> res;
        for(auto & seg : als) {
            res.push_back(&seg.contig());
        }
        return {*start_, res};
    }

    size_t size() const {
        return als.size();
    }

    size_t len() const {
        size_t res = 0;
        for(auto & seg : als) {
            res += seg.size();
        }
        return res;
    }

    bool operator==(const GraphAlignment<htype> &other) const {
        return start_ == other.start_ && als == other.als;
    }

    bool operator!=(const GraphAlignment<htype> &other) const {
        return start_ != other.start_ || als != other.als;
    }
};

template<typename htype>
struct EdgePosition {
    Vertex<htype> *start;
    Edge<htype> *edge;
    size_t pos;

    EdgePosition(Vertex<htype> &_start, Edge<htype> &_edge, size_t _pos) : start(_start), edge(*_edge), pos(_pos) {
        VERIFY(pos > 0);
    }

    EdgePosition() : start(nullptr), edge(nullptr), pos(0) {
    }

    GraphAlignment<htype> align(const Sequence &seq) {
        GraphAlignment<htype> res(*start, {*edge, pos, pos});
        Edge<htype> *cedge = edge;
        size_t epos = pos;
        for(size_t cpos = 0; cpos < seq.size(); cpos++) {
            unsigned char c = seq[cpos];
            if(epos == cedge->size()) {
                Vertex<htype> &v = *cedge->end();
                if(v.hasOutgoing(c)) {
                    cedge = &v.getOutgoing(c);
                    res.addStep(cedge);
                    epos = 1;
                } else {
                    return {};
                }
            } else {
                if(cedge->seq[epos] == c) {
                    res.addStep();
                    epos += 1;
                } else {
                    return {};
                }
            }
        }
        return std::move(res);
    }

    std::vector<EdgePosition> step() {
        if(pos == edge->size()) {
            std::vector<EdgePosition> res;
            Vertex<htype> &v = *edge->end();
            for(Edge<htype> &next : v.getOutgoing()) {
                res.emplace_back(v, next, 1);
            }
            return std::move(res);
        } else {
            return {start, *edge, pos + 1};
        }
    }

    unsigned char lastNucl() const {
        return edge->seq[pos - 1];
    }
};



template<typename htype>
class SparseDBG {
public:
    struct EdgePosition {
        Edge<htype> *edge;
        Vertex<htype> *start;
        size_t pos;

        EdgePosition(Edge<htype> &edge, Vertex<htype> &start, size_t pos) : edge(&edge), start(&start), pos(pos) {}
        EdgePosition(const EdgePosition &other) : edge(other.edge), start(other.start), pos(other.pos) {}

        EdgePosition RC() {
            Vertex<htype> &s = *start;
            Edge<htype> &rc_edge = s.rcEdge(*edge);
            return EdgePosition(rc_edge, edge->end()->rc(), edge->size() - pos);
        }
    };
    typedef std::unordered_map<htype, Vertex<htype>, alt_hasher<htype>> vertex_map_type;
    typedef std::unordered_map<htype, EdgePosition, alt_hasher<htype>> anchor_map_type;
private:
//    TODO: replace with perfect hash map? It is parallel, maybe faster and compact.
    vertex_map_type v;
    anchor_map_type anchors;
    const RollingHash<htype> hasher_;

    std::vector<KWH<htype>> extractVertexPositions(const Sequence &seq) const {
        std::vector<KWH<htype>> res;
        KWH<htype> kwh(hasher_, seq, 0);
        while(true) {
            if(v.find(kwh.hash()) != v.end()) {
                res.emplace_back(kwh);
            }
            if (!kwh.hasNext())
                break;
            kwh = kwh.next();
        }
        return std::move(res);
    }

//    Be careful since hash does not define vertex. Rc vertices share the same hash
    Vertex<htype> & innerAddVertex(htype h) {
        return v.emplace(h, Vertex<htype>(h)).first->second;
    }

public:

    template<class Iterator>
    SparseDBG(Iterator begin, Iterator end, RollingHash<htype> _hasher) : hasher_(_hasher) {
        while(begin != end) {
            htype hash = *begin;
            if(v.find(hash) == v.end())
                addVertex(hash);
            ++begin;
        }
    }

    SparseDBG(RollingHash<htype> _hasher) : hasher_(_hasher) {
    }

    SparseDBG(SparseDBG<htype> &&other) noexcept : hasher_(other.hasher_) {
        std::swap(v, other.v);
        std::swap(anchors, other.anchors);
    }

    SparseDBG(const SparseDBG<htype> &other) noexcept = delete;

    SparseDBG &operator=(SparseDBG<htype> &&other) = delete;

    bool containsVertex(const htype &hash) const {
        return v.find(hash) != v.end();
    }

    void checkConsistency(size_t threads, logging::Logger &logger) {
        logger.info() << "Checking consistency" << std::endl;
        std::function<void(std::pair<const htype, Vertex<htype>> &)> task =
                [this](std::pair<const htype, Vertex<htype>> & pair) {
                    const Vertex<htype> & vert = pair.second;
                    vert.checkConsistency();
                    vert.rc().checkConsistency();
                };
        processObjects(v.begin(), v.end(), logger, threads, task);
        logger.info() << "Consistency check success" << std::endl;
    }

    void checkSeqFilled(size_t threads, logging::Logger &logger) {
        logger.info() << "Checking vertex sequences" << std::endl;
        std::function<void(std::pair<const htype, Vertex<htype>> &)> task =
                [&logger](std::pair<const htype, Vertex<htype>> & pair) {
                    const Vertex<htype> & vert = pair.second;
                    if(vert.seq.empty() || vert.rc().seq.empty()) {
                        logger.info() << "Sequence not filled " << pair.first << std::endl;
                        VERIFY(false);
                    }
                    if(!vert.isCanonical()) {
                        logger.info() << "Canonical vertex marked not canonical " << pair.first << std::endl;
                        VERIFY(false);
                    }
                    if(vert.rc().isCanonical()) {
                        logger.info() << "Noncanonical vertex marked canonical " << pair.first << std::endl;
                        VERIFY(false);
                    }
                };
        processObjects(v.begin(), v.end(), logger, threads, task);
        logger.info() << "Vertex sequence check success" << std::endl;
    }

    const RollingHash<htype> &hasher() const {
        return hasher_;
    }

    void addVertex(htype h) {
        innerAddVertex(h);
    }

    Vertex<htype> &addVertex(const KWH<htype> &kwh) {
        Vertex<htype> &newVertex = innerAddVertex(kwh.hash());
        Vertex<htype> &res = kwh.isCanonical() ? newVertex : newVertex.rc();
        res.setSequence(kwh.getSeq());
        return res;
    }

    Vertex<htype> &addVertex(const Sequence &seq) {
        return addVertex(KWH<htype>(hasher_, seq, 0));
    }

    Vertex<htype> &bindTip(Vertex<htype> & start, Edge<htype> &tip) {
        Sequence seq = start.seq + tip.seq;
        Vertex<htype> &end = addVertex(seq.Subseq(seq.size() - hasher().k));
        tip.bindTip(start, end);
        return end;
    }

    Vertex<htype> &getVertex(const KWH<htype> &kwh) {
        VERIFY(v.find(kwh.hash()) != v.end());
        if(kwh.isCanonical()) {
            return v[kwh.hash()];
        } else {
            return v[kwh.hash()].rc();
        }
    }

    Vertex<htype> &getVertex(htype hash) {
        return v.find(hash)->second;
    }

    const Vertex<htype> &getVertex(const KWH<htype> &kwh) const {
        VERIFY(v.find(kwh.hash()) != v.end());
        if(kwh.isCanonical()) {
            return v.find(kwh.hash())->second;
        } else {
            return v.find(kwh.hash())->second.rc();
        }
    }

    void fillAnchors(size_t w, logging::Logger &logger, size_t threads) {
        logger.info() << "Adding anchors from long edges for alignment" << std::endl;
        ParallelRecordCollector<std::pair<const htype, EdgePosition>> res(threads);
        std::function<void(std::pair<const htype, Vertex<htype>> &)> task = [&res, w, this] (std::pair<const htype, Vertex<htype>> &iter) {
            Vertex<htype> &vertex = iter.second;
            for(Edge<htype> &edge : vertex.getOutgoing()) {
                if(edge.size() > w) {
                    Sequence seq = vertex.seq + edge.seq;
//                    Does not run for the first and last kmers.
                    for(KWH<htype> kmer(this->hasher_, seq, 1); kmer.hasNext(); kmer = kmer.next()) {
                        if (kmer.pos % w == 0) {
                            EdgePosition ep(edge, vertex, kmer.pos);
                            if(kmer.isCanonical())
                                res.emplace_back(kmer.hash(), ep);
                            else {
                                res.emplace_back(kmer.hash(), ep.RC());
                            }
                        }
                    }
                }
            }
            for(Edge<htype> &edge : vertex.rc().getOutgoing()) {
                if(edge.size() > w) {
                    Sequence seq = vertex.rc().seq + edge.seq;
//                    Does not run for the first and last kmers.
                    for(KWH<htype> kmer(this->hasher_, seq, 1); kmer.hasNext(); kmer = kmer.next()) {
                        if (kmer.pos % w == 0) {
                            EdgePosition ep(edge, vertex.rc(), kmer.pos);
                            if(kmer.isCanonical())
                                res.emplace_back(kmer.hash(), ep);
                            else {
                                res.emplace_back(kmer.hash(), ep.RC());
                            }
                        }
                    }
                }
            }
        };
        processObjects(begin(), end(), logger, threads, task);
        for(auto & tmp : res) {
            anchors.emplace(tmp);
        }
        logger.info() << "Added " << anchors.size() << " anchors" << std::endl;
    }

    bool isAnchor(htype hash) const {
        return anchors.find(hash) != anchors.end();
    }

    EdgePosition getAnchor(const KWH<htype> &kwh) {
        if(kwh.isCanonical())
            return anchors.find(kwh.hash())->second;
        else
            return anchors.find(kwh.hash())->second.RC();
    }

//    GraphAlignment<htype> partialAlign(const Sequence & seq) {
//        std::vector<KWH<htype>> kmers = extractVertexPositions(seq);
//        std::vector<Segment<Edge<htype>>> res;
//        if(kmers.size() == 0) {
//            KWH<htype> kwh(hasher_, seq, 0);
//            while(true) {
//                if(isAnchor(kwh.hash())) {
//                    EdgePosition pos = getAnchor(kwh);
//                    if(kwh.pos >= pos.pos || pos.pos + seq.size() - kwh.pos > pos.edge->size() + hasher_.k) {
//                        return {nullptr, std::move(res)};
//                    }
//                    Segment<Edge<htype>> seg(*pos.edge, pos.pos - kwh.pos, pos.pos + seq.size() - kwh.pos - hasher_.k);
//                    if (seg.seq() == seq)
//                        return {pos.start, std::vector<Segment<Edge<htype>>>({seg})};
//                    else
//                        return {nullptr, std::move(res)};
//                }
//                if (!kwh.hasNext()) {
//                    return {nullptr, std::move(res)};
//                }
//                kwh = kwh.next();
//            }
//        }
//        Vertex<htype> *prestart = &getVertex(kmers.front());
//        if (kmers.front().pos > 0) {
//            const Vertex<htype> &rcstart = prestart->rc();
//            if(rcstart.hasOutgoing(seq[kmers.front().pos - 1] ^ 3)) {
//                const Edge<htype> &rcedge = rcstart.getOutgoing(seq[kmers.front().pos - 1] ^ 3);
//                const Edge<htype> &edge = rcstart.rcEdge(rcedge);
//                if (edge.size() >= kmers.front().pos) {
//                    Segment<Edge<htype>> seg(edge, edge.size() - kmers.front().pos, edge.size());
//                    if(rcedge.seq.Subseq(0, kmers.front().pos == !(seq.Subseq(0, kmers.front().pos))) {
//                        res.emplace_back(seg);
//                        prestart = &rcedge.end()->rc();
//                    }
//                }
//            }
//        }
//        for(const KWH<htype> & kmer : kmers) {
//            if (kmer.pos + hasher_.k < seq.size()) {
//                Vertex<htype> &vertex = getVertex(kmer);
//                if(vertex.hasOutgoing(seq[kmer.pos + hasher_.k])) {
//                    const Edge<htype> &edge = vertex.getOutgoing(seq[kmer.pos + hasher_.k]);
//                    Segment<Edge<htype>> seg(edge, 0, std::min(seq.size() - kmer.pos - hasher_.k, edge.size()));
//                    if(seg.seq() == seq.Subseq(kmer.pos + hasher_.k, std::min(kmer.pos + hasher_.k + seg.size(), seq.size()))){
//                        if (res.empty())
//                            prestart = &vertex;
//                        res.emplace_back(seg);
//                    }
//                }
//            }
//        }
//        return {prestart, std::move(res)};
//    }

    GraphAlignment<htype> align(const Sequence & seq) {
        std::vector<KWH<htype>> kmers = extractVertexPositions(seq);
        std::vector<Segment<Edge<htype>>> res;
        if(kmers.size() == 0) {
            KWH<htype> kwh(hasher_, seq, 0);
            while(true) {
                if(isAnchor(kwh.hash())) {
                    EdgePosition pos = getAnchor(kwh);
                    VERIFY(kwh.pos < pos.pos);
                    VERIFY(pos.pos + seq.size() - kwh.pos <= pos.edge->size() + hasher_.k);
                    Segment<Edge<htype>> seg(*pos.edge, pos.pos - kwh.pos, pos.pos + seq.size() - kwh.pos - hasher_.k);
                    return {pos.start, std::vector<Segment<Edge<htype>>>({seg})};
                }
                if (!kwh.hasNext()) {
#pragma omp critical
                    {
                        std::cout << "Error: could not align sequence " << seq.size() << std::endl;
                        std::cout << seq << std::endl;
                        abort();
                    };
                    return {nullptr, std::move(res)};
                }
                kwh = kwh.next();
            }
        }
        Vertex<htype> *prestart = &getVertex(kmers.front());
        if (kmers.front().pos > 0) {
            Vertex<htype> &rcstart = prestart->rc();
            if(!rcstart.hasOutgoing(seq[kmers.front().pos - 1] ^ 3)) {
                std::cout << "No outgoing for start" << std::endl << seq << std::endl <<
                        kmers.front().pos << " " << seq[kmers.front().pos - 1] << std::endl
                        << kmers.front().getSeq() << std::endl;
                VERIFY(false);
            }
            Edge<htype> &rcedge = rcstart.getOutgoing(seq[kmers.front().pos - 1] ^ 3);
            prestart = &rcedge.end()->rc();
            Edge<htype> &edge = rcstart.rcEdge(rcedge);
            VERIFY(edge.size() >= kmers.front().pos);
            Segment<Edge<htype>> seg(edge, edge.size() - kmers.front().pos, edge.size());
            res.emplace_back(seg);
        }
        for(const KWH<htype> & kmer : kmers) {
            if (kmer.pos + hasher_.k < seq.size()) {
                Vertex<htype> &vertex = getVertex(kmer);
                if(!vertex.hasOutgoing(seq[kmer.pos + hasher_.k])) {
                    std::cout << "No outgoing for middle" << std::endl << seq << std::endl <<
                              kmer.pos << " " << size_t(seq[kmer.pos + hasher_.k]) << std::endl
                              << kmer.getSeq() << std::endl;
                    std::cout << vertex.hash() << " " << vertex.outDeg() << " " << vertex.inDeg() << std::endl;
                    for (const Edge<htype> & e : vertex.getOutgoing()) {
                        std::cout << e.seq << std::endl;
                    }
                    VERIFY(false);
                }
                Edge<htype> &edge = vertex.getOutgoing(seq[kmer.pos + hasher_.k]);
                Segment<Edge<htype>> seg(edge, 0, std::min(seq.size() - kmer.pos - hasher_.k, edge.size()));
                res.emplace_back(seg);
            }
        }
        return {prestart, std::move(res)};
    }

    std::vector<std::pair<const Edge<htype> *, size_t>> carefulAlign(const Sequence & seq) const {
        std::vector<KWH<htype>> kmers = extractVertexPositions(seq);
        std::vector<std::pair<const Edge<htype>*, size_t>> res;
        if(kmers.size() == 0) {
            return res;
        }
        std::vector<const Vertex<htype> *> vertices;
        for(size_t i = 0; i < kmers.size(); i++) {
            vertices.push_back(&getVertex(kmers[i]));
        }
        for(size_t i = 0; i + 1 < kmers.size(); i++) {
            Sequence relevant= seq.Subseq(kmers[i].pos, kmers[i + 1].pos + hasher_.k);
            const Edge<htype> * best = nullptr;
            size_t best_score = 0;
            for(const Edge<htype> &edge : vertices[i]->getOutgoing()) {
                if (edge.end() == vertices[i + 1]) {
                    size_t score = std::min(
                        std::max(   hasher_.k,
                                    edge.common(relevant.Subseq(hasher_.k)) +
                                        vertices[i]->rcEdge(edge).common((!relevant).Subseq(hasher_.k))),
                        edge.size());
                    if(score > best_score) {
                        best_score = score;
                        best = &edge;
                    }
                }
            }
            if(best != nullptr) {
                res.emplace_back(best, best_score);
            }
        }
        return std::move(res);
    }

    void printEdge(std::ostream &os, Vertex<htype> & start, Edge<htype> &edge, bool output_coverage) {
        Vertex<htype> &end = *edge.end();
        os << "\"";
        if (!start.isCanonical())
            os << "-";
        os << start.hash() % 100000 << "\" -> \"";
        if (!end.isCanonical())
            os << "-";
        if (output_coverage)
            os << end.hash() % 100000  << "\" [label=\"" << edge.size() << "(" << edge.getCoverage() << ")\"]\n";
        else
            os << end.hash() % 100000  << "\" [label=\"" << edge.size() << "                                                                                                                                                          \"]\n";
    }
    void printDot(std::ostream &os, bool output_coverage) {
        os << "digraph {\nnodesep = 0.5;\n";
        for(std::pair<const htype, Vertex<htype>> & it : this->v) {
            Vertex<htype> &start = it.second;
            for(Edge<htype> &edge : start.getOutgoing()) {
                Vertex<htype> &end = *edge.end();
                printEdge(os, start, edge, output_coverage);
                if(v.find(end.hash()) == v.end()) {
                    printEdge(os, end.rc(), start.rcEdge(edge), output_coverage);
                }
            }
            for(Edge<htype> &edge : start.rc().getOutgoing()) {
                Vertex<htype> &end = *edge.end();
                printEdge(os, start.rc(), edge, output_coverage);
                if(v.find(end.hash()) == v.end()) {
                    printEdge(os, end.rc(), start.rc().rcEdge(edge), output_coverage);
                }
            }
        }
        os << "}\n";
    }


    void processRead(const Sequence & seq) {
        std::vector<KWH<htype>> kmers = extractVertexPositions(seq);
        if(kmers.size() == 0) {
            std::cout << seq << std::endl;
        }
        VERIFY(kmers.size() > 0);
        std::vector<Vertex<htype> *> vertices;
        for(size_t i = 0; i < kmers.size(); i++) {
            vertices.emplace_back(&getVertex(kmers[i]));
            if(i == 0 || vertices[i] != vertices[i - 1]){
                vertices.back()->setSequence(kmers[i].getSeq());
                vertices.back()->incCoverage();
            }
        }
        for(size_t i = 0; i + 1 < vertices.size(); i++) {
//            TODO: if too memory heavy save only some of the labels
            VERIFY(kmers[i].pos + hasher_.k <= seq.size())
            if (i > 0 && vertices[i] == vertices[i - 1] && vertices[i] == vertices[i + 1] &&
                (kmers[i].pos - kmers[i - 1].pos == kmers[i + 1].pos - kmers[i].pos) &&
                kmers[i + 1].pos - kmers[i].pos < hasher_.k) {
                continue;
            }
            vertices[i]->addEdge(Edge<htype>(vertices[i + 1], Sequence(seq.Subseq(kmers[i].pos + hasher_.k,
                    kmers[i + 1].pos + hasher_.k).str())));
            vertices[i + 1]->rc().addEdge(Edge<htype>(&vertices[i]->rc(),
                    !Sequence(seq.Subseq(kmers[i].pos, kmers[i + 1].pos).str())));
        }
        if (kmers.front().pos > 0) {
            vertices.front()->rc().addEdge(Edge<htype>(nullptr, !(seq.Subseq(0, kmers[0].pos))));
        }
        if (kmers.back().pos + hasher_.k < seq.size()) {
            vertices.back()->addEdge(Edge<htype>(nullptr, seq.Subseq(kmers.back().pos + hasher_.k, seq.size())));
        }
    }

    size_t size() const {
        return v.size();
    }

    void printStats(logging::Logger &logger) {
        std::vector<size_t> arr(10);
        size_t isolated = 0;
        size_t isolatedSize = 0;
        size_t n11 = 0;
        size_t n01 = 0;
        size_t ltips = 0;
        size_t ntips = 0;
        size_t e = 0;
        std::vector<size_t> inout(25);
        for(auto &val: v) {
            Vertex<htype> &tmp = val.second;
            if (tmp.inDeg() == 0 && tmp.outDeg() == 1) {
                std::vector<Edge<htype> *> path = tmp.walkForward(tmp.getOutgoing()[0]);
                if (path.back()->end() != nullptr && path.back()->end()->outDeg() == 0 && path.back()->end()->inDeg() == 1) {
                    isolated += 1;
                    for(Edge<htype> * edge : path) {
                        isolatedSize += edge->size();
                    }
                    isolatedSize += hasher().k;
                }
            }
            if (tmp.inDeg() == 1 && tmp.outDeg() == 0) {
                std::vector<Edge<htype>*> path = tmp.rc().walkForward(tmp.rc().getOutgoing()[0]);
                if (path.back()->end() != nullptr && path.back()->end()->outDeg() == 0 && path.back()->end()->inDeg() == 1) {
                    isolated += 1;
                    for(Edge<htype> * edge : path) {
                        isolatedSize += edge->size();
                    }
                    isolatedSize += hasher().k;
                }
            }
            e == tmp.outDeg() + tmp.inDeg();
            arr[std::min(arr.size() - 1, tmp.outDeg())] += 1;
            arr[std::min(arr.size() - 1, tmp.inDeg())] += 1;
            inout[std::min<size_t>(4u, tmp.outDeg()) * 5 + std::min<size_t>(4u, tmp.inDeg())] += 1;
            inout[std::min<size_t>(4u, tmp.inDeg()) * 5 + std::min<size_t>(4u, tmp.outDeg())] += 1;
            if(tmp.outDeg() == 1 && tmp.inDeg() == 0) {
                Vertex<htype> & tmp1 = tmp.getOutgoing()[0].end()->rc();
                VERIFY(tmp.getOutgoing()[0].end() != nullptr);
            }
            if(tmp.outDeg() == 0 && tmp.inDeg() == 1) {
                Vertex<htype> &tmp1 = tmp.rc().getOutgoing()[0].end()->rc();
                VERIFY(tmp.rc().getOutgoing()[0].end()!= nullptr);
            }
            if (tmp.inDeg() == 1 && tmp.outDeg() == 1) {
                n11 += 1;
            }
            if (tmp.inDeg() + tmp.outDeg() == 1) {
                n01 += 1;
                Edge<htype> tip_edge = tmp.outDeg() == 1 ? tmp.getOutgoing()[0] : tmp.rc().getOutgoing()[0];
                if (tip_edge.end()->outDeg() > 1) {
                    ltips += tip_edge.size();
                    ntips += 1;
                }

            }
            for(const Edge<htype> &edge : tmp.getOutgoing()) {
                e += 1;
            }
            for(const Edge<htype> &edge : tmp.rc().getOutgoing()) {
                e += 1;
            }
        }
        logger << "Graph statistics:" << std::endl;
        logger << "Total edges: " << e / 2 << std::endl;
        logger << "Total vertices: " << v.size() << std::endl;
        logger << "Number of end vertices: " << n01 << std::endl;
        logger << "Number of unbranching vertices: " << n11 << std::endl;
        logger << "Number of isolated edges " << isolated << " " << isolatedSize << std::endl;
        logger << "Distribution of degrees:" << std::endl;
        for(size_t i = 0; i < arr.size(); i++) {
            logger << i << " " << arr[i] << std::endl;
        }
        logger << "Distribution of in/out degrees:" << std::endl;
        for(size_t i = 0; i < inout.size(); i++) {
            logger << inout[i] << " ";
            if(i % 5 == 4)
                logger << std::endl;
        }
    }

    void printCoverageStats(logging::Logger &logger) const {
        std::vector<size_t> cov(100);
        std::vector<size_t> cov_tips(100);
        std::vector<std::vector<size_t>> cov_ldist(100);
        std::vector<std::vector<size_t>> cov_ldist_tips(100);
        for(size_t i = 0; i < cov_ldist.size(); i++) {
            cov_ldist[i].resize(30);
        }
        for(size_t i = 0; i < cov_ldist.size(); i++) {
            cov_ldist_tips[i].resize(30);
        }
        std::vector<size_t> covLen(100);
        std::vector<size_t> covLen_tips(100);
        for(auto &val: v) {
            const Vertex<htype> &tmp = val.second;
            for(const Edge<htype> &edge : tmp.getOutgoing()) {
                cov[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                covLen[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
                cov_ldist[std::min(size_t(edge.getCoverage()), cov.size() - 1)][std::min(edge.size() / 50, cov_ldist[0].size() - 1)] += 1;
                if(edge.getTipSize() < hasher_.k * 2) {
                    cov_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                    covLen_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
                    cov_ldist_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)][std::min(edge.size() / 50, cov_ldist[0].size() - 1)] += 1;
                }
            }
            for(const Edge<htype> &edge : tmp.rc().getOutgoing()) {
                cov[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                covLen[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
                cov_ldist[std::min(size_t(edge.getCoverage()), cov.size() - 1)][std::min(edge.size() / 50, cov_ldist[0].size() - 1)] += 1;
                if(edge.getTipSize() < hasher_.k * 2) {
                    cov_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += 1;
                    covLen_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)] += edge.size();
                    cov_ldist_tips[std::min(size_t(edge.getCoverage()), cov.size() - 1)][std::min(edge.size() / 50, cov_ldist[0].size() - 1)] += 1;
                }
            }
        }
//        logger.info() << "Distribution of coverages:" << std::endl;
//        for(size_t i = 0; i < cov.size(); i++) {
//            logger << i << " " << cov[i] << " " << covLen[i];
//            for(size_t val : cov_ldist[i]) {
//                logger << " " << val;
//            }
//            logger << std::endl;
//            logger << i << " " << cov_tips[i] << " " << covLen_tips[i];
//            for(size_t val : cov_ldist_tips[i]) {
//                logger << " " << val;
//            }
//            logger << std::endl;
//        }
    }

    typename vertex_map_type::iterator begin() {
        return v.begin();
    }

    typename vertex_map_type::iterator end() {
        return v.end();
    }

    typename vertex_map_type::const_iterator begin() const {
        return v.begin();
    }

    typename vertex_map_type::const_iterator end() const {
        return v.end();
    }

    void removeIsolated() {
        vertex_map_type newv;
        for(auto & item : v) {
            if(item.second.outDeg() != 0 || item.second.inDeg() != 0) {
                newv.emplace(item.first, std::move(item.second));
            }
        }
        std::swap(v, newv);
    }

    void removeMarked() {
        vertex_map_type newv;
        for(auto & item : v) {
            if(!item.second.marked()) {
                newv.emplace(item.first, std::move(item.second));
            }
        }
        std::swap(v, newv);
    }

    void printFasta(std::ostream &out) const {
        size_t cnt = 0;
        for(const auto &it : v) {
            const Vertex<htype> &vertex = it.second;
            VERIFY(!vertex.seq.empty());
            for(size_t i = 0; i < vertex.outDeg(); i++) {
                const Edge<htype> & edge = vertex.getOutgoing()[i];
                Sequence tmp = vertex.seq + edge.seq;
                Vertex<htype> &end = *edge.end();
                out << ">" << cnt << "_" << vertex.hash() << int(vertex.isCanonical()) <<
                       "_" << end.hash() << int(end.isCanonical()) << "_" << edge.size() << "_" << edge.getCoverage() << std::endl;
                cnt++;
                out << tmp.str() << "\n";
            }
            const Vertex<htype> &rcvertex = vertex.rc();
            for(size_t i = 0; i < rcvertex.outDeg(); i++) {
                const Edge<htype> & edge = rcvertex.getOutgoing()[i];
                Sequence tmp = rcvertex.seq + edge.seq;
                Vertex<htype> &end = *edge.end();
                out << ">" << cnt << "_" << rcvertex.hash() << int(rcvertex.isCanonical()) <<
                       "_" << end.hash() << int(end.isCanonical()) << "_" << edge.size() << "_" << edge.getCoverage() << std::endl;
                cnt++;
                out << tmp.str() << "\n";
            }
        }
    }

    void printGFA(std::ostream &out, bool calculate_coverage) const {
        out << "H\tVN:Z:1.0" << std::endl;
        size_t cnt = 0;
        for(const auto &it : v) {
            for (const Vertex<htype> *pv : {&it.second, &it.second.rc()}) {
                const Vertex<htype> &vertex = *pv;
                VERIFY(!vertex.seq.empty());
                for(const Edge<htype> & edge : vertex.getOutgoing()) {
                    if(vertex.isCanonical(edge)) {
                        if(calculate_coverage)
                            out << "S\t" << vertex.edgeId(edge) << "\t" << vertex.seq << edge.seq
                                << "\tKC:i:" << edge.intCov() << std::endl;
                        else
                            out << "S\t" << vertex.edgeId(edge) << "\t" << vertex.seq << edge.seq << std::endl;
                    }
                }
            }
        }
        for(const auto &it : v) {
            const Vertex<htype> &vertex =it.second;
            for(const Edge<htype> & out_edge : vertex.getOutgoing()) {
                std::string outid = vertex.edgeId(out_edge);
                bool outsign = vertex.isCanonical(out_edge);
                for(const Edge<htype> & inc_edge : vertex.rc().getOutgoing()) {
                    std::string incid = vertex.rc().edgeId(inc_edge);
                    bool incsign = !vertex.rc().isCanonical(inc_edge);
                    out << "L\t" << incid << "\t" << (incsign ? "+" : "-") << "\t" << outid << "\t"
                        << (outsign ? "+" : "-") << "\t" << hasher_.k << "M" <<std::endl;
                }
            }
        }
    }

    template<class Iterator>
    void fillSparseDBGEdges(Iterator begin, Iterator end, logging::Logger &logger, size_t threads, const size_t min_read_size) {
        typedef typename Iterator::value_type ContigType;
        logger.info() << "Starting to fill edges" << std::endl;
        std::function<void(ContigType &)> task = [this, min_read_size](ContigType & contig) {
            Sequence seq = contig.makeSequence();
            if(seq.size() >= min_read_size)
                processRead(seq);
        };
        processRecords(begin, end, logger, threads, task);
        logger.info() << "Sparse graph edges filled." << std::endl;
    }

    static SparseDBG<htype> loadDBGFromFasta(const io::Library &lib, RollingHash<htype> & hasher, logging::Logger &logger, size_t threads) {
        logger.info() << "Loading graph from fasta" << std::endl;
        io::SeqReader reader(lib);
        ParallelRecordCollector<Sequence> sequences(threads);
        ParallelRecordCollector<htype> vertices(threads);
        std::function<void(StringContig &)> collect_task = [&sequences, &vertices, hasher] (StringContig &contig){
            Sequence seq = contig.makeSequence();
            KWH<htype> start(hasher, seq, 0);
            KWH<htype> end(hasher, !seq, 0);
            vertices.add(start.hash());
            vertices.add(end.hash());
            sequences.add(seq);
        };
        processRecords(reader.begin(), reader.end(), logger, threads, collect_task);
        SparseDBG<htype> res(vertices.begin(), vertices.end(), hasher);
        reader.reset();
        res.fillSparseDBGEdges(reader.begin(), reader.end(), logger, threads, hasher.k + 1);
        logger.info() << "Finished loading graph" << std::endl;
        return std::move(res);
    }
};

template<class htype>
class Component {
private:
    SparseDBG<htype> & graph;
    std::unordered_set<htype, alt_hasher<htype>> v;
    struct EdgeRec {
        Vertex<htype> * start;
        Vertex<htype> * end;
        size_t size;
        size_t cov;
    };

    size_t outDeg(const Vertex<htype> &vert, size_t min_cov) const {
        size_t res = 0;
        for(const Edge<htype> &edge : vert.getOutgoing()) {
            if(edge.getCoverage() >= min_cov) {
                res += 1;
            }
        }
        return res;
    }

    bool isUnbranching(const Vertex<htype> &vert, size_t min_cov) const {
        return v.find(vert.hash()) != v.end() && outDeg(vert, min_cov) == 1 && outDeg(vert.rc(), min_cov) == 1;
    }

    Edge<htype> &getOut(Vertex<htype> &vert, size_t min_cov) {
        for(Edge<htype> &edge : vert.getOutgoing()) {
            if(edge.getCoverage() >= min_cov) {
                return edge;
            }
        }
        VERIFY(false);
    }

    Path<htype> unbranching(Vertex<htype> &vert, Edge<htype> &edge, size_t minCov) {
        std::vector<Edge<htype>*> res;
        res.push_back(&edge);
        Vertex<htype> *cur = edge.end();
        while (cur != &vert && isUnbranching(*cur, minCov)) {
            res.push_back(&getOut(*cur, minCov));
            cur = res.back()->end();
        }
        return Path<htype>(vert, res);
    }

public:
    template<class I>
    Component(SparseDBG<htype> &_graph, I begin, I end) : graph(_graph), v(begin, end) {
    }

    template<class I>
    static Component<htype> neighbourhood(SparseDBG<htype> &graph, I begin, I end, size_t radius, size_t min_coverage = 0) {
        std::unordered_set<htype, alt_hasher<htype>> v;
        std::priority_queue<std::pair<size_t, htype>> queue;
        while(begin != end) {
            queue.emplace(0, *begin);
            ++begin;
        }
        while(!queue.empty()) {
            std::pair<size_t, htype> val = queue.top();
            queue.pop();
            if(v.find(val.second) != v.end())
                continue;
            v.insert(val.second);
            if(val.first > radius)
                continue;
            Vertex<htype> &vert = graph.getVertex(val.second);
            for(Edge<htype> & edge : vert.getOutgoing()) {
                if(edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
            for(Edge<htype> & edge : vert.rc().getOutgoing()) {
                if(edge.getCoverage() >= min_coverage)
                    queue.emplace(val.first + edge.size(), edge.end()->hash());
            }
        }
        return Component<htype>(graph, v.begin(), v.end());
    }

    void printEdge(std::ostream &os, Vertex<htype> & start, Edge<htype> &edge) {
        Vertex<htype> &end = *edge.end();
        os << "\"";
        if (!start.isCanonical())
            os << "-";
        os << start.hash() % 10000 << "\" -> \"";
        if (!end.isCanonical())
            os << "-";
        os << end.hash() % 10000  << "\" [label=\"" << edge.size() << "(" << edge.getCoverage() << ")\"]\n";
    }

    void printEdge(std::ostream &os, Path<htype> & path) {
        size_t len = 0;
        size_t cov = 0;
        for(size_t i = 0; i < path.size(); i++) {
            len += path[i].size();
            cov += path[i].intCov();
        }
        Vertex<htype> &start = path.start();
        Vertex<htype> &end = *path.back().end();
        os << "\"";
        if (!start.isCanonical())
            os << "-";
        os << start.hash() % 10000 << "\" -> \"";
        if (!end.isCanonical())
            os << "-";
        os << end.hash() % 10000  << "\" [label=\"" << len << "(" << double(cov) / len << ")\"]\n";
    }

    size_t size() const {
        return v.size();
    }

    void printDot(std::ostream &os, size_t min_cov = 0) {
        os << "digraph {\nnodesep = 0.5;\n";
        for(htype vid : v) {
            Vertex<htype> &start = graph.getVertex(vid);
            for(Edge<htype> &edge : start.getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Vertex<htype> &end = *edge.end();
                printEdge(os, start, edge);
                if(v.find(end.hash()) == v.end()) {
                    printEdge(os, end.rc(), start.rcEdge(edge));
                }
            }
            for(Edge<htype> &edge : start.rc().getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Vertex<htype> &end = *edge.end();
                printEdge(os, start.rc(), edge);
                if(v.find(end.hash()) == v.end()) {
                    printEdge(os, end.rc(), start.rc().rcEdge(edge));
                }
            }
        }
        os << "}\n";
    }
    void printCompressedDot(std::ostream &os, size_t min_cov = 0) {
        os << "digraph {\nnodesep = 0.5;\n";
        for(htype vid : v) {
            Vertex<htype> &start = graph.getVertex(vid);
            if(isUnbranching(start, min_cov))
                continue;
            for(Edge<htype> &edge : start.getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Path<htype> path = unbranching(start, edge, min_cov);
                printEdge(os, path);
                Vertex<htype> &end = *path.back().end();
                if(v.find(end.hash()) == v.end()) {
                    Path<htype> rcpath = path.RC();
                    printEdge(os, rcpath);
                }
            }
            for(Edge<htype> &edge : start.rc().getOutgoing()) {
                if(edge.getCoverage() < min_cov)
                    continue;
                Path<htype> path = unbranching(start.rc(), edge, min_cov);
                printEdge(os, path);
                Vertex<htype> &end = *path.back().end();
                if(v.find(end.hash()) == v.end()) {
                    Path<htype> rcpath = path.RC();
                    printEdge(os, rcpath);
                }
            }
        }
        os << "}\n";
    }
};

template<typename htype, class Iterator>
void fillCoverage(SparseDBG<htype> &sdbg, logging::Logger &logger, Iterator begin, Iterator end, size_t threads,
                        const RollingHash<htype> &hasher, const size_t min_read_size) {
    typedef typename Iterator::value_type ContigType;
    logger.info() << "Starting to fill edge coverages" << std::endl;
    ParallelRecordCollector<size_t> lens(threads);
    std::function<void(ContigType &)> task = [&sdbg, &lens, min_read_size](ContigType & contig) {
        Sequence seq = std::move(contig.makeSequence());
        if(seq.size() >= min_read_size) {
            GraphAlignment<htype> path = sdbg.align(seq);
            lens.add(path.size());
            for(Segment<Edge<htype>> &seg : path) {
                seg.contig().incCov(seg.size());
            }
            path = sdbg.align(!seq);
            for(Segment<Edge<htype>> &seg : path) {
                seg.contig().incCov(seg.size());
            }
        }
    };
    processRecords(begin, end, logger, threads, task);
    logger.info() << "Edge coverage calculated." << std::endl;
    std::vector<size_t> lens_distr(1000);
    for(size_t l : lens) {
        lens_distr[std::min(l, lens_distr.size() - 1)] += 1;
    }
//    logger.info() << "Distribution of path sizes." << std::endl;
//    for(size_t i = 0; i < lens_distr.size(); i++)
//        std::cout << i << " " << lens_distr[i] << std::endl;
}

template<typename htype>
SparseDBG<htype> constructSparseDBGFromReads(logging::Logger & logger, const io::Library &reads_file, size_t threads, const RollingHash<htype> &hasher,
                                       const std::vector<htype> &hash_list, const size_t w) {
    logger.info() << "Starting construction of sparse de Bruijn graph" << std::endl;
    SparseDBG<htype> sdbg(hash_list.begin(), hash_list.end(), hasher);
    logger.info() << "Vertex map constructed." << std::endl;
    io::SeqReader reader(reads_file, (hasher.k + w) * 20, (hasher.k + w) * 4);
    sdbg.fillSparseDBGEdges(reader.begin(), reader.end(), logger, threads, w + hasher.k - 1);
    return std::move(sdbg);
}


template<typename htype>
void tieTips(logging::Logger &logger, SparseDBG<htype> &sdbg, size_t w, size_t threads) {
    logger.info() << " Collecting tips " << std::endl;
//    TODO reduce memory consumption!! A lot of duplicated k-mer storing
    ParallelRecordCollector<Sequence> old_edges(threads);
    ParallelRecordCollector<htype> new_minimizers(threads);
    std::function<void(std::pair<const htype, Vertex<htype>> &)> task =
            [&sdbg, &old_edges, &new_minimizers](std::pair<const htype, Vertex<htype>> & pair) {
        Vertex<htype> &rec = pair.second;
        VERIFY(!rec.seq.empty());
        for (size_t i = 0; i < rec.getOutgoing().size(); i++) {
            const auto &ext = rec.getOutgoing()[i];
            Sequence seq = rec.seq + ext.seq;
            old_edges.add(seq);
            if (ext.end() == nullptr) {
                KWH<htype> kwh(sdbg.hasher(), seq, ext.size());
                new_minimizers.emplace_back(kwh.hash());
            }
        }
        const Vertex<htype> &rec1 = rec.rc();
        for (size_t i = 0; i < rec1.getOutgoing().size(); i++) {
            const auto &ext = rec1.getOutgoing()[i];
            Sequence seq = rec1.seq + ext.seq;
            old_edges.add(seq);
            if (ext.end() == nullptr) {
                KWH<htype> kwh(sdbg.hasher(), seq, ext.size());
                new_minimizers.emplace_back(kwh.hash());
            }
        }
        rec.clear();
    };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);

//#pragma omp parallel default(none) shared(sdbg, old_edges, new_minimizers, logger)
//    {
//#pragma omp single
//        {
//            for (auto &it: sdbg) {
//                Vertex<htype> &rec = it.second;
//                VERIFY(!rec.seq.empty());
//#pragma omp task default(none) shared(sdbg, old_edges, new_minimizers, rec, logger)
//                {
//                    for (size_t i = 0; i < rec.getOutgoing().size(); i++) {
//                        const auto &ext = rec.getOutgoing()[i];
//                        Sequence seq = rec.seq + ext.seq;
//                        old_edges.add(seq);
//                        if (ext.end() == nullptr) {
//                            KWH<htype> kwh(sdbg.hasher(), seq, ext.size());
//                            new_minimizers.emplace_back(kwh.hash());
//                        }
//                    }
//                    const Vertex<htype> &rec1 = rec.rc();
//                    for (size_t i = 0; i < rec1.getOutgoing().size(); i++) {
//                        const auto &ext = rec1.getOutgoing()[i];
//                        Sequence seq = rec1.seq + ext.seq;
//                        old_edges.add(seq);
//                        if (ext.end() == nullptr) {
//                            KWH<htype> kwh(sdbg.hasher(), seq, ext.size());
//                            new_minimizers.emplace_back(kwh.hash());
//                        }
//                    }
//                    rec.clear();
//                }
//            }
//        }
//    }
    logger.info() << "Added " << new_minimizers.size() << " artificial minimizers from tips." << std::endl;
    logger.info() << "Collected " << old_edges.size() << " old edges." << std::endl;
    for(auto it = new_minimizers.begin(); it != new_minimizers.end(); ++it) {
        sdbg.addVertex(*it);
    }
    logger.info() << "New minimizers added to sparse graph." << std::endl;
    logger.info() << "Refilling graph with edges." << std::endl;
    sdbg.fillSparseDBGEdges(old_edges.begin(), old_edges.end(), logger, threads, sdbg.hasher().k + 1);
    logger.info() << "Finished fixing sparse de Bruijn graph." << std::endl;
}

template<typename htype>
void UpdateVertexTips(Vertex<htype> &rec, ParallelRecordCollector<Vertex<htype> *> &queue) {
    bool ok = true;
    for (const Edge<htype> &edge : rec.getOutgoing()) {
        if (edge.getTipSize() == size_t(-1)) {
            edge.updateTipSize();
        }
        if (edge.getTipSize() == size_t(-1)) {
            ok = false;
        }
    }
    if(ok && rec.inDeg() == 1) {
        queue.add(&(rec.rc().getOutgoing()[0].end()->rc()));
    }
}

template<typename htype>
void findTips(logging::Logger &logger, SparseDBG<htype> &sdbg, size_t threads) {
    logger.info() << " Finding tips " << std::endl;
//    TODO reduce memory consumption!! A lot of duplicated k-mer storing
    ParallelRecordCollector<Vertex<htype> *> queue(threads);
#pragma omp parallel default(none) shared(sdbg, logger, queue)
    {
#pragma omp single
        {
            for (auto &it: sdbg) {
                Vertex<htype> &rec = it.second;
                VERIFY(!rec.seq.empty());
#pragma omp task default(none) shared(sdbg, rec, logger, queue)
                {
                    UpdateVertexTips(rec, queue);
                    UpdateVertexTips(rec.rc(), queue);
                }
            }
        }
    }
    logger.info() << "Found initial tips. Looking for iterative tips" << std::endl;
    size_t cnt = 0;
    while(!queue.empty()) {
        logger.info() << "Iteration " << cnt << ". Queue size " << queue.size() << std::endl;
        std::vector<Vertex<htype> *> prev_queue = queue.collectUnique();
        queue.clear();
#pragma omp parallel default(none) shared(sdbg, logger, prev_queue, queue)
        {
#pragma omp single
            {
                for (auto &it: prev_queue) {
                    Vertex<htype> &rec = *it;
                    VERIFY(!rec.seq.empty());
#pragma omp task default(none) shared(sdbg, rec, logger, queue)
                    {
                        UpdateVertexTips(rec, queue);
                    }
                }
            }
        }
    }
    logger.info() << "Tip finding finished" << std::endl;
}


template<typename htype>
void mergeLoop(Vertex<htype> &start, std::vector<Edge<htype>*> &path) {
    VERIFY(*path.back()->end() == start)
    if(path.size() % 2 == 0 && *path[path.size() / 2]->end() == start.rc()) {
        path =std::vector<Edge<htype>*>(path.begin(), path.begin() + path.size() / 2);
    }
    Sequence newSeq(start.pathSeq(path));
    size_t cov = 0;
    for(const Edge<htype> *e : path) {
        if (e->end()->hash() != start.hash()) {
            e->end()->mark();
            e->end()->rc().mark();
        }
        cov += e->intCov();
    }
    Edge<htype> &new_edge = start.addEdgeLockFree(Edge<htype>(path.back()->end(), newSeq.Subseq(start.seq.size())));
    new_edge.incCov(cov - new_edge.intCov());
    Edge<htype> &rc_new_edge = path.back()->end()->rc().addEdgeLockFree(Edge<htype>(&start.rc(), (!newSeq).Subseq(start.seq.size())));
    rc_new_edge.incCov(cov - rc_new_edge.intCov());
}

template<class htype>
void MergeEdge(SparseDBG<htype> &sdbg, Vertex<htype> &start, Edge<htype> &edge) {
    std::vector<Edge<htype>*> path = start.walkForward(edge);
    Vertex<htype> &end = path.back()->end()->rc();
    if (path.size() > 1 && end.hash() >= start.hash()) {
        VERIFY(start.seq.size() > 0)
        VERIFY(end.seq.size() > 0);
        Sequence newSeq(start.pathSeq(path));
        if (start != end)
            end.lock();
        size_t cov = 0;
        for(size_t i = 0; i + 1 < path.size(); i++) {
//            path[i].end()->clear();
            path[i]->end()->mark();
            path[i]->end()->rc().mark();
            cov += path[i]->intCov();
        }
        cov += path.back()->intCov();
        Edge<htype> &new_edge = start.addEdgeLockFree(Edge<htype>(&end.rc(), newSeq.Subseq(start.seq.size())));
        Edge<htype> &rc_new_edge = end.addEdgeLockFree(Edge<htype>(&start.rc(), (!newSeq).Subseq(start.seq.size())));
        new_edge.incCov(cov - new_edge.intCov());
        rc_new_edge.incCov(cov - rc_new_edge.intCov());
        if (start != end)
            end.unlock();
    }
}

template<class htype>
void mergeLinearPaths(logging::Logger & logger, SparseDBG<htype> &sdbg, size_t threads) {
    logger.info() << "Merging linear unbranching paths" << std::endl;
    std::function<void(std::pair<const htype, Vertex<htype>> &)> task =
            [&sdbg](std::pair<const htype, Vertex<htype>> & pair) {
                Vertex<htype> &start = pair.second;
                if (!start.isJunction())
                    return;
                start.lock();
                for (Edge<htype> &edge: start.getOutgoing()) {
                    MergeEdge(sdbg, start, edge);
                }
                start.unlock();
                start.rc().lock();
                for (Edge<htype> &edge: start.rc().getOutgoing()) {
                    MergeEdge(sdbg, start.rc(), edge);
                }
                start.rc().unlock();
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);

//#pragma omp parallel default(none) shared(sdbg)
//    {
//#pragma omp single
//        {
//            for (auto &it: sdbg) {
//                Vertex<htype> &start = it.second;
//                if (!start.isJunction())
//                    continue;
//#pragma omp task default(none) shared(start, sdbg)
//                {
//                    start.lock();
//                    for (const Edge<htype> &edge: start.getOutgoing()) {
//                        MergeEdge(sdbg, start, edge);
//                    }
//                    start.unlock();
//                    start.rc().lock();
//                    for (const Edge<htype> &edge: start.rc().getOutgoing()) {
//                        MergeEdge(sdbg, start.rc(), edge);
//                    }
//                    start.rc().unlock();
//                }
//            }
//        }
//    }
    logger.info() << "Finished merging linear unbranching paths" << std::endl;
}

template<class htype>
void mergeCyclicPaths(logging::Logger & logger, SparseDBG<htype> &sdbg, size_t threads) {
    logger.info() << "Merging cyclic paths" << std::endl;
    ParallelRecordCollector<htype> loops(threads);
    std::function<void(std::pair<const htype, Vertex<htype>> &)> task =
            [&sdbg, &loops](std::pair<const htype, Vertex<htype>> & pair) {
                Vertex<htype> &start = pair.second;
                if(start.isJunction() || start.marked()) {
                    return;
                }
                std::vector<Edge<htype>*> path = start.walkForward(start.getOutgoing()[0]);
                VERIFY(*path.back()->end() == start);
                bool ismin = true;
                for (const Edge<htype> *e : path) {
                    if (e->end()->hash() < start.hash()) {
                        ismin = false;
                        break;
                    }
                }
                if(ismin) {
                    loops.emplace_back(start.hash());
                }
                start.unlock();
            };
    processObjects(sdbg.begin(), sdbg.end(), logger, threads, task);
    logger.info() << "Found " << loops.size() << " perfect loops" << std::endl;
    for(htype loop : loops) {
        Vertex<htype> &start = sdbg.getVertex(loop);
        std::vector<Edge<htype>*> path = start.walkForward(start.getOutgoing()[0]);
        mergeLoop(start, path);
    }
    logger.info() << "Finished merging cyclic paths" << std::endl;
}

template<class htype>
void mergeAll(logging::Logger & logger, SparseDBG<htype> &sdbg, size_t threads) {
    logger.info() << "Merging unbranching paths" << std::endl;
    mergeLinearPaths(logger, sdbg, threads);
//    sdbg.checkConsistency(threads, logger);
    mergeCyclicPaths(logger, sdbg, threads);
//    sdbg.checkConsistency(threads, logger);
    logger.info() << "Removing isolated vertices" << std::endl;
    sdbg.removeMarked();
    logger.info() << "Finished removing isolated vertices" << std::endl;
//    sdbg.checkConsistency(threads, logger);
}

template<class htype>
void CalculateCoverage(const std::experimental::filesystem::path &dir, const RollingHash<htype> &hasher, const size_t w,
                  const io::Library &lib, size_t threads, logging::Logger &logger, SparseDBG<htype> &dbg) {
    logger.info() << "Calculating edge coverage." << std::endl;
    io::SeqReader reader(lib);
    fillCoverage(dbg, logger, reader.begin(), reader.end(), threads, hasher, w + hasher.k - 1);
    std::ofstream os;
    os.open(dir / "coverages.save");
    os << dbg.size() << std::endl;
    for (std::pair<const htype, Vertex<htype>> &pair : dbg) {
        Vertex<htype> &v = pair.second;
        os << v.hash() << " " << v.outDeg() << " " << v.inDeg() << std::endl;
        for (const Edge<htype> &edge : v.getOutgoing()) {
            os << size_t(edge.seq[0]) << " " << edge.intCov() << std::endl;
        }
        for (const Edge<htype> &edge : v.rc().getOutgoing()) {
            os << size_t(edge.seq[0]) << " " << edge.intCov() << std::endl;
        }
    }
    dbg.printCoverageStats(logger);
    os.close();
}

template<typename htype>
std::experimental::filesystem::path alignLib(logging::Logger &logger, SparseDBG<htype> &dbg, const io::Library &align_lib, const RollingHash<htype> &hasher,
              const size_t w, const std::experimental::filesystem::path &dir, size_t threads) {
    logger.info() << "Aligning reads" << std::endl;
    ParallelRecordCollector<std::string> alignment_results(threads);
    std::string acgt = "ACGT";

    std::function<void(StringContig &)> task = [&dbg, &alignment_results, &hasher, w, acgt](StringContig & contig) {
        Contig read = contig.makeContig();
        if(read.size() < w + hasher.k - 1)
            return;
        Path<htype> path = dbg.align(read.seq).path();
        std::stringstream ss;
        ss << read.id << " " << path.start().hash() << int(path.start().isCanonical()) << " ";
        for (size_t i = 0; i < path.size(); i++) {
            ss << acgt[path[i].seq[0]];
        }
        alignment_results.emplace_back(ss.str());
        Contig rc_read = read.RC();
        Path<htype> rc_path = dbg.align(rc_read.seq).path();
        std::stringstream rc_ss;
        rc_ss << rc_read.id << " " << rc_path.start().hash() << int(rc_path.start().isCanonical()) << " ";
        for (size_t i = 0; i < rc_path.size(); i++) {
            rc_ss << acgt[rc_path[i].seq[0]];
        }
        alignment_results.emplace_back(rc_ss.str());
    };
    std::experimental::filesystem::path alignments_file = dir / "alignments.txt";
    std::ofstream os(alignments_file);
    io::SeqReader reader(align_lib);
    processRecords(reader.begin(), reader.end(), logger, threads, task);
    for(std::string & rec : alignment_results) {
        os << rec << "\n";
    }
    os.close();
    logger.info() << "Finished read alignment. Results are in " << (dir / "alignments.txt") << std::endl;
    return alignments_file;
}