//
// Created by anton on 17.07.2020.
//
#define _GLIBCXX_PARALLEL
#include "dbg_construction.hpp"
#include "dbg_disjointigs.hpp"
#include "minimizer_selection.hpp"
#include "sparse_dbg.hpp"
#include "rolling_hash.hpp"
#include "hash_utils.hpp"
#include "sequences/seqio.hpp"
#include "common/dir_utils.hpp"
#include "common/cl_parser.hpp"
#include "common/logging.hpp"
#include "hash_utils.hpp"
#include <iostream>
#include <queue>
#include <omp.h>
#include <unordered_set>
#include <wait.h>

using logging::Logger;

//
//class htype128 {
//    unsigned __int_128 val;
//public:
//    template<class T>
//    htype128(cosnt T &_val) : val(_val) {
//    }
//
//    htype128 operator+(const htype128& other) const {
//        return htype128(val + other.val);
//    }
//
//    htype128 operator+(const htype128& other) const {
//        return htype128(val + other.val);
//    }
//
//    htype128 operator*(const htype128& other) const {
//        return htype128(val * other.val);
//    }
//
//    htype128 operator/(const htype128& other) const {
//        return htype128(val / other.val);
//    }
//
//    bool operator<(const htype128& other) const {
//        return hval < other.val;
//    }
//    bool operator>(const htype128& other) const {
//        return hval > other.val;
//    }
//    bool operator<=(const htype128& other) const {
//        return hval <= other.val;
//    }
//    bool operator>=(const htype128& other) const {
//        return hval >= other.val;
//    }
//    bool operator==(const htype128& other) const {
//        return hval == other.val;
//    }
//    bool operator!=(const htype128& other) const {
//        return hval != other.val;
//    }
//};
//
//std::ostream &operator<<(std::ostream &os, htype128 val) {
//    while(val > 0) {
//        os << size_t(val % 10);
//        val /= 10;
//    }
//    return os;
//}

struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const
    {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

void analyseGenome(SparseDBG<htype128> &dbg, const std::string &ref_file, const std::experimental::filesystem::path &path_dump,
                   const std::experimental::filesystem::path &mult_dump, Logger &logger) {
    logger.info() << "Reading reference" << std::endl;
    std::vector<StringContig> ref = io::SeqReader(ref_file).readAll();
    logger.info() << "Finished reading reference. Starting alignment" << std::endl;
    std::vector<Segment<Edge<htype128>>> path;
    std::ofstream os;
    os.open(path_dump);
    size_t cur = 0;
    for(StringContig & contig : ref) {
        auto tmp = dbg.align(contig.makeSequence());
        os << "New chromosome " << contig.id << std::endl;
        for(size_t i = 0; i < tmp.size(); i++) {
            const Segment<Edge<htype128>> &seg = tmp[i];
            os << "[" << cur << ", " << cur + seg.size() << "] -> [" << seg.left << ", " << seg.right <<"] ";
            os << tmp.getVertex(i).hash() << tmp.getVertex(i).isCanonical() << " "
               << tmp.getVertex(i + 1).hash() << tmp.getVertex(i + 1).isCanonical() << " "
               << seg.size() << " " << seg.contig().getCoverage() << std::endl;
            cur += seg.size();
        }
        logger.info() << "Aligned chromosome " << contig.id << " . Path length " << tmp.size() << std::endl;
        path.insert(path.end(), tmp.begin(), tmp.end());
    }
    os.close();
    logger.info() << "Reference path consists of " << path.size() << " edges" << std::endl;
//    std::vector<size_t> fill_distr(11);
//    std::vector<size_t> fill_distr_len(11);
//    for(const std::pair<const Edge<htype128> *, size_t> & val : path) {
//        size_t ind = val.second * 10 / val.first->size();
//        fill_distr[ind] += 1;
//        fill_distr_len[ind] += 1;
//    }
//    logger.info() << "Edge filling distribution" <<std::endl;
//    logger << fill_distr << std::endl << fill_distr_len << std::endl;
    size_t max_cov = 50;
    std::vector<size_t> cov(max_cov + 1);
    std::vector<size_t> cov_len(max_cov + 1);
    std::vector<size_t> cov_bad(max_cov + 1);
    std::vector<size_t> cov_bad_len(max_cov + 1);
    std::vector<size_t> cov_good(max_cov + 1);
    std::vector<size_t> cov_good_len(max_cov + 1);
    std::unordered_map<Edge<htype128> const *, size_t> eset;
    for(size_t i = 0; i < path.size(); i++) {
        const Edge<htype128> &edge = path[i].contig();
        eset[&edge] += 1;
    }
    std::ofstream os_mult;
    os_mult.open(mult_dump);
    for(auto & it : eset) {
        os_mult << it.second << " " << it.first->getCoverage() << " " << it.first->size() << std::endl;
    }
    os_mult.close();
    for(auto & pair : dbg) {
        Vertex<htype128> &vert = pair.second;
        for (Edge<htype128> &edge : vert.getOutgoing()) {
            size_t cov_val = std::min(max_cov, size_t(edge.getCoverage()));
            if (eset.find(&edge) == eset.end() && eset.find(&vert.rcEdge(edge)) == eset.end()) {
                cov_bad[cov_val] += 1;
                cov_bad_len[cov_val] += edge.size();
            } else {
                cov_good[cov_val] += 1;
                cov_good_len[cov_val] += edge.size();
            }
            cov[cov_val] += 1;
            cov_len[cov_val] += edge.size();
        }
    }
    logger.info() << "All coverages" << std::endl;
    logger << cov << std::endl << cov_len << std::endl;
    logger.info() << "Coverages of edges in genome path" << std::endl;
    logger << cov_good << std::endl << cov_good_len << std::endl;
    logger.info() << "Coverages of edges outside genome path" << std::endl;
    logger << cov_bad << std::endl << cov_bad_len << std::endl;
}

void LoadCoverage(const std::experimental::filesystem::path &fname, Logger &logger, SparseDBG<htype128> &dbg) {
    logger.info() << "Loading edge coverages." << std::endl;
    std::ifstream is;
    is.open(fname);
    size_t n;
    is >> n;
    for (size_t i = 0; i < n; i++) {
        htype128 vid;
        is >> vid;
        Vertex<htype128> *v = &dbg.getVertex(vid);
        size_t inDeg, outDeg;
        is >> outDeg >> inDeg;
        for (size_t j = 0; j < inDeg + outDeg; j++) {
            if (j == outDeg)
                v = &v->rc();
            size_t next;
            is >> next;
            Edge<htype128> &edge = v->getOutgoing(char(next));
            size_t cov;
            is >> cov;
            edge.incCov(cov);
        }
    }
    is.close();
    logger.info() << "Finished loading edge coverages." << std::endl;
}
std::string constructMessage() {
    std::stringstream ss;
    ss << "JumboDB version 1.0\n\n";
    ss << "Usage dbg [options] -o <output-dir> -k <int>\n\n";
    ss << "Basic options:\n";
    ss << "  -o <file_name> (or --output-dir <file_name>)  Name of output folder. Resulting graph will be stored there.\n";
    ss << "  -k <int>                                      Value of k (vertex size) to be used for de Bruijn graph construction. k should be odd (otherwise k + 1 is used instead).\n";
    ss << "  --reads <file_name>                           Name of file that contains reads in fasta or fastq format. This option can be used any number of times in the same command line resulting in collecting reads from multiple files.\n";
    ss << "  -h (or --help)                                Print this help message.\n";
    ss << "\nAdvanced options:\n";
    ss << "  -t <int> (or --threads <int>)                 Number of threads. The default value is 16.\n";
    ss << "  -w <int> (or --window <int>`)                 The window size to be used for sparse de Bruijn graph construction. The default value is 2000. Note that all reads of length less than k + w are ignored during graph construction.\n";
    ss << "  --compress                                    Compress all homolopymers in reads.\n";
    ss << "  --coverage                                    Calculate edge coverage of edges in the constructed de Bruijn graph.\n";
    return ss.str();
}
int main(int argc, char **argv) {
    CLParser parser({"vertices=none", "unique=none", "coverages=none", "segments=none", "dbg=none", "output-dir=",
                     "threads=16", "k-mer-size=", "window=2000", "base=239", "debug", "disjointigs=none", "reference=none",
                     "correct", "simplify", "coverage", "cov-threshold=2", "tip-correct", "crude-correct", "initial-correct",
                     "compress", "help"},
                    {"reads", "align"},
                    {"h=help", "o=output-dir", "t=threads", "k=k-mer-size","w=window"},
                    constructMessage());
    parser.parseCL(argc, argv);
    if (parser.getCheck("help")) {
        std::cout << parser.message() << std::endl;
        return 0;
    }
    if (!parser.check().empty()) {
        std::cout << "Incorrect parameters:" << std::endl;
        std::cout << parser.check() << "\n" << std::endl;
        std::cout << parser.message() << std::endl;
        return 1;
    }
    if(parser.getCheck("compress"))
        StringContig::needs_compressing = true;
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    Logger logger;
    logger.addLogFile(ls.newLoggerFile());
    for(size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    size_t k = std::stoi(parser.getValue("k-mer-size"));
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    RollingHash<htype128> hasher(k, std::stoi(parser.getValue("base")));
    const size_t w = std::stoi(parser.getValue("window"));
    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::Library reads_lib = lib;
    if (parser.getValue("reference") != "none") {
        logger.info() << "Added reference to graph construction. Careful, some edges may have coverage 0" << std::endl;
        lib.insert(lib.begin(), std::experimental::filesystem::path(parser.getValue("reference")));
    }
    size_t threads = std::stoi(parser.getValue("threads"));
    omp_set_num_threads(threads);

    std::string disjointigs_file = parser.getValue("disjointigs");
    std::string vertices_file = parser.getValue("vertices");
    std::string dbg_file = parser.getValue("dbg");
    SparseDBG<htype128> dbg = dbg_file == "none" ?
            DBGPipeline(logger, hasher, w, lib, dir, threads, disjointigs_file, vertices_file) :
          SparseDBG<htype128>::loadDBGFromFasta({std::experimental::filesystem::path(dbg_file)}, hasher, logger, threads);

    bool calculate_coverage = parser.getCheck("coverage") || parser.getCheck("simplify") ||
            parser.getCheck("correct") || parser.getValue("segments") != "none" ||
            parser.getValue("reference") != "none" || parser.getCheck("tip-correct") ||
            parser.getCheck("crude-correct") || parser.getCheck("initial-correct");

    if (!parser.getListValue("align").empty() || calculate_coverage) {
        dbg.fillAnchors(w, logger, threads);
    }

    if (calculate_coverage && !parser.getCheck("initial-correct")) {
        if (parser.getValue("coverages") == "none") {
            CalculateCoverage(dir, hasher, w, reads_lib, threads, logger, dbg);
        } else {
            LoadCoverage(parser.getValue("coverages"), logger, dbg);
        }
    }

    if(parser.getValue("dbg") == "none") {
        logger.info() << "Printing graph to fasta file " << (dir / "graph.fasta") << std::endl;
        std::ofstream edges;
        edges.open(dir / "graph.fasta");
        dbg.printFasta(edges);
        edges.close();
        logger.info() << "Printing graph to gfa file " << (dir / "graph.gfa") << std::endl;
        std::ofstream gfa;
        gfa.open(dir / "graph.gfa");
        dbg.printGFA(gfa, calculate_coverage);
        gfa.close();
        logger.info() << "Printing graph to dot file " << (dir / "graph.dot") << std::endl;
        std::ofstream dot;
        dot.open(dir / "graph.dot");
        dbg.printDot(dot, calculate_coverage);
        dot.close();
    }

    if (parser.getCheck("tip-correct")) {
        logger.info() << "Removing tips from reads" << std::endl;
        std::experimental::filesystem::path out = dir / "tip_correct.fasta";
        ParallelRecordCollector<Contig> alignment_results(threads);

        std::function<void(StringContig &)> task = [&dbg, &alignment_results, &hasher, w, &logger](StringContig & contig) {
            Contig read = contig.makeContig();
            if(read.size() < w + hasher.k - 1)
                return;
            GraphAlignment<htype128> gal = dbg.align(read.seq);
            if (gal.size() > 0 && gal.front().contig().getCoverage() < 2 && gal.start().inDeg() == 0 && gal.start().outDeg() == 1) {
                gal = gal.subalignment(1, gal.size());
            }
            if (gal.size() > 0 && gal.back().contig().getCoverage() < 2 && gal.finish().outDeg() == 0 && gal.finish().inDeg() == 1) {
                gal = gal.subalignment(0, gal.size() - 1);
            }
            for(Segment<Edge<htype128>> seg : gal) {
                if (seg.contig().getCoverage() < 2)
                    return;
            }
            if (gal.size() > 0) {
                alignment_results.emplace_back(gal.Seq(), read.id);
            }
        };
        std::ofstream os(dir / "tip_correct.fasta");
        io::SeqReader reader(reads_lib);
        processRecords(reader.begin(), reader.end(), logger, threads, task);
        for(Contig & rec : alignment_results) {
            os << ">" << rec.id << "\n" << rec.seq << "\n";
        }
        os.close();
    }


    if (!parser.getListValue("align").empty()) {
        io::Library align_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("align"));
        alignLib(logger, dbg, align_lib, hasher, w, dir, threads);
    }

    if (parser.getValue("segments") != "none") {
        logger.info() << "Drawing components" << std::endl;
        io::SeqReader segs(parser.getValue("segments"));
        ensure_dir_existance(dir / "pictures");
        size_t cnt = 0;
        for(StringContig scontig : segs) {
            Contig s = scontig.makeContig();
            std::vector<std::pair<const Edge<htype128> *, size_t>> path = dbg.carefulAlign(s.seq);
            std::vector<htype128> hashs;
            hashs.reserve(path.size());
            for(auto & tmp : path) {
                hashs.push_back(tmp.first->end()->hash());
            }
            Component<htype128> comp = Component<htype128>::neighbourhood(dbg, hashs.begin(), hashs.end(), 600, 2);
            std::ofstream os;
            os.open(dir/ "pictures" / (logging::itos(cnt) + ".dot"));
            comp.printCompressedDot(os, 2);
            os.close();
            logger.info() << cnt << " " << comp.size() << std::endl;
            cnt += 1;
        }
    }
    if (parser.getValue("reference") != "none") {
        analyseGenome(dbg, parser.getValue("reference"), dir / "ref.info", dir / "mult.info", logger);
    }
    if (parser.getCheck("simplify")) {
        logger.info() << "Removing low covered edges" << std::endl;
        size_t threshold = std::stoull(parser.getValue("cov-threshold"));
        std::vector<Sequence> edges;
        std::vector<htype128> vertices_again;
        for(auto & it : dbg) {
            Vertex<htype128> &vert = it.second;
            bool add = false;
            for(Edge<htype128> & edge : vert.getOutgoing()) {
                if (edge.getCoverage() >= threshold) {
                    edges.push_back(vert.seq + edge.seq);
                    add = true;
                }
            }
            for(Edge<htype128> & edge : vert.rc().getOutgoing()) {
                if (edge.getCoverage() >= threshold){
                    edges.push_back(vert.rc().seq + edge.seq);
                    add = true;
                }
            }
            if (add)
                vertices_again.push_back(vert.hash());
        }
        SparseDBG<htype128> simp_dbg(vertices_again.begin(), vertices_again.end(), hasher);
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
        std::ofstream dot1;
        dot1.open(dir / "simp_graph1.dot");
        simp_dbg.printDot(dot1, calculate_coverage);
        dot1.close();
        mergeAll(logger, simp_dbg, threads);
        std::ofstream simp_os;
        simp_os.open(dir / "simp_graph.fasta");
        simp_dbg.printFasta(simp_os);
        simp_os.close();
        std::ofstream dot;
        dot.open(dir / "simp_graph.dot");
        simp_dbg.printDot(dot, calculate_coverage);
        dot.close();
    }
    logger.info() << "DBG construction finished" << std::endl;
    return 0;
}
