#include <iostream>
#include <unordered_map>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include <succinct/mapper.hpp>


#include "index_types.hpp"
#include "wand_data_compressed.hpp"
#include "wand_data_raw.hpp"
#include "queries.hpp"
#include "util.hpp"

namespace {
void printUsage(const std::string &programName) {
  std::cerr << "Usage: " << programName
            << " index_type query_algorithm index_filename [--wand wand_data_filename]"
            << " [--compressed-wand] [--query query_filename]" << std::endl;
}
} // namespace

using namespace ds2i;

template<typename Functor>
void op_perftest(Functor query_func, // XXX!!!
                 std::vector<std::pair<uint32_t, ds2i::term_id_vec>> const &queries,
                 size_t runs) {
    using namespace ds2i;

    std::map<uint32_t, double> query_times;

    for (size_t run = 0; run <= runs; ++run) {
        for (auto const &query: queries) {
            
            auto tick = get_time_usecs();
            uint64_t result = query_func(query.second);
            do_not_optimize_away(result);
            double elapsed = double(get_time_usecs() - tick);
            if (run != 0) { // first run is not timed
                auto itr = query_times.find(query.first);
                if(itr != query_times.end()) {
                    itr->second += elapsed;
                } else {
                    query_times[query.first] = elapsed;
                }
            }
        }
    }

    // Take mean of the timings and dump per-query
    for(auto& timing : query_times) {
      timing.second = timing.second / (runs-1);
      std::cout << timing.first << ";" << (timing.second / 1000.0) << std::endl;
    }

}

template<typename IndexType, typename WandType>
void perftest(const char *index_filename,
              const char *wand_data_filename,
              std::vector<std::pair<uint32_t, ds2i::term_id_vec>> const &queries,
              std::string const &type,
              std::string const &query_type,
              const uint64_t m_k = 0) {
    using namespace ds2i;
    IndexType index;
    logger() << "Loading index from " << index_filename << std::endl;
    boost::iostreams::mapped_file_source m(index_filename);
    succinct::mapper::map(index, m);



    logger() << "Warming up posting lists" << std::endl;
    std::unordered_set<term_id_type> warmed_up;
    for (auto const &q: queries) {
        for (auto t: q.second) {
            if (!warmed_up.count(t)) {
                index.warmup(t);
                warmed_up.insert(t);
            }
        }
    }

    WandType wdata;

    std::vector<std::string> query_types;
    boost::algorithm::split(query_types, query_type, boost::is_any_of(":"));
    boost::iostreams::mapped_file_source md;
    if (wand_data_filename) {
        md.open(wand_data_filename);
        succinct::mapper::map(wdata, md, succinct::mapper::map_flags::warmup);
    }

    uint64_t k = m_k;
    if (k == 0) {
      k = configuration::get().k;
    }
    
    logger() << "Performing " << type << " queries" << std::endl;
    for (auto const &t: query_types) {
        logger() << "Query type: " << t << std::endl;
        std::function<uint64_t(ds2i::term_id_vec)> query_fun;
        if (t == "and") {
            query_fun = [&](ds2i::term_id_vec query) { return and_query<false>()(index, query); };
        } else if (t == "and_freq") {
            query_fun = [&](ds2i::term_id_vec query) { return and_query<true>()(index, query); };
        } else if (t == "or") {
            query_fun = [&](ds2i::term_id_vec query) { return or_query<false>()(index, query); };
        } else if (t == "or_freq") {
            query_fun = [&](ds2i::term_id_vec query) { return or_query<true>()(index, query); };
        } else if (t == "wand" && wand_data_filename) {
            query_fun = [&](ds2i::term_id_vec query) { return wand_query<WandType>(wdata, k)(index, query); };
        } else if (t == "block_max_wand" && wand_data_filename) {
            query_fun = [&](ds2i::term_id_vec query) { return block_max_wand_query<WandType>(wdata, k)(index, query); };
        } else if (t == "ranked_or" && wand_data_filename) {
            query_fun = [&](ds2i::term_id_vec query) { return ranked_or_query<WandType>(wdata, k)(index, query); };
        } else if (t == "maxscore" && wand_data_filename) {
            query_fun = [&](ds2i::term_id_vec query) { return maxscore_query<WandType>(wdata, k)(index, query); };
        } else {
            logger() << "Unsupported query type: " << t << std::endl;
            break;
        }
        op_perftest(query_fun, queries, 4);
    }


}

typedef wand_data<bm25, wand_data_raw<bm25>> wand_raw_index;
typedef wand_data<bm25, wand_data_compressed<bm25, uniform_score_compressor>> wand_uniform_index;

int main(int argc, const char **argv) {
    using namespace ds2i;

    std::string programName = argv[0];
    if (argc < 3) {
    printUsage(programName);
    return 1;
    }

    std::string type = argv[1];
    std::string query_type = argv[2];
    const char *index_filename = argv[3];
    const char *wand_data_filename = nullptr;
    const char *query_filename = nullptr;
    uint64_t m_k = 0;
    bool compressed = false;
    std::vector<std::pair<uint32_t, term_id_vec>> queries;

    for (int i = 4; i < argc; ++i) {
        std::string arg = argv[i];

        if(arg == "--wand"){
            wand_data_filename = argv[++i];;
        }

        if(arg == "--compressed-wand"){
            compressed = true;
        }

        if (arg == "--query") {
            query_filename = argv[++i];
        }

        if (arg == "--k") {
          m_k = std::stoull(argv[++i]);
        }
    }

    term_id_vec q;
    uint32_t qid;
    if(query_filename){
        std::filebuf fb;
        if (fb.open(query_filename, std::ios::in)) {
            std::istream is(&fb);
            while (read_query(q, qid, is)) queries.emplace_back(qid, q);
        }
    } else {
        while (read_query(q, qid)) queries.emplace_back(qid, q);
    }

    /**/
    if (false) {
#define LOOP_BODY(R, DATA, T)                                                            \
        } else if (type == BOOST_PP_STRINGIZE(T)) {                                      \
            if (compressed) {                                                            \
                 perftest<BOOST_PP_CAT(T, _index), wand_uniform_index>                   \
                 (index_filename, wand_data_filename, queries, type, query_type, m_k);   \
            } else {                                                                     \
                perftest<BOOST_PP_CAT(T, _index), wand_raw_index>                        \
                (index_filename, wand_data_filename, queries, type, query_type, m_k);    \
            }                                                                            \
    /**/

BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, DS2I_INDEX_TYPES);
#undef LOOP_BODY

    } else {
        logger() << "ERROR: Unknown type " << type << std::endl;
    }

}
