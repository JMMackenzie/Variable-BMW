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
            << " index_type query_algorithm index_filename --map map_filename [--output out_name] [--wand wand_data_filename]"
            << " [--compressed-wand] [--query query_filename]" << std::endl;
}
} // namespace

using namespace ds2i;

template<typename Functor>
void op_dump_trec(Functor query_func, // XXX!!!
                 std::vector<ds2i::term_id_vec> const &queries,
                 std::vector<std::string>& id_map,
                 std::ofstream& output) {
    using namespace ds2i;
    
    // Run queries
    size_t arbitrary_id = 1;
    for (auto const &query: queries) {
      std::vector<std::pair<float, uint64_t>> top_k;
      top_k = query_func(query);
      for (size_t n = 0; n < top_k.size(); ++n) {
        output << arbitrary_id << " "
               << "Q0" << " "
               << id_map[top_k[n].second] << " "
               << n+1 << " "
               << top_k[n].first << " "
               << "VBMW" << std::endl; 
      }
      ++arbitrary_id;
    }

}

template<typename IndexType, typename WandType>
void effectivenesstest(const char *index_filename,
              const char *wand_data_filename,
              std::vector<ds2i::term_id_vec> const &queries,
              std::string const &type,
              std::string const &query_type,
              const char *map_filename,
              const char *output_filename) {
    using namespace ds2i;
    IndexType index;
    logger() << "Loading index from " << index_filename << std::endl;
    boost::iostreams::mapped_file_source m(index_filename);
    succinct::mapper::map(index, m);

    std::vector<std::string> doc_map;
    logger() << "Loading map file from " << map_filename << std::endl;
    std::ifstream map_in(map_filename);
    std::string t_docid;
    while (map_in >> t_docid) {
      doc_map.emplace_back(t_docid); 
    }
    logger() << "Loaded " << doc_map.size() << " DocID's" << std::endl; 
  
    logger() << "Warming up posting lists" << std::endl;
    std::unordered_set<term_id_type> warmed_up;
    for (auto const &q: queries) {
        for (auto t: q) {
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

    std::ofstream output_handle(output_filename);
 
    uint64_t k = configuration::get().k;
    logger() << "Performing " << type << " queries" << std::endl;
    for (auto const &t: query_types) {
        logger() << "Query type: " << t << std::endl;
       //  std::function<uint64_t(ds2i::term_id_vec)> query_fun;
        
std::function<std::vector<std::pair<float, uint64_t>>(ds2i::term_id_vec)> query_fun;
        if (t == "wand" && wand_data_filename) {
            query_fun = [&](ds2i::term_id_vec query) {
              // Default returns count of top-k, but we want the vector
              auto tmp = wand_query<WandType>(wdata, k);
              tmp(index, query); 
              return tmp.topk();
            };
        } else if (t == "block_max_wand" && wand_data_filename) {
            query_fun = [&](ds2i::term_id_vec query) {
              auto tmp = block_max_wand_query<WandType>(wdata, k);
              tmp(index, query);
              return tmp.topk();
            };
        } else if (t == "ranked_or" && wand_data_filename) {
            query_fun = [&](ds2i::term_id_vec query) { 
              auto tmp = ranked_or_query<WandType>(wdata, k);
              tmp(index, query);
              return tmp.topk();
          };
        } else if (t == "maxscore" && wand_data_filename) {
            query_fun = [&](ds2i::term_id_vec query) { 
              auto tmp = maxscore_query<WandType>(wdata, k);
              tmp(index, query); 
              return tmp.topk();
            };
        } else {
            logger() << "Unsupported query type: " << t << std::endl;
            break;
        }
        op_dump_trec(query_fun, queries, doc_map, output_handle);
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
    const char *map_filename = nullptr;
    const char *out_filename = nullptr;
    bool compressed = false;
    std::vector<term_id_vec> queries;

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

        if (arg == "--map") {
            map_filename = argv[++i]; 
        }

        if (arg == "--output") {
          out_filename = argv[++i];
        }
    }

    if (out_filename == nullptr || map_filename == nullptr) {
      std::cerr << "ERROR: Must provide map and output file. Quitting."
                << std::endl;
      return EXIT_FAILURE;
    }

    term_id_vec q;
    if(query_filename){
        std::filebuf fb;
        if (fb.open(query_filename, std::ios::in)) {
            std::istream is(&fb);
            while (read_query(q, is)) queries.push_back(q);
        }
    } else {
        while (read_query(q)) queries.push_back(q);
    }

    /**/
    if (false) {
#define LOOP_BODY(R, DATA, T)                                                       \
        } else if (type == BOOST_PP_STRINGIZE(T)) {                                 \
            if (compressed) {                                                       \
                 effectivenesstest<BOOST_PP_CAT(T, _index), wand_uniform_index>              \
                 (index_filename, wand_data_filename, queries, type, query_type, map_filename, out_filename);   \
            } else {                                                                \
                effectivenesstest<BOOST_PP_CAT(T, _index), wand_raw_index>                   \
                (index_filename, wand_data_filename, queries, type, query_type, map_filename, out_filename);    \
            }                                                                       \
    /**/

BOOST_PP_SEQ_FOR_EACH(LOOP_BODY, _, DS2I_INDEX_TYPES);
#undef LOOP_BODY

    } else {
        logger() << "ERROR: Unknown type " << type << std::endl;
    }

}
