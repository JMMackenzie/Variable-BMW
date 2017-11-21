/* This binary generates a ds2i query file from a wand query file of the form
 * qid;term1 term2 term3 .. termk
 */

#include <unordered_map>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>

void usage() {
  std::cerr << "Usage: ./this <dict_file> <query_file>\n";
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {

  if (argc != 3) {
    usage();
  }

  std::ifstream in_dict (argv[1]);
  std::ifstream in_query (argv[2]);

  std::unordered_map<std::string, uint32_t> term_map;
  std::string in_term;
  uint32_t id, y, z;
  while (in_dict >> in_term >> id >> y >> z) {
    term_map[in_term] = id;
    //std::cerr << in_term << " " << id << " " << y << " " << z << "\n";
  }
  std::cerr << "Read: " << term_map.size() << " entries." << std::endl;

  std::string query_str;
  while (std::getline(in_query,query_str)) {       
    //std::cerr << query_str << std::endl; 
    auto id_sep_pos = query_str.find(';');
    auto qryid_str = query_str.substr(0,id_sep_pos);
    auto qry_id = std::stoull(qryid_str);
    std::cerr << "Processing query ID = " << qry_id << std::endl;
    auto qry_content = query_str.substr(id_sep_pos+1);
    std::istringstream qry_content_stream(qry_content);
    std::cout << qry_id;
    for(std::string qry_token; std::getline(qry_content_stream,qry_token,' ');) {
      auto id_itr = term_map.find(qry_token);
      if(id_itr != term_map.end()) {
        std::cout << " " << id_itr->second; 
      }
      else {
        std::cerr << "Error: Could not find term '" << qry_token 
                  << "' in the lexicon.\n";
      }
    }
    std::cout << std::endl;
  }
}

