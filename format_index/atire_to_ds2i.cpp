#include <iostream>
#include <chrono>
#include <vector>
#include <unordered_map>
#include <map>
#include <fstream>
#include <utility>
#include <algorithm>

// ATIRE Includes are required in makefile
#include "ant_param_block.h"
#include "search_engine.h"
#include "search_engine_btree_leaf.h"
#include "btree_iterator.h"
#include "memory.h"

#include <string>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>

char *ATIRE_DOCUMENT_FILE_START = "~documentfilenamesstart";
char *ATIRE_DOCUMENT_FILE_END = "~documentfilenamesfinish";

int main(int argc, char **argv)
{
	ANT_ANT_param_block params(argc, argv);
	long last_param = params.parse();

	if (last_param == argc)
	{
		std::cout << "USAGE: " << argv[0];
		std::cout << " -findex <atire_index> output_basename\n";
		return EXIT_FAILURE;
	}
	using clock = std::chrono::high_resolution_clock;

  std::string basename = argv[argc-1];

	auto build_start = clock::now();

  std::ofstream doc_file(basename+".docs", std::ios::binary);
  std::ofstream freq_file(basename+".freqs", std::ios::binary);
  std::ofstream size_file(basename+".sizes", std::ios::binary);
  std::ofstream doc_id_file(basename+".docids");
  std::ofstream dict_file(basename+".lexicon");

  // load stuff
  ANT_memory memory;
  ANT_search_engine search_engine(&memory);
  search_engine.open(params.index_filename);

  // Write first sequence: no docs in collection
  uint32_t t = 1;
  doc_file.write((char *)&t, sizeof(uint32_t));
  t = search_engine.document_count();
  doc_file.write((char *)&t, sizeof(uint32_t));

  // Also write no. docs in as the sequence header for the doc lengths file
  size_file.write((char *)&t, sizeof(uint32_t));
 
  std::cout << "Writing document lengths and ID's." << std::endl; 
  // write the lengths and names
  {
    
    long long start = search_engine.get_variable(ATIRE_DOCUMENT_FILE_START);
    long long end = search_engine.get_variable(ATIRE_DOCUMENT_FILE_END);
    unsigned long bsize = end - start;
    char *buffer = (char *)malloc(bsize);
    auto filenames = search_engine.get_document_filenames(buffer, &bsize);

    uint64_t uniq_terms = search_engine.get_unique_term_count();
    double mean_length;
    auto lengths = search_engine.get_document_lengths(&mean_length);
    {
      for (long long i = 0; i < search_engine.document_count(); i++)
      {
        uint32_t length = lengths[i];
        size_file.write((char *)&length, sizeof(uint32_t));
        doc_id_file << filenames[i] << std::endl;
      }
    }

    free(buffer);
  }
  // write dictionary
  {
    std::cout << "Writing dictionary. " << std::endl;

    ANT_search_engine_btree_leaf leaf;
    ANT_btree_iterator iter(&search_engine);

    size_t j = 0;
    for (char *term = iter.first(NULL); term != NULL; term = iter.next()) {
      iter.get_postings_details(&leaf);
      dict_file << term << " " << j << " "
        << leaf.local_document_frequency << " "
        << leaf.local_collection_frequency << " "
        << "\n";
      j++;
    }
  }

  // write inverted files
  {
    ANT_search_engine_btree_leaf leaf;
    ANT_btree_iterator iter(&search_engine);
    ANT_impact_header impact_header;
    ANT_compression_factory factory;

    ANT_compressable_integer *raw;
    long long impact_header_size = ANT_impact_header::NUM_OF_QUANTUMS * sizeof(ANT_compressable_integer) * 3;
    ANT_compressable_integer *impact_header_buffer = (ANT_compressable_integer *)malloc(impact_header_size);
    auto postings_list_size = search_engine.get_postings_buffer_length();
    auto raw_list_size = sizeof(*raw) * (search_engine.document_count() + ANT_COMPRESSION_FACTORY_END_PADDING);
    unsigned char *postings_list = (unsigned char *)malloc((size_t)postings_list_size);
    raw = (ANT_compressable_integer *)malloc((size_t)raw_list_size);

    std::vector<std::pair<uint32_t, uint32_t>> post; 
    std::cout << "Writing postings lists." << std::endl;
    uint64_t term_count = 0;

     for (char *term = iter.first(NULL); term != NULL; term_count++, term = iter.next())
    {
	// don't capture ~ terms, they are specific to ATIRE
      if (*term == '~')
        break;

      iter.get_postings_details(&leaf);
      postings_list = search_engine.get_postings(&leaf, postings_list);

      auto the_quantum_count = ANT_impact_header::get_quantum_count(postings_list);
      auto beginning_of_the_postings = ANT_impact_header::get_beginning_of_the_postings(postings_list);
      factory.decompress(impact_header_buffer, postings_list + ANT_impact_header::INFO_SIZE, the_quantum_count * 3);

      if (term_count % 100000 == 0) {
      /* if (true) { */
        std::cout << term << " @ " << leaf.postings_position_on_disk << " (cf:" << leaf.local_collection_frequency << ", df:" << leaf.local_document_frequency << ", q:" << the_quantum_count << ")" << std::endl;
		fflush(stdout);
      }

      long long docid, max_docid, sum;
      ANT_compressable_integer *impact_header = (ANT_compressable_integer *)impact_header_buffer;
      ANT_compressable_integer *current, *end;

      max_docid = sum = 0;
      ANT_compressable_integer *impact_value_ptr = impact_header;
      ANT_compressable_integer *doc_count_ptr = impact_header + the_quantum_count;
      ANT_compressable_integer *impact_offset_start = impact_header + the_quantum_count * 2;
      ANT_compressable_integer *impact_offset_ptr = impact_offset_start;

      post.clear();
      post.reserve(leaf.local_document_frequency);


      while (doc_count_ptr < impact_offset_start) {
        factory.decompress(raw, postings_list + beginning_of_the_postings + *impact_offset_ptr, *doc_count_ptr);
        docid = -1;
        current = raw;
        end = raw + *doc_count_ptr;
        while (current < end) {
          docid += *current++;
          post.emplace_back(docid, *impact_value_ptr);
        }
        impact_value_ptr++;
        impact_offset_ptr++;
        doc_count_ptr++;
      }

      // The above will result in sorted by impact first, so re-sort by docid
      std::sort(std::begin(post), std::end(post));

      uint32_t list_length = post.size();
      doc_file.write((char *)&list_length, sizeof(uint32_t));
      freq_file.write((char *)&list_length, sizeof(uint32_t));

      // Dump to file now
      for (size_t cnt = 0; cnt < post.size(); ++cnt) {
        uint32_t doc_id = post[cnt].first;
        uint32_t term_freq = post[cnt].second;
        doc_file.write((char *)&doc_id, sizeof(uint32_t));
        freq_file.write((char *)&term_freq, sizeof(uint32_t));
      }

    }
  }
  doc_file.close();
  freq_file.close();
  size_file.close();
  doc_id_file.close();
  dict_file.close();
 
	auto build_stop = clock::now();
	auto build_time_sec = std::chrono::duration_cast<std::chrono::seconds>(build_stop-build_start);
	std::cout << "Index transformed in " << build_time_sec.count() << " seconds." << std::endl;

	return EXIT_SUCCESS;
}
