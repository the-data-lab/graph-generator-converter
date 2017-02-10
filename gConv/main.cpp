#include "csr_map.h"
#include "csr_vector.h"
#include "csr_array.h"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>


int main(int argc, char* argv[]) {


	index_t *strt_pos, *adj_card;
	data_t *adj_list;
	index_t edge_count, vert_count;

	csr_array<data_t, index_t, g500_t>
	(argc, argv, adj_list, strt_pos, 
	adj_card, edge_count, vert_count);

	return 0;
}
