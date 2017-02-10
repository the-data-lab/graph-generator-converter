
/*
 * Copyright 2016 The George Washington University
 * Written by Pradeep Kumar 
 * Directed by Prof. Howie Huang
 *
 * https://www.seas.gwu.edu/~howie/
 * Contact: iheartgraph@gmail.com
 *
 * 
 * Please cite the following paper:
 * 
 * Pradeep Kumar and H. Howie Huang. 2016. G-Store: High-Performance Graph Store for Trillion-Edge Processing. In Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis (SC '16).
 
 *
 * This file is part of G-Store.
 *
 * G-Store is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * G-Store is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with G-Store.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <omp.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include <algorithm>
#include <errno.h>
#include <cmath>
#include <fstream>
#include <libaio.h>
#include <math.h>
#include <asm/mman.h>
#include "wtime.h"
#include "gConv.h"

void gConv::main(int argc, char * argv[])
{
    int o;
    uint32_t scale;
    int c = 0;
    string edgefile;
    string part_file;
    while ((o = getopt (argc, argv, "s:o:hi:j:c:m:v:a:")) != -1) {
        switch(o) {
            case 's': //scale
                scale = atoi(optarg);
                vert_count = (1L << scale);
                scale -= bit_shift1;
                p = (1 << scale);
                p_s = (vert_count >> bit_shift2) + (0 != (vert_count & part_mask2_2));
                read_part = new bitmap_t(p_s);
                read_part_next = new bitmap_t(p_s);
                read_part_next->reset();
                break;
            case 'v'://vert count
                //vert_count = atoi(optarg);
				sscanf(optarg, "%ld", &vert_count);
                p = (vert_count >> bit_shift1) + (0 != (vert_count & part_mask1_2));
                //p_s = (vert_count >> bit_shift2) + (0 != (vert_count & part_mask2_2));
                p_s = p_p*p;
                read_part = new bitmap_t(p_s);
                read_part_next = new bitmap_t(p_s);
                read_part_next->reset();
                break;
            case 'i':
                edgefile = optarg;
                break;
            case 'o':
                part_file = optarg;
                break;
            case 'h':
                cout << "Coming soon" << endl;
                return;
            case 'c':
                c = atoi(optarg);
                break;
            default:
               cout << "Unknown argument" << endl ; 
        }
    }

    cout << "Partition Count = " << p << endl;
    cout << "Partition Count smaller = " << p_s << endl;
    double start, end;
  
    switch(c) {
	case 1:
		start = mywtime();
		proc_csr(edgefile, part_file);
        save_csr(part_file);
		end = mywtime();
		cout << "CSR Time = " << end - start << endl;
		return;
    default:
        return;
    } 

   
}

void grid::pre_csr(string edgefile, gedge_t* edges, index_t nedges)
{
    _beg_pos = (index_t*) calloc(sizeof(index_t), vert_count + 1);  
	
	#ifndef HALF_GRID
	_beg_pos_in = (index_t*) calloc(sizeof(index_t), vert_count + 1);  
	#endif

   	vert_degree = (degree_t*)calloc(sizeof(degree_t), vert_count);
    
	#pragma omp parallel num_threads(NUM_THDS) 
    {
        gedge_t  edge;
        vertex_t v0, v1;
        
        #pragma omp for
        for(index_t k = 0; k < nedges; ++k) {
            edge = edges[k];
            if (edge.is_self_loop()) continue;
	    
            v0 = edge.get_v0();
            v1 = edge.get_v1();
            
            __sync_fetch_and_add(_beg_pos + v0, 1);
			__sync_fetch_and_add(vert_degree + v0, 1);
			
			#ifdef HALF_GRID
            __sync_fetch_and_add(_beg_pos + v1, 1);
			__sync_fetch_and_add(vert_degree + v1, 1);
			#else 
            __sync_fetch_and_add(_beg_pos_in + v1, 1);
			#endif
        }
    }
    
    //Calculate the CSR beg_pos
    index_t prefix_sum = 0;
	index_t curr_value = 0;
	#ifndef HALF_GRID
    index_t prefix_sum_in = 0;
	index_t curr_value_in = 0;
	#endif
    for (index_t ipart = 0; ipart < vert_count; ++ipart) {
		curr_value = _beg_pos[ipart];
        _beg_pos[ipart] = prefix_sum;
        prefix_sum += curr_value;
		#ifndef HALF_GRID
		curr_value_in = _beg_pos_in[ipart];
        _beg_pos_in[ipart] = prefix_sum_in;
        prefix_sum_in += curr_value_in;
		#endif
    }

    _beg_pos[vert_count] = prefix_sum;
	#ifndef HALF_GRID
    _beg_pos_in[vert_count] = prefix_sum_in;
    cout << "Total edges = " << prefix_sum_in << endl;
	#endif
    cout << "Total edges = " << prefix_sum << endl;

}
void grid::proc_csr(string edgefile, string part_file)
{
    //read the binary edge file
    int fid_edge = open(edgefile.c_str(), O_RDONLY);
    struct stat st_edge;
    fstat(fid_edge, &st_edge);
    assert(st_edge.st_size != 0);
    index_t nedges = st_edge.st_size/sizeof(gedge_t);
    gedge_t* edges;
    /*
    edges = (gedge_t*)mmap(0, st_edge.st_size, PROT_READ, 
                                    MAP_PRIVATE, fid_edge, 0);
    madvise(edges, st_edge.st_size, MADV_SEQUENTIAL);
    */
    edges = (gedge_t*) malloc(st_edge.st_size);
    FILE* f = fopen(edgefile.c_str(), "rb");
    fread(edges, sizeof(gedge_t), nedges, f);

    double start = mywtime();

    pre_csr(edgefile, edges, nedges);

    index_t* count_edge = (index_t*) calloc(sizeof(index_t), vert_count);  
    _adj = (uint32_t*)calloc(sizeof(uint32_t), _beg_pos[vert_count]);
    
	#ifndef HALF_GRID
	uint32_t* _adj_in = (uint32_t*)calloc(sizeof(uint32_t), _beg_pos_in[vert_count]);
    index_t* count_edge_in = (index_t*) calloc(sizeof(index_t), vert_count);  
	#endif

    //---classify the edges in the grid
    #pragma omp parallel num_threads(NUM_THDS) 
    { 
        gedge_t  edge;
        vertex_t v0, v1;
        
        index_t n, m;
        
        #pragma omp for
		for(index_t k = 0; k < nedges; ++k) {
			edge = edges[k];
			if (edge.is_self_loop()) continue;
	    
			v0 = edge.get_v0();
			v1 = edge.get_v1();
			
			n = _beg_pos[v0];
			m = __sync_fetch_and_add(count_edge + v0, 1);
			_adj[n + m] = v1;

			#ifdef HALF_GRID
			n = _beg_pos[v1];
			m = __sync_fetch_and_add(count_edge + v1, 1);
			_adj[n + m] = v0;
			#else
			n = _beg_pos_in[v1];
			m = __sync_fetch_and_add(count_edge_in + v1, 1);
			_adj_in[n + m] = v0;
			#endif	
		}
    }
   
    //munmap (edges, st_edge.st_size);
	//free(_adj);
	//free(_beg_pos);
	//#ifndef HALF_GRID
	//free(_adj_in);	
	//#endif

    close(fid_edge);
    fclose(f); 
    double end = mywtime();
    cout << "CSR conversion Time = " << end -  start << endl;
    cout << "classifcation done" << endl; 
}

void grid::compress_degree() 
{
	sdegree_t big_degree = 0;
	
	#pragma omp parallel for reduction(+:big_degree)
	for (vertex_t i = 0; i < vert_count; ++i) {
		big_degree += (vert_degree[i] > 32767);
	}
    
	++big_degree;
	svert_degree = (sdegree_t*)calloc(sizeof(sdegree_t), vert_count);
	bvert_degree = (bdegree_t*)calloc(sizeof(bdegree_t), big_degree);

	bdegree_count = big_degree;
	cout << bdegree_count << endl;
	big_degree = 1;
	bvert_degree[0] = 6000;//put any number
	
    #pragma omp parallel for 
	for (vertex_t i = 0; i < vert_count; ++i) {
		if ((vert_degree[i]) > 32767) {
			sdegree_t m = __sync_fetch_and_add(&big_degree, 1);
			bvert_degree[m] = vert_degree[i];
			svert_degree[i] = -m;
		} else {
			svert_degree[i] = vert_degree[i];
		}
	}
}

int 
main(int argc, char *argv[]) 
{
    g = new gConv;
    g->init(argc, argv);
    return 0;
}

void grid::save_csr(string edgefile)
{
    cout << "save CSR start" << endl;
    string file = edgefile + ".adj";
    FILE* f = fopen(file.c_str(), "wb");
    assert(f != 0);
    fwrite(_adj, sizeof(vertex_t), _beg_pos[vert_count], f);
    fclose(f);

    
    file = edgefile + ".beg_pos";
    f = fopen(file.c_str(), "wb");
    assert(f != 0);
    fwrite(_beg_pos, sizeof(index_t), vert_count + 1, f);
    fclose(f);
    cout << _beg_pos[vert_count] << endl;

    #ifndef HAL_GRID
    file = edgefile + ".adj_in";
    f = fopen(file.c_str(), "wb");
    assert(f != 0);
    fwrite(_adj_in, sizeof(vertex_t), _beg_pos_in[vert_count], f);
    fclose(f);
    
    file = edgefile + ".beg_pos_in";
    f = fopen(file.c_str(), "wb");
    assert(f != 0);
    fwrite(_beg_pos_in, sizeof(index_t), vert_count + 1, f);
    fclose(f);
    cout << _beg_pos_in[vert_count] << endl;
    #endif

    save_degree_files(edgefile);
} 

void grid::save_degree_files(string edgefile)
{
    file = edgefile + ".degree";
    f = fopen(file.c_str(), "wb");
    assert(f != 0);
    fwrite(vert_degree, sizeof(degree_t), vert_count, f);
    fclose(f);
	
    cout << "Compressing degree start" << endl;
	compress_degree();
	cout << "Compressing done" << endl;
    
	file = edgefile + ".sdegree";
    f = fopen(file.c_str(), "wb");
    assert(f != 0);
    fwrite(svert_degree, sizeof(sdegree_t), vert_count, f);
    fclose(f);
	
	file = edgefile + ".bdegree";
    f = fopen(file.c_str(), "wb");
    assert(f != 0);
    fwrite(bvert_degree, sizeof(bdegree_t), bdegree_count, f);
    cout << "saving done" << endl;
    fclose(f);
}
