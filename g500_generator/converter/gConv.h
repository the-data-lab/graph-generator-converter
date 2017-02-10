
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


#ifndef __GRID_H__
#define __GRID_H__
#include <stdint.h>
#include <string.h>

using namespace std;

#define handle_error(msg) \
                   do { perror(msg); exit(EXIT_FAILURE); } while (0)


typedef uint64_t index_t;

typedef int32_t degree_t;
typedef int32_t bdegree_t;
typedef int16_t sdegree_t;



#ifdef GENEARATOR_USE_48BITS
typedef uint64_t vertex_t;
struct gedge_t {
private:
    uint32_t v0_low;
    uint32_t v1_low;
    uint32_t high; /* v1 in high half, v0 in low half */

public:    
    inline bool is_self_loop() {
        return ((v0_low == v1_low) && 
		(high & 0xFFFF) == (high >> 16));
    }
    
    bool operator == (gedge_t e) {
        return ((v0_low == e.v0_low) && 
		(v1_low == e.v1_low) &&
		(high & 0xFFFF) == (high >> 16));
    }

    inline vertex_t get_v0()
    {
	    return (v0_low | ((int64_t)((int16_t)(high & 0xFFFF)) << 32));
    }
    inline vertex_t get_v1()
    {
	    return (v1_low | ((int64_t)((int16_t)(high >> 16)) << 32));
    }
    inline void set_v0(vertex_t a_v0) { assert(0); }
    inline void set_v1(vertex_t a_v1) { assert(0); } 
};
#elif  GENEARATOR_USE_64BITS
typedef uint64_t vertex_t;

struct gedge_t {
private:
    vertex_t v0;
    vertex_t_t v1;

public:
    inline bool is_self_loop() { return (v0 == v1); }
    
    bool operator == (gedge_t e) {
        return ((v0 == e.v0) && (v1 == e.v1));
    }

    inline vertex_t get_v0() { return v0; }
    inline vertex_t get_v1() { return v1; }
    inline void set_v0(vertex_t a_v0) { v0 = a_v0; }
    inline void set_v1(vertex_t a_v1) { v1 = a_v1; } 
};

#else
typedef uint32_t vertex_t;

struct gedge_t {
private:
    vertex_t v0;
    vertex_t_t v1;

public:
    inline bool is_self_loop() { return (v0 == v1); }
    
    bool operator == (gedge_t e) {
        return ((v0 == e.v0) && (v1 == e.v1));
    }

    inline vertex_t get_v0() { return v0; }
    inline vertex_t get_v1() { return v1; }
    inline void set_v0(vertex_t a_v0) { v0 = a_v0; }
    inline void set_v1(vertex_t a_v1) { v1 = a_v1; } 
};

#endif

class gConv {
public:
	
	void pre_csr(string edgefile, gedge_t* edges, index_t nedges);
    void proc_csr(string edgefile, string part_file);
	
	void compress_degree();
    void init(int argc, char* argv[]);
    void save_csr(string edgefile);
    void save_degree_files(string edgefile);

    
    

public:
    vertex_t    vert_count;
    index_t*    _beg_pos;
    vertex_t*   _adj; 
    
    #ifndef HALF_GRID
    index_t*    _beg_pos_in;
    vertex_t*   _adj_in; 
    #endif
    
	vertex_t    bdegree_count;
	degree_t*   vert_degree;
	sdegree_t*  svert_degree;
	bdegree_t*  bvert_degree;

};

#endif
