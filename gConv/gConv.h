
/*
 * Copyright 2016 The George Washington University
 * Written by Pradeep Kumar 
 * Directed by Prof. Howie Huang
 *
 * https://www.seas.gwu.edu/~howie/
 * Contact: iheartgraph@gmail.com
 *
 */ 


#ifndef __GCONV_H__
#define __GCONV_H__
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
    vertex_t v1;

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
	
    void proc_csr(string edgefile, string part_file);
    void proc_csr_dir(string edgefile, string part_file);
    void proc_csr_rank(string edgefile, string part_file, int rank_by);
    void init(int argc, char* argv[]);
	
private:	
	void pre_csr(gedge_t* edges, index_t nedges);
	void pre_csr_rank(string edgefile, gedge_t* edges, index_t nedges, int rank_by);
    void proc_csr_rankbydegree(string csrfile, string part_file, int rank_by);
	void compress_degree();
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
