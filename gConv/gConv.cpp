
/*
 * Copyright 2016 The George Washington University
 * Written by Pradeep Kumar 
 * Directed by Prof. Howie Huang
 *
 * https://www.seas.gwu.edu/~howie/
 * Contact: iheartgraph@gmail.com
 *
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
#include <dirent.h>
#include <asm/mman.h>
#include "wtime.h"
#include "gConv.h"

//XXX
#define NUM_THDS 56



gConv* g;

inline off_t fsize(const char *filename) {
    struct stat st; 
    if (stat(filename, &st) == 0)
        return st.st_size;
    return -1; 
}

void text_to_bin(string textfile, string ofile)
{
	int fd;
	char* ss_head;
	char* ss;

	size_t file_size = fsize(textfile.c_str());
	fd = open( textfile.c_str(), O_CREAT|O_RDWR, 00777);

	ss_head = (char*)mmap(NULL,file_size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);

	size_t head_offset=0;
    
    //remove initial comments, 
    //generally present in real-world graphs downloaded from internet.
	while(ss_head[head_offset]=='%'){
		while(ss_head[head_offset]!='\n'){
			head_offset++;
		}
		head_offset++;
	}
	ss = ss_head + head_offset;
	file_size -= head_offset;

	size_t curr=0;
	size_t next=0;

	//step 4: write adjacent list 
	vertex_t v0;
    vertex_t v1;
	size_t offset =0;
	next = 0;
	curr = 0;

	//int fd4 = open( ofile.c_str(), O_CREAT|O_RDWR,00777 );
	//ftruncate(fd4, edge_count*sizeof(gedge_t));
	//gedge_t* adj = (gedge_t*)mmap(NULL,edge_count*sizeof(gedge_t),
    //                                PROT_READ|PROT_WRITE,MAP_SHARED,fd4,0);
	//gedge_t* adj = (gedge_t*)malloc(edge_count*sizeof(gedge_t));
	gedge_t* adj = (gedge_t*)malloc(2*file_size);//allocating more memory, but should be ok
	
	while(next < file_size) {
	    char* sss = ss+curr;
	    v0 = atoi(sss);//first end of pair
	    while((ss[next]!=' ')&&(ss[next]!='\n')&&(ss[next]!='\t')){
		    next++;
	    }
	    while((ss[next]==' ')||(ss[next]=='\n')||(ss[next]=='\t')){
		    next++;
	    }
	    curr = next;

	    char* sss1=ss+curr;
	    v1 = atoi(sss1);//second end of pair
         
	    adj[offset].set_v0(v0);
	    adj[offset].set_v1(v1);
		

	    while((ss[next]!=' ')&&(ss[next]!='\n')&&(ss[next]!='\t')){
		    next++;
	    }
	    while((ss[next]==' ')||(ss[next]=='\n')||(ss[next]=='\t')){
		    next++;
	    }
	    curr = next;

	    offset++;
	}
	
	munmap(ss_head,sizeof(char)*file_size);
	//munmap( adj,sizeof(vertex_t)*edge_count );
	close(fd);
	FILE* fd4 = fopen(ofile.c_str(), "wb");
    assert(fd4 != 0);
    fwrite(adj, sizeof(gedge_t), offset, fd4);
    
    fclose(fd4);
}

void text_to_bin_manyfiles(string idir, string odir)
{
    struct dirent *ptr;
    DIR *dir;
    int file_count = 0;
    string ifile[1024];
    string file;
    string ofile;
    
    //Read graph file
    dir = opendir(idir.c_str());
    while (NULL != (ptr = readdir(dir))) {
        if (ptr->d_name[0] == '.') continue;
        
        ifile[file_count] =  string(ptr->d_name);
        file_count++;
        
    }
    closedir(dir);
    #pragma omp parallel for schedule(dynamic, 2)
    for (int i = 0; i < file_count; i++) {
        file = idir + "/" + ifile[i];
        ofile = odir + "/" + ifile[i] + ".edge"; 
        text_to_bin(file, ofile);
    }
    
}

void text_to_bin_id_translate(string textfile)
{
	int fd;
	char* ss_head;
	char* ss;

	size_t file_size = fsize(textfile.c_str());
	fd = open( textfile.c_str(), O_CREAT|O_RDWR, 00777);

	ss_head = (char*)mmap(NULL,file_size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);

	size_t head_offset=0;
	while(ss_head[head_offset]=='%'){
		while(ss_head[head_offset]!='\n'){
			head_offset++;
		}
		head_offset++;
	}
	ss = &ss_head[head_offset];
	file_size -= head_offset;

	size_t curr=0;
	size_t next=0;

	//step 1. vert_count,edge_count,
	size_t edge_count=0;
	size_t vert_count;
	uint32_t v_max = 0;
	uint32_t v_min = 999999;//as infinity
	vertex_t a;
	while(next<file_size){
		char* sss=ss+curr;
		a = atoi(sss);

		if(v_max<a){
			v_max = a;
		}
		if(v_min>a){
			v_min = a;
		}

		while((ss[next]!=' ')&&(ss[next]!='\n')&&(ss[next]!='\t')){
			next++;
		}
		while((ss[next]==' ')||(ss[next]=='\n')||(ss[next]=='\t')){
			next++;
		}
		curr = next;
        edge_count++; 
	}
	edge_count /=2;
	vert_count = v_max - v_min + 1;

	cerr<<"max vertex id: "<<v_max<<endl;
	cerr<<"min vertex id: "<<v_min<<endl;

	cerr<<"edge count: "<<edge_count<<endl;
	cerr<<"vert count: "<<vert_count<<endl;

	//step 4: write adjacent list 
	uint32_t v0;
    uint32_t v1;
	size_t offset =0;
	next = 0;
	curr = 0;

    string edgefile = textfile + ".edge";    
	//int fd4 = open( edgefile.c_str(), O_CREAT|O_RDWR,00777 );
	//ftruncate(fd4, edge_count*sizeof(gedge_t));
	//gedge_t* adj = (gedge_t*)mmap(NULL,edge_count*sizeof(gedge_t),
    //                                PROT_READ|PROT_WRITE,MAP_SHARED,fd4,0);
	gedge_t* adj = (gedge_t*)malloc(edge_count*sizeof(gedge_t));
	
	while(next < file_size) {
	    char* sss = ss+curr;
	    v0 = atoi(sss)-v_min;//first end of pair
	    while((ss[next]!=' ')&&(ss[next]!='\n')&&(ss[next]!='\t')){
		    next++;
	    }
	    while((ss[next]==' ')||(ss[next]=='\n')||(ss[next]=='\t')){
		    next++;
	    }
	    curr = next;

	    char* sss1=ss+curr;
	    v1 = atoi(sss1)-v_min;//second end of pair
         
	    adj[offset].set_v0(v0);
	    adj[offset].set_v1(v1);
		

	    while((ss[next]!=' ')&&(ss[next]!='\n')&&(ss[next]!='\t')){
		    next++;
	    }
	    while((ss[next]==' ')||(ss[next]=='\n')||(ss[next]=='\t')){
		    next++;
	    }
	    curr = next;

	    offset++;
	}
	
	munmap( ss,sizeof(char)*file_size );
	//munmap( adj,sizeof(vertex_t)*edge_count );
	close(fd);
	FILE* fd4 = fopen(edgefile.c_str(), "wb");
    assert(fd4 != 0);
    fwrite(adj, sizeof(gedge_t), edge_count, fd4);
    
    fclose(fd4);
}

void gConv::init(int argc, char * argv[])
{
    int o;
    uint32_t scale;
    int c = 1000;
    string edgefile;
    string part_file;
    while ((o = getopt (argc, argv, "s:o:hi:j:c:m:v:a:")) != -1) {
        switch(o) {
            case 's': //scale
                scale = atoi(optarg);
                vert_count = (1L << scale);
                break;
            case 'v'://vert count
                //vert_count = atoi(optarg);
				sscanf(optarg, "%ld", &vert_count);
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

    double start, end;
  
    switch(c) {
    case 0:
		start = mywtime();
		text_to_bin(edgefile, part_file);//input is text file
		end = mywtime();
		cout << "Time = " << end - start << endl;
		return;
    case 1:
		start = mywtime();
		//Both arguments are directory
        text_to_bin_manyfiles(edgefile, part_file);//input/output are dir
		end = mywtime();
		cout << "Time = " << end - start << endl;
		return;
	case 2:
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

void gConv::pre_csr(string edgefile, gedge_t* edges, index_t nedges)
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
void gConv::proc_csr(string edgefile, string part_file)
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

void gConv::compress_degree() 
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

void gConv::save_csr(string edgefile)
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

    #ifndef HALF_GRID
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

void gConv::save_degree_files(string edgefile)
{
    string file = edgefile + ".degree";
    FILE* f = fopen(file.c_str(), "wb");
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
