#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <list>
#include "make_graph.h"
#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <cstring>
#include <algorithm>

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>
#include <stdio.h>
#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "make_graph.h"

using std::endl;

#ifndef TIME_H
#define TIME_H


inline double get_time(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1.e-6;

}
#endif

#define DUMP_ADJ_BIN
#define DUMP_BEG_POS_BIN
#define DUMP_CSR
//#define DUMP_RMAT

typedef unsigned long	vertex_t;
typedef unsigned long	index_t;
typedef unsigned long data_t;
typedef long g500_t;

off_t fsize(const char *filename) {
    struct stat st;

    if (stat(filename, &st) == 0)
        return st.st_size;

//    fprintf(stderr, "Cannot determine size of %s: %s\n",
//            filename, strerror(errno));

    return -1;
}

//using 3 64bits to store 4 48bits verts.
//the first 3 is stored lower bits
//the 4th is stored in the higher bits across 3 64bits.
#ifdef COMPRESS
void adj_store(index_t off, vertex_t vert, vertex_t *adj_list){
  index_t base_off=(off>>2)*3;
  index_t group_off=off&0x03;
  if(group_off<3)
	{
		/*clear your bits*/
		adj_list[base_off+group_off]&=0xffff000000000000;
		
		/*set your bits*/
		adj_list[base_off+group_off]|=vert;
	}
  else
	{//group_off==3
		/*clear your bits*/
		adj_list[base_off]&=0x0000ffffffffffff;
		adj_list[base_off+1]&=0x0000ffffffffffff;
		adj_list[base_off+2]&=0x0000ffffffffffff;
    
		/*set your bits*/
		adj_list[base_off]|=((vert&0xffff)<<48);
    adj_list[base_off+1]|=((vert&0xffff0000)<<32);
    adj_list[base_off+2]|=((vert&0xffff00000000)<<16);
  }
}

vertex_t adj_load(index_t off, const vertex_t *adj_list){
  vertex_t res;
  index_t base_off=(off>>2)*3;
  index_t group_off=off&0x03;
  if(group_off<3) res=adj_list[base_off+group_off]&0xffffffffffff;
  else{
    res=(adj_list[base_off]>>48)+
        ((adj_list[base_off+1]&0xffff000000000000)>>32)+
        ((adj_list[base_off+2]&0xffff000000000000)>>16);
  }
  return res; 
}
#else
void adj_store(index_t off, vertex_t vert, vertex_t *adj_list){
	adj_list[off]=vert;
}

vertex_t adj_load(index_t off, const vertex_t *adj_list){
  return adj_list[off]; 
}
#endif

template<
typename packed_edge,
typename index_t,
typename data_t>
int 
transformer
(
	packed_edge *tuples,
	index_t tuple_count,
	data_t *adj_list,
	index_t *beg_pos,
	index_t *offset
){
	double report_criteria=0.0;	
	index_t i;

	double tm_taken=get_time();
//	#pragma omp parallel for \
//	private(i) num_threads(48)
	for(i=0; i<tuple_count; i++)
	{
		int tid=omp_get_thread_num();
		if(tid==0 && (i*1.0)/tuple_count > report_criteria){
			std::cout<<(100.0*i)/tuple_count<<"%\n";
			report_criteria += 0.1;
		}
		
		/*store one endpoint*/
		index_t ptr_id = (index_t)(tuples[i].v0);
		adj_store(beg_pos[ptr_id]+offset[ptr_id],
							(data_t)(tuples[i].v1), adj_list);
//		if((data_t)(tuples[i].v1)>1073741823)
//		{
//			std::cout<<"Wrong \n";
//			exit(-1);
//		}
//		#pragma omp atomic
		offset[ptr_id]++;
		ptr_id = (index_t)(tuples[i].v1);
		adj_store(beg_pos[ptr_id]+offset[ptr_id],
							(data_t)(tuples[i].v0), adj_list);
		
//		if((data_t)(tuples[i].v0)>1073741823)
//		{
//			std::cout<<"Wrong \n";
//			exit(-1);
//		}
		
//		#pragma omp atomic
		offset[ptr_id]++;
	}

	tm_taken = get_time() - tm_taken;
	std::cout<<"Time: "<<tm_taken<<" s\n"
		<<"Rate: "<<(tuple_count<<1)/tm_taken * 1.e-6<<"Medges/s\n";
	
	return 0;
}


template< 	typename data_t, typename index_t,
			typename g500_t>
bool csr_array( int 		argc,
				char**		argv,
				data_t*		&adj_list,
				index_t*	&beg_pos,
				index_t*	&adj_card,
				index_t		&edge_count,
				index_t		&vert_count)
{
	g500_t log_numverts;
	double start, tm_taken;
	FILE *fid;
	char filenames[1024];
	g500_t nedges;
	packed_edge* g500_tuples;
	double report_criteria = 0.0;
	std::cout<<"./exe log_vtx (default 16) edge_factor"
			<<" (default 16)\n";	
	
	log_numverts 			= 16; /* In base 2 */
	g500_t edge_factor	= 16;
	g500_t dump_count = 1;
	
	if(argc != 3) { 
        std::cout<<"Wrong input\n";
        exit(-1);
    }

	if(argc >= 2) log_numverts 	= atoi(argv[1]);
	if(argc >= 3) edge_factor	= atoi(argv[2]);
	

	vert_count	= 1 << log_numverts;
	
    /* Start of graph generation timing */
	start = omp_get_wtime();
	make_graph(log_numverts, edge_factor << log_numverts, \
				1, 2, &nedges, &g500_tuples);
	tm_taken = omp_get_wtime() - start;
	/* End of graph generation timing */

	edge_count = nedges;

	fprintf(stderr,"%"PRIu64"edge%s generated in %fs (%f Medges/s)\n",\
			edge_count,(edge_count == 1 ? "" : "s"),tm_taken,\
			1. * edge_count / tm_taken * 1.e-6);
	
	adj_card	= new index_t[vert_count];
	beg_pos		= new index_t[vert_count+1];

	/*double edge count regarding graph as undirected*/
	beg_pos[vert_count] = edge_count << 1;
	index_t     *offset	= new index_t[vert_count];
	
	memset(adj_card, 0, sizeof(index_t)*vert_count);
	memset(offset,   0, sizeof(index_t)*vert_count);
	
	start = get_time();
	std::cout<<"================Count==================\n";
	for(index_t i=0; i< edge_count; i++) {
		if((i*1.0)/edge_count > report_criteria){
			std::cout<<(100.0*i)/edge_count<<"%\n";
			report_criteria += 0.1;
		}
		
		/*assume it is an undirected edge*/
		adj_card[(index_t)(g500_tuples[i].v0)]++;
		adj_card[(index_t)(g500_tuples[i].v1)]++;
	}
	
	/*prefix sum based on out-degree*/
	std::cout<<"Prefix-sum ... \n";
	beg_pos[0] = 0;
	for(index_t i=0; i< vert_count; i++)
		beg_pos[i+1] = beg_pos[i] + adj_card[i];
	
	tm_taken = get_time() - start;
	std::cout << "Time: "<<tm_taken<<" s" << endl
		      << "Rate: "<< (edge_count << 1)/tm_taken * 1.e-6<<"Medges/s\n";
	
	#ifdef DUMP_BEG_POS_BIN
	/*dump the beg_pos into file*/
	sprintf(filenames,"beg_%lu_%lu.bin",log_numverts,edge_factor);
	fid = fopen(filenames,"wb");
	fwrite(beg_pos, sizeof(index_t), vert_count + 1, fid);
	fclose(fid);
	#endif
	
    
    #if 0
	/*dump adj_list: partitioned by vert count*/
    index_t bdump_sz = 0;
	index_t bdump_ptr = 0;
	index_t bstep_sz = vert_count / dump_count;
	char bfilenamec[1024];
	
    for(index_t i=0; i < dump_count; i++) {
		sprintf(bfilenamec,"beg_%lu_%lu.%d.bin", log_numverts, edge_factor, i);
		fid = fopen(bfilenamec,"wb");
		if (i == dump_count - 1) {
            bdump_sz = vert_count - bdump_ptr;
        } else {
            bdump_sz = bstep_sz;
        }
		
		std::cout <<"Beg-end-vert vs beg-pos-end-pos: " 
                  << bdump_ptr
                  <<" "
                  <<bdump_ptr+bdump_sz-1
                  <<"\n"
		          <<beg_pos[bdump_ptr]
                  <<" "
                  <<beg_pos[bdump_ptr+bdump_sz]
                  <<"\n";
		
        int ret = fwrite (beg_pos + bdump_ptr, sizeof(index_t), bdump_sz + 1, fid);
		fclose(fid);
		
        if(ret != bdump_sz + 1) {
            std::cout<<"Wrong dump1\n";exit(-1);
        }
		bdump_ptr += bdump_sz;
	}
    #endif

	/*dump half the tuple-list into binary file*/
	index_t half_tuple = nedges/2;
	if(half_tuple*2 != nedges) {
        std::cout<<"Wrong\n"; exit(-1);
    }

	fid = fopen("half-tuple-binary.dat","wb");
	fwrite(g500_tuples, sizeof(packed_edge), half_tuple,fid);
	fclose(fid);
	
	/*copy the other half to smaller array*/
	packed_edge *in_edges = new packed_edge[nedges - half_tuple];
	memcpy(in_edges, g500_tuples + half_tuple,
				sizeof(packed_edge)*(nedges-half_tuple));

	
    #ifdef DUMP_TUPLE
	/*dump tuple list*/
	sprintf(filenames,"tuple_%lu_%lu.dat",log_numverts,edge_factor);
	fid = fopen(filenames,"w");
	for(index_t i = 0; i < nedges; i++) {
		fprintf(fid, "%lu,%lu\n", (vertex_t)(g500_tuples[i].v0),
						(vertex_t)(g500_tuples[i].v1));
	}
	fclose(fid);
	#endif

	free(g500_tuples);
	
    /*alloc space*/
	#ifdef COMPRESS
	index_t edge_comp_count = ((edge_count*2)>>2)*3+((edge_count*2)&0x03);
	#else
	index_t edge_comp_count= edge_count*2;
	#endif
	
    adj_list	= new data_t[edge_comp_count];
	memset(adj_list, 0, sizeof(data_t)*edge_comp_count);	
	
	/*store the first half into adj_list -- 48 bits*/
	transformer<packed_edge,index_t,data_t> (in_edges, half_tuple, adj_list, beg_pos, offset);

	/*load the dumped half into memory*/
	fid = fopen("half-tuple-binary.dat","rb");
	
    //off_t size_offset = fsize("half-tuple-binary.dat");
	index_t ret = fread(in_edges, sizeof(packed_edge), half_tuple, fid);
	
    if(ret != half_tuple) {
        std::cout<<"Wrong\n";exit(-3);
    }

	fclose(fid);
	
	/*remove half-tuple-binary.dat file*/
	if(remove("half-tuple-binary.dat")!=0)
	{
		std::cout<<"Error deleting half-tuple-binary.dat\n";
		exit(-1);
	}

	/*store the second half into adj_list -- 48 bits*/
	transformer<packed_edge,index_t,data_t>
	(in_edges,half_tuple,adj_list,beg_pos,offset);

    #pragma omp parallel for  num_threads(54) 
	for(index_t i = 0; i < vert_count; i++) {
		index_t a = beg_pos[i];
		index_t b = beg_pos[i + 1];
		//quickSort(adj_list, a, b);
        std::sort(adj_list+a, adj_list+b);
	}
  
	#ifdef DUMP_ADJ_BIN
	/*dump adj_list: binary -> partitioned by edge count*/
/*	index_t dump_sz=edge_comp_count/dump_count;
	index_t dump_ptr=0;
	char filenamec[1024];
	for(index_t i=0;i<dump_count;i++)
	{
		sprintf(filenamec,"csr_%lu_%lu.%d.bin",log_numverts,edge_factor,i);
		fid=fopen(filenamec,"wb");
		int ret=fwrite(adj_list+dump_ptr,sizeof(data_t),dump_sz,fid);
		fclose(fid);
		if(ret!=dump_sz){std::cout<<"Wrong dump\n";exit(-1);}
		dump_ptr+=dump_sz;

		/ *next round is the last dump* /
		if(i==dump_count-2) dump_sz=edge_comp_count-dump_sz*(dump_count-1);
	}
	*/	

	/*dump adj_list: partitioned by vert count*/
	index_t dump_sz = 0;
	index_t dump_ptr = 0;
	index_t step_sz = vert_count / dump_count;
	char filenamec[1024];
	
    for(index_t i=0;i<dump_count;i++) {
		sprintf(filenamec,"csr_%lu_%lu.%d.bin",log_numverts,edge_factor,i);
		fid = fopen(filenamec,"wb");
		
        if( i == dump_count - 1) dump_sz = edge_comp_count - dump_ptr;
		else dump_sz = beg_pos[(i+1)*step_sz] - dump_ptr;
		
		std::cout<<"beg_pos vs end-pos: "<<dump_ptr<<" "<<dump_ptr+dump_sz<<"\n";

		int ret = fwrite(adj_list+dump_ptr,sizeof(data_t),dump_sz,fid);
		fclose(fid);
		if( ret != dump_sz) {
            std::cout<<"Wrong dump2\n";
            exit(-1);
        }
		dump_ptr+=dump_sz;
	}
	#endif

    //exit(-1);
	#ifdef DUMP_CSR
	/*dump adj_list: text format*/
    sprintf(filenames,"csr_%lu_%lu.dat",log_numverts,edge_factor);
	std::ofstream outfile(filenames);

    for(index_t i=0; i<vert_count; i++) {
        outfile << i << " " << adj_card[i] << " ";
        for(index_t j = beg_pos[i]; j < beg_pos[i+1]; j++) {
            outfile<<adj_load(j,adj_list)<<" ";
		}
        outfile<<"\n";
    }
    outfile.close();
	#endif
	
    #ifdef DUMP_RMAT
	/*dump adj_list: text format*/
    sprintf(filenames,"rmat_%lu_%lu.dat",log_numverts,edge_factor);
	std::ofstream outfile1(filenames);
    outfile1 << "AdjacencyGraph" << endl;
    outfile1 << vert_count << endl;
    outfile1 << beg_pos[vert_count] << endl;
    
    for(index_t i=0; i<vert_count; i++) {
        outfile1<< beg_pos[i] <<endl;
    }
  
    for(index_t i=0; i<vert_count; i++) {
        for(index_t j = beg_pos[i]; j < beg_pos[i+1]; j++) {
            outfile1 << adj_load(j,adj_list) << endl;
		}
    }
    outfile1.close();
    #endif

	/*output edge.csv*/
	#ifdef DUMP_EDGE
	sprintf(filenames,"edge.csv");
	fid = fopen(filenames,"w");
	
    for(index_t i = 0; i < vert_count; i++) {
        for(index_t j=beg_pos[i]; j<beg_pos[i+1]; j++) {
			fprintf(fid,"%lu,%lu\n",i,adj_load(j,adj_list));
		}
	}
	fclose(fid);
	#endif
	
    return true;
}
