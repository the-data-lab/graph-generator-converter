
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
//#include <libaio.h>
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

void print_csr(string part_file)
{
    string file_etable = part_file + ".adj_rankbydegree";
    string file_vtable = part_file + ".beg_pos_rankbydegree";
    FILE* f_vtable = fopen(file_vtable.c_str(), "rb");
    FILE* f_etable = fopen(file_etable.c_str(), "rb");
    
    //Number of vertex
    int fid_vtable = open(file_vtable.c_str(), O_RDONLY);
    struct stat st_vert;
    fstat(fid_vtable, &st_vert);
    assert(st_vert.st_size != 0);
    vertex_t v_count = st_vert.st_size/sizeof(index_t) - 1;
    close(fid_vtable);
    cout << "v_count:" << v_count << endl;
    
    index_t* vtable = (index_t*)calloc(sizeof(index_t), v_count + 1);
    fread(vtable, sizeof(index_t), v_count+1, f_vtable);
    fclose(f_vtable);


    vertex_t* etable= (vertex_t*)calloc(sizeof(vertex_t), vtable[v_count]);
    fread(etable, sizeof(index_t), vtable[v_count], f_etable);
    fclose(f_etable);

    for (vertex_t v = 0; v < v_count; ++v) {
        cout << "V = " << v << endl;
        vertex_t degree = vtable[v+1] - vtable[v];
        vertex_t* adj = etable + vtable[v];
        for (vertex_t d = 0; d < degree; d++) {
            cout << adj[d] << " "; 
        }
        cout << endl;
    }

}

void print_edgetuple(string edgefile)
{

}

void text_to_bin(string textfile, string ofile)
{
	int fd;
	char* ss_head;
	char* ss;
    vertex_t min_v = 100000;//infinity
    vertex_t max_v = 0;

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
	gedge_t* adj = (gedge_t*)malloc(file_size*4);//allocating more memory, but should be ok
	
	while(next < file_size) {
	    char* sss = ss+curr;
	    v0 = atoi(sss);//first end of pair
        min_v = min(min_v, v0);
        max_v = max(max_v, v0);
	    while((ss[next]!=' ')&&(ss[next]!='\n')&&(ss[next]!='\t')){
		    next++;
	    }
	    while((ss[next]==' ')||(ss[next]=='\n')||(ss[next]=='\t')){
		    next++;
	    }
	    curr = next;

	    char* sss1=ss+curr;
	    v1 = atoi(sss1);//second end of pair
        min_v = min(min_v, v1);
        max_v = max(max_v, v1);
         
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
    free(adj); 
    fclose(fd4);
    cout << "min_v = " << min_v << endl;
    cout << "max_v = " << max_v << endl;
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
    cout << file_count << endl;
    closedir(dir);
    #pragma omp parallel  
    {
        string file, ofile; 
        #pragma omp for 
        for (int i = 0; i < file_count; i++) {
            file = idir + "/" + ifile[i];
            ofile = odir + "/" + ifile[i] + ".edge";
            cout << "started " << file << endl;
            text_to_bin(file, ofile);
            cout << "done " << file << endl;
        }
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

    string edgefile;
    string part_file;
    int o;
    uint32_t scale;
    int c = 1000;
	//rank by: 0: no ranking, 1: rank by degree
	int rank = 0;
    
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
			case 'r':
				rank = atoi(optarg);
            default:
               cout << "Unknown argument" << endl ; 
        }
    }

    double start, end;
    cout << vert_count << endl; 
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
        //directory, output file name
		proc_csr_dir(edgefile, part_file);
        save_csr(part_file);
		end = mywtime();
		cout << "CSR Time = " << end - start << endl;
		return;
        return;
	case 3:
		start = mywtime();
		proc_csr(edgefile, part_file);
        save_csr(part_file);
		end = mywtime();
		cout << "CSR Time = " << end - start << endl;
		return;
    case 4: //clean of duplication, self loop 
		start = mywtime();
        //csr file as input, csr file as output 
		proc_csr_clean(edgefile, part_file);
        //save_csr(part_file);
		end = mywtime();
		cout << "Clean CSR Time = " << end - start << endl;
		return;
    case 5:
		start = mywtime();
		proc_csr_rankbig(edgefile, part_file, rank);
		end = mywtime();
		cout << "Rank by degree Time = " << end - start << endl;
		return;
	case 6:
		start = mywtime();
		proc_csr_rank(edgefile, part_file, rank);
        save_csr(part_file);
		end = mywtime();
		cout << "CSR Time = " << end - start << endl;
		return;

    case 7:
        print_csr(edgefile);
        return;
    default:
        return;
    } 

   
}

void gConv::pre_csr_rank(string edgefile, gedge_t* edges, index_t nedges, int rank_by)
{
    _beg_pos = (index_t*) calloc(sizeof(index_t), vert_count + 1);  
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
            
			if (vert_degree[v1] >= vert_degree[v0]) {
				__sync_fetch_and_add(vert_degree + v0, 1);
			
			} else {
				__sync_fetch_and_add(vert_degree + v1, 1);
			}
        }
    }
		
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
		
			if (vert_degree[v1] >= vert_degree[v0]) {
				__sync_fetch_and_add(_beg_pos + v0, 1);
			} else {
				__sync_fetch_and_add(_beg_pos + v1, 1);
			}
		}
	}
    
	//Calculate the CSR beg_pos
    index_t prefix_sum = 0;
	index_t curr_value = 0;
    
	for (index_t ipart = 0; ipart < vert_count; ++ipart) {
		curr_value = _beg_pos[ipart];
        _beg_pos[ipart] = prefix_sum;
        prefix_sum += curr_value;
    }

    _beg_pos[vert_count] = prefix_sum;
    cout << "Total edges = " << prefix_sum << endl;
}

//CSR undirected to CSR rank by degree input
//First clean and remove duplication, self loop
void gConv::proc_csr_clean(string csrfile, string part_file) 
{
    string file_etable = csrfile + ".adj";
    string file_vtable = csrfile + ".beg_pos"; 

    //New files
    string file_newetable = part_file + ".adj_clean";
    string file_newvtable = part_file + ".beg_pos_clean";
    string file_newdegree = part_file + ".degree_clean";
    FILE* f_newvtable = fopen(file_newvtable.c_str(), "wb");
    FILE* f_newetable = fopen(file_newetable.c_str(), "wb");
    FILE* f_newdegree = fopen(file_newdegree.c_str(), "wb");

    //Number of vertex
    int fid_vtable = open(file_vtable.c_str(), O_RDONLY);
    struct stat st_vert;
    fstat(fid_vtable, &st_vert);
    assert(st_vert.st_size != 0);
    vertex_t v_count = st_vert.st_size/sizeof(index_t) - 1;
    close(fid_vtable);
    cout << "v_count:" << v_count << endl;

    //memory allocation
    index_t  max_ecount = (1L << 28);
    index_t  max_vcount = (1L << 26);
    index_t* vtable     = (index_t*)calloc(sizeof(index_t), v_count + 1);
    index_t* new_vtable = (index_t*)calloc(sizeof(index_t), v_count + 1);
    vertex_t* etable    = (vertex_t*)calloc(sizeof(vertex_t), max_ecount);
    vertex_t* new_etable= (vertex_t*)calloc(sizeof(vertex_t), max_ecount);
    vertex_t* new_degree= (vertex_t*)calloc(sizeof(vertex_t), v_count);
    
    //read beg pos file
    FILE* f_etable = fopen(file_etable.c_str(), "rb");
    FILE* f_vtable = fopen(file_vtable.c_str(), "rb");
    fread(vtable, sizeof(index_t), v_count+1, f_vtable);
    fclose(f_vtable);
    
    index_t     e_count = 0;
    vertex_t    init_v  = 0;
    vertex_t    degree  = 0;
    index_t     prefix  = 0;
    index_t     new_ecount = 0;
    vertex_t*   adj_list;
    vertex_t*   new_adjlist;
    

    for (vertex_t v = 0; v < v_count; ++v) {
        degree = vtable[v+1] - vtable[v];
        if ((e_count + degree < max_ecount) && (v - init_v < max_vcount)) {
            e_count += degree;
            continue;
        }

        //read edge file
        fread(etable, sizeof(vertex_t), e_count, f_etable);
        
        index_t e_prefix = vtable[init_v];
        
        cout << "V:" << v << ":" << e_count << endl;
    

        //sort, remove duplication and self loop
        #pragma omp parallel for private(degree, adj_list) reduction(+:new_ecount) 
        for (vertex_t u = init_v; u < v; ++u) {
            adj_list = etable + vtable[u] - e_prefix;
            degree = vtable[u+1] - vtable[u];
            sort(adj_list, adj_list + degree);
            
            if (degree != 0)  {
                if (u != adj_list[0]) {
                    new_vtable[u]++;
                    new_ecount++;
                }

                for (vertex_t d = 1; d < degree; ++d) {
                    if ((adj_list[d-1] != adj_list[d]) && (u != adj_list[d])) {
                        new_vtable[u]++;
                        new_ecount++;
                    }
                }
            }
        }
        //prepare new beg pos
        for (vertex_t u = init_v; u < v; ++u) {
            degree = new_vtable[u];
            new_vtable[u] = prefix;
            prefix += degree; 
        }
        memset(new_degree, 0, sizeof(vertex_t)*(v - init_v));
        index_t new_eprefix = new_vtable[init_v];

        //prepare new adj list
        #pragma omp parallel for private(degree, adj_list, new_adjlist) 
        for (vertex_t u = init_v; u < v; ++u) {
            degree = vtable[u+1] - vtable[u];
            adj_list = etable + vtable[u] - e_prefix;
            new_adjlist = new_etable + new_vtable[u] - new_eprefix;

            if (degree != 0) {
                if (u != adj_list[0]) {
                    new_adjlist[new_degree[u]] = adj_list[0];
                    new_degree[u]++;
                }
            
                for (vertex_t d = 1; d < degree; ++d) {
                    if ((adj_list[d-1] != adj_list[d]) && (u != adj_list[d])) {
                        new_adjlist[new_degree[u]] = adj_list[d];
                        new_degree[u]++;
                    }
                }
            }
        }

        //write the new adj list
        fwrite(new_etable, sizeof(vertex_t), new_ecount, f_newetable);
        cout << "Wrote :" << v << ":" << new_ecount << endl;

        //update the variable for bookkeeping.
        init_v = v;
        e_count = 0;
        new_ecount = 0;

        //Handle the current vertex again.
        degree = vtable[v+1] - vtable[v];
        if (e_count + degree < max_ecount) {
            e_count += degree;
        } else {
            cout << "increase the memory size" << endl;
            assert(0);
        }
    }
    //last part not covered
    vertex_t v = v_count;
    
    //read edge file
    fread(etable, sizeof(vertex_t), e_count, f_etable);
    
    index_t e_prefix = vtable[init_v];
    
    cout << "V:" << v << ":" << e_count << endl;


    //sort, remove duplication and self loop
    #pragma omp parallel for private(degree, adj_list) reduction(+:new_ecount) 
    for (vertex_t u = init_v; u < v; ++u) {
        adj_list = etable + vtable[u] - e_prefix;
        degree = vtable[u+1] - vtable[u];
        sort(adj_list, adj_list + degree);
        
        if (degree != 0)  {
            if (u != adj_list[0]) {
                new_vtable[u]++;
                new_ecount++;
            }

            for (vertex_t d = 1; d < degree; ++d) {
                if ((adj_list[d-1] != adj_list[d]) && (u != adj_list[d])) {
                    new_vtable[u]++;
                    new_ecount++;
                }
            }
        }
    }

    //prepare new beg pos
    for (vertex_t u = init_v; u < v; ++u) {
        degree = new_vtable[u];
        new_vtable[u] = prefix;
        prefix += degree; 
    }
    new_vtable[v] = prefix;//special line 
    index_t new_eprefix = new_vtable[init_v];

    //prepare new adj list
    #pragma omp parallel for private(degree, adj_list, new_adjlist) 
    for (vertex_t u = init_v; u < v; ++u) {
        degree = vtable[u+1] - vtable[u];
        adj_list = etable + vtable[u] - e_prefix;
        new_adjlist = new_etable + new_vtable[u] - new_eprefix;
        
        if(degree != 0) {
            if (u != adj_list[0]) {
                new_adjlist[new_degree[u]] = adj_list[0];
                new_degree[u]++;
            }
        
            for (vertex_t d = 1; d < degree; ++d) {
                if ((adj_list[d-1] != adj_list[d]) || (u == adj_list[d])) {
                    new_adjlist[new_degree[u]] = adj_list[d];
                    new_degree[u]++;
                }
            }
        }
    }

    //write the new adj list
    fwrite(new_etable, sizeof(vertex_t), new_ecount, f_newetable);
    cout << "Wrote :" << v << ":" << new_ecount << endl;

    //Everything done
    fclose(f_newetable); 
    //write the new vertex table    
    fwrite(new_vtable, sizeof(index_t), v_count + 1, f_newvtable);
    fclose(f_newvtable);
    fwrite(new_degree, sizeof(vertex_t), v_count, f_newdegree);
    fclose(f_newdegree);
}

//CSR cleaned to CSR rank by degree
void gConv::proc_csr_rankbig(string csrfile, string part_file, int rank) 
{
    string file_etable = csrfile + ".adj_clean";
    string file_vtable = csrfile + ".beg_pos_clean"; 

    //New files
    string file_newetable = part_file + ".adj_rankbydegree";
    string file_newvtable = part_file + ".beg_pos_rankbydegree";
    FILE* f_newvtable = fopen(file_newvtable.c_str(), "wb");
    FILE* f_newetable = fopen(file_newetable.c_str(), "wb");

    //Number of vertex
    int fid_vtable = open(file_vtable.c_str(), O_RDONLY);
    struct stat st_vert;
    fstat(fid_vtable, &st_vert);
    assert(st_vert.st_size != 0);
    vertex_t v_count = st_vert.st_size/sizeof(index_t) - 1;
    close(fid_vtable);
    cout << "v_count:" << v_count << endl;

    //memory allocation
    index_t  max_ecount = (1L << 28);
    index_t  max_vcount = (1L << 26);
    index_t* vtable     = (index_t*)calloc(sizeof(index_t), v_count + 1);
    index_t* new_vtable = (index_t*)calloc(sizeof(index_t), v_count + 1);
    vertex_t* etable    = (vertex_t*)calloc(sizeof(vertex_t), max_ecount);
    vertex_t* new_etable= (vertex_t*)calloc(sizeof(vertex_t), max_ecount);
    vertex_t* new_degree= (vertex_t*)calloc(sizeof(vertex_t), v_count);
    
    //read beg pos file
    FILE* f_etable = fopen(file_etable.c_str(), "rb");
    FILE* f_vtable = fopen(file_vtable.c_str(), "rb");
    fread(vtable, sizeof(index_t), v_count+1, f_vtable);
    fclose(f_vtable);
    
    index_t     e_count = 0;
    vertex_t    init_v  = 0;
    vertex_t    degree  = 0;
    index_t     prefix  = 0;
    index_t     new_ecount = 0;
    vertex_t*   adj_list;
    vertex_t*   new_adjlist;
    

    for (vertex_t v = 0; v < v_count; ++v) {
        degree = vtable[v+1] - vtable[v];
        if ((e_count + degree < max_ecount) && (v - init_v < max_vcount)) {
            e_count += degree;
            continue;
        }

        //read edge file
        fread(etable, sizeof(vertex_t), e_count, f_etable);
        
        index_t e_prefix = vtable[init_v];
        
        cout << "V:" << v << ":" << e_count << endl;
    

        //do rank by degree
        #pragma omp parallel for private(degree, adj_list) reduction(+:new_ecount) 
        for (vertex_t u = init_v; u < v; ++u) {
            adj_list = etable + vtable[u] - e_prefix;
            degree = vtable[u+1] - vtable[u];
            vertex_t degree2 = 0;
            
            for (vertex_t d = 0; d < degree; ++d) {
                vertex_t nebr = adj_list[d];
                degree2 = vtable[nebr+1] - vtable[nebr];
                if ((degree < degree2) ||(degree == degree2 && u < nebr)) {
                    new_vtable[u]++;
                    new_ecount++;
                }
            }
        }

        //prepare new beg pos
        for (vertex_t u = init_v; u < v; ++u) {
            degree = new_vtable[u];
            new_vtable[u] = prefix;
            prefix += degree; 
        }
        
        index_t new_eprefix = new_vtable[init_v];

        //prepare new adj list
        #pragma omp parallel for private(degree, adj_list, new_adjlist) 
        for (vertex_t u = init_v; u < v; ++u) {
            degree = vtable[u+1] - vtable[u];
            adj_list = etable + vtable[u] - e_prefix;
            new_adjlist = new_etable + new_vtable[u] - new_eprefix;
            vertex_t degree2 = 0;

            for (vertex_t d = 0; d < degree; ++d) {
                vertex_t nebr = adj_list[d];
                degree2 = vtable[nebr+1] - vtable[nebr];
                if ((degree < degree2) ||(degree == degree2 && u < nebr)) {
                    new_adjlist[new_degree[u]] = adj_list[d];
                    new_degree[u]++;
                }
            }
        }

        //write the new adj list
        fwrite(new_etable, sizeof(vertex_t), new_ecount, f_newetable);
        cout << "Wrote :" << v << ":" << new_ecount << endl;

        //update the variable for bookkeeping.
        init_v = v;
        e_count = 0;
        new_ecount = 0;

        //Handle the current vertex again.
        degree = vtable[v+1] - vtable[v];
        if (e_count + degree < max_ecount) {
            e_count += degree;
        } else {
            cout << "increase the memory size" << endl;
            assert(0);
        }
    }
    //last part not covered
    vertex_t v = v_count;
        
    //read edge file
    fread(etable, sizeof(vertex_t), e_count, f_etable);
    
    index_t e_prefix = vtable[init_v];
    
    cout << "V:" << v << ":" << e_count << endl;


    //do rank by degree
    #pragma omp parallel for private(degree, adj_list) reduction(+:new_ecount) 
    for (vertex_t u = init_v; u < v; ++u) {
        adj_list = etable + vtable[u] - e_prefix;
        degree = vtable[u+1] - vtable[u];
        vertex_t degree2 = 0;
        
        for (vertex_t d = 0; d < degree; ++d) {
            vertex_t nebr = adj_list[d];
            degree2 = vtable[nebr+1] - vtable[nebr];
            if ((degree < degree2) ||(degree == degree2 && u < nebr)) {
                new_vtable[u]++;
                new_ecount++;
            }
        }
    }

    //prepare new beg pos
    for (vertex_t u = init_v; u < v; ++u) {
        degree = new_vtable[u];
        new_vtable[u] = prefix;
        prefix += degree; 
    }
    new_vtable[v] = prefix;//special line 
    index_t new_eprefix = new_vtable[init_v];

    //prepare new adj list
    #pragma omp parallel for private(degree, adj_list, new_adjlist) 
    for (vertex_t u = init_v; u < v; ++u) {
        degree = vtable[u+1] - vtable[u];
        adj_list = etable + vtable[u] - e_prefix;
        new_adjlist = new_etable + new_vtable[u] - new_eprefix;
        vertex_t degree2 = 0;

        for (vertex_t d = 0; d < degree; ++d) {
            vertex_t nebr = adj_list[d];
            degree2 = vtable[nebr+1] - vtable[nebr];
            if ((degree < degree2) ||(degree == degree2 && u < nebr)) {
                new_adjlist[new_degree[u]] = adj_list[d];
                new_degree[u]++;
            }
        }
    }

    //write the new adj list
    fwrite(new_etable, sizeof(vertex_t), new_ecount, f_newetable);
    cout << "Wrote :" << v << ":" << new_ecount << endl;
    

    //Everything done
    fclose(f_newetable); 
    //write the new vertex table    
    fwrite(new_vtable, sizeof(index_t), v_count + 1, f_newvtable);
    fclose(f_newvtable);
}

void gConv::proc_csr_rank(string edgefile, string part_file, int rank_by)
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

    pre_csr_rank(edgefile, edges, nedges, rank_by);

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
		
			if (vert_degree[v1] >= vert_degree[v0]) {
				n = _beg_pos[v0];
				m = __sync_fetch_and_add(count_edge + v0, 1);
				_adj[n + m] = v1;
			} else {	
				n = _beg_pos[v1];
				m = __sync_fetch_and_add(count_edge + v1, 1);
				_adj[n + m] = v0;
			}
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

void gConv::pre_csr(gedge_t* edges, index_t nedges)
{
	#pragma omp parallel 
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
    
    _beg_pos = (index_t*) calloc(sizeof(index_t), vert_count + 1);  
	#ifndef HALF_GRID
	_beg_pos_in = (index_t*) calloc(sizeof(index_t), vert_count + 1);  
	#endif

	vert_degree = (degree_t*)calloc(sizeof(degree_t), vert_count);

    //calculate the degree and beg_pos
    pre_csr(edges, nedges);

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

void gConv::proc_csr_dir(string idir, string part_file)
{
    struct dirent *ptr;
    DIR *dir;
    int file_count = 0;
    string edgefile;
    
    double start = mywtime();
    
    _beg_pos = (index_t*) calloc(sizeof(index_t), vert_count + 1);  
    #ifndef HALF_GRID
    _beg_pos_in = (index_t*) calloc(sizeof(index_t), vert_count + 1);  
    #endif

    vert_degree = (degree_t*)calloc(sizeof(degree_t), vert_count);
    
    //Read graph files
    dir = opendir(idir.c_str());
    while (NULL != (ptr = readdir(dir))) {
        if (ptr->d_name[0] == '.') continue;
        
        file_count++;
        edgefile = idir + "/" + string(ptr->d_name);
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
        //---------------calculate the degree and beg_pos-------------
        edges = (gedge_t*) malloc(st_edge.st_size);
        FILE* f = fopen(edgefile.c_str(), "rb");
        assert(f != 0);
        fread(edges, sizeof(gedge_t), nedges, f);
        
        #pragma omp parallel 
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
        
        free(edges);
        close(fid_edge); 
        fclose(f);
        
    }
    closedir(dir);
    
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
    cout << file_count << endl;
    
    //----------- Make CSR file -------------------
    index_t* count_edge = (index_t*) calloc(sizeof(index_t), vert_count);  
    _adj = (vertex_t*)calloc(sizeof(vertex_t), _beg_pos[vert_count]);
    
	#ifndef HALF_GRID
	uint32_t* _adj_in = (uint32_t*)calloc(sizeof(uint32_t), _beg_pos_in[vert_count]);
    index_t* count_edge_in = (index_t*) calloc(sizeof(index_t), vert_count);  
	#endif

    //---classify the edges in the grid
    //Read graph files
    dir = opendir(idir.c_str());
    while (NULL != (ptr = readdir(dir))) {
        if (ptr->d_name[0] == '.') continue;
        
        file_count++;
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
        
        //calculate the degree and beg_pos
        edges = (gedge_t*) malloc(st_edge.st_size);
        FILE* f = fopen(edgefile.c_str(), "rb");
        assert(f != 0);
        fread(edges, sizeof(gedge_t), nedges, f);

        #pragma omp parallel  
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
        free(edges);
        fclose(f);
        close(fid_edge);
    }
    closedir(dir);
   
    //munmap (edges, st_edge.st_size);
	//free(_adj);
	//free(_beg_pos);
	//#ifndef HALF_GRID
	//free(_adj_in);	
	//#endif

    double end = mywtime();
    cout << "CSR conversion Time = " << end -  start << endl;
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
	
    /*
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
    */
}
