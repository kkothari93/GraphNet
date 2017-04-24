#include <iostream>
#include "input.h"


#define DIM 2

using namespace std;

//TODO: Have to add support for other types of elements

inline int get_num_vertices(int elem_type){
	switch(elem_type){
		case 1: return 2;
		case 2: return 3;
		case 3: return 4;
		case 4: return 4;
		case 5: return 8;
		case 7: return 5;
		default: return -1;
	}

}

inline bool does_file_exist(string& fname){
	ifstream infile(fname);
	return infile.good();
}

void mapping(int& edge_counter, int elem_type, int n_vertices, stringstream& input, int* edges){
	int* local_nodes = new int[n_vertices];
	int garbage[3];
	input>>garbage[0]>>garbage[1]>>garbage[2];
	for(int n = 0; n<n_vertices; n++){
		input>>local_nodes[n];
	}
	switch(elem_type){
		case 2: 
		// [0,1], [0,2], [1, 2]
		edges[edge_counter*2] = local_nodes[0] -1 ;
		edges[edge_counter*2+1] = local_nodes[1] - 1 ;
		edges[edge_counter*2+2] = local_nodes[0] - 1 ;
		edges[edge_counter*2+3] = local_nodes[2] - 1 ;
		edges[edge_counter*2+4] = local_nodes[1] - 1 ;
		edges[edge_counter*2+5] = local_nodes[2] - 1 ;
		edge_counter += 3; 
		//TODO: write mapping code

	}
	delete[] local_nodes;
}

void mapping(int& edge_counter, int elem_type){
	switch(elem_type){
		case 2:
			edge_counter += 3;
		//TODO: add other mappings
	}

}

void read_n(int& n_nodes, int& n_elems, string& fname){

	string line;
	ifstream source;
	stringstream in_pos;
	source.open(fname);
	bool read_num_nodes = false, read_n_elems = false;
	
	n_nodes = 0;
	n_elems = 0;
	do{
		getline(source, line);
		if(line.find("$Nodes") != string::npos){
			getline(source, line);
			n_nodes = stoi(line);
			cout<<"The number of nodes are "<<n_nodes<<endl;
			read_num_nodes = true;
		}
	}while(!read_num_nodes);
	do{
		getline(source, line);
		if(line.find("$Elements") != string::npos){
			int id, elem_type;
			int num_vertices, num_lines;
			getline(source, line);
			num_lines = stoi(line);
			int c = 0;
			for(int i=0; i<num_lines; i++){
				getline(source, line);
				//cout<<line<<endl;
				in_pos<<line;
				in_pos>>id;
				in_pos>>elem_type;
				//cout<<"element type is "<<elem_type<<endl;
				num_vertices = get_num_vertices(elem_type);
				if(num_vertices>0){
					mapping(n_elems, elem_type);
				}
				in_pos.str(std::string());
				in_pos.clear();
			}
			read_n_elems = true;
			}
			//TODO: get elements
	}while(!read_n_elems);
	source.close();
}

void take_input(float* R, int* edges, int n_nodes, int n_elems, string& fname) {
	
	string line;
	ifstream source;
	stringstream in_pos;
	source.open(fname);

	bool read_num_nodes = false, read_n_elems = false;

	do{
		getline(source, line);
		if(line.find("$Nodes") != string::npos){
			getline(source, line);
			float r[3]; int id;
			for(int i=0; i<n_nodes; i++){
				getline(source, line);
				in_pos<<line;

				in_pos>>id>>r[0]>>r[1]>>r[2];
				for(int d=0; d<DIM; d++){
					R[(id-1)*DIM + d] = r[d];
				}
				in_pos.str(std::string());
				in_pos.clear();
			}
			read_num_nodes = true;
		}

	}while(!read_num_nodes);

	do{
		getline(source, line);
		//cout<<line<<"\n";
		if(line.find("$Elements") != string::npos){
			int id, elem_type;
			int num_vertices, num_lines;
			getline(source, line);
			num_lines = stoi(line);
			int c = 0;
			for(int i=0; i<num_lines; i++){
				getline(source, line);
				//cout<<line<<endl;
				in_pos<<line;
				in_pos>>id;
				in_pos>>elem_type;
				//cout<<"element type is "<<elem_type<<endl;
				num_vertices = get_num_vertices(elem_type);
				if(num_vertices>0){
					mapping(c, elem_type, num_vertices, in_pos, edges);
				}
				in_pos.str(std::string());
				in_pos.clear();
			}
			read_n_elems = true;
			if(c==n_elems){
				cout<<"Everything's alright here!\n";
			}
			else{
				cout<<"Damn! Something's not right here...\n";
				cout<<"Read "<<n_elems;
				cout<<", Got "<<c<<endl;
			}
			}
			//TODO: get elements
	}while(!read_n_elems);
	source.close();
}
