#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include "MyQueue.h"

using namespace std;

#define 	MAXMSG 		100000		// the largest block for one message
#define  	UNREACH		-1			// the prev field is UNREACH, means this node is unreachable
#define 	FAR			1000		// the initial distance from the source to current node
#define		STOP		-1			// the stop signal in distribution phase

/* Read graph file into one machine and keep them in the data structure */
void ReadGraph(char *file);

/* Proc 0 send its data to other ones */
void DataTransmit();

/* Dijkstra algorithm in parallel */
void Dijkstra();

/* Write results into the file */
void WriteResults(char *file);

/* global variables */
int commSize;				// the total number of the machine
int commRank;				// the rank of this one machine
int vertex;					// the number of vertex in the graph
int maxneigh;				// the maximum number of one vertex's neighbor in the graph
int source;					// the source node for the shortest path
int *addr;					// mark the location of each vertex
vector<Adj> graph;			// the vector of the neighbor list structure
MyQueue *queue;				// the path heap(queue) for each machine

/* All of the test func */
void Test_graph();
void Test_queue();

int main(int argc, char **argv){

	// Analysis those arguments
	if(argc < 3){
		cout << "Error: arguements <input graph file> <output res file> and <source node>" << endl;
		exit(-1);
	}

	if(argc < 4){
		cout << "Please input the single source for shortest path" << endl;
		cin >> source;
	}
	else
		source = atoi(argv[3]);

	/* Initialize MPI */
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

	/* Root proc read the graph file */
	if(commRank == 0) ReadGraph(argv[1]);
	//if(commRank == 0) printf("All done\n");
	/* Root proc distribute those data to other one uniformly */
	DataTransmit();

	/* Create shortest path heap for each processor */
	queue = new MyQueue(graph);
//	if(commRank == 0) Test_graph();

	/* Dijkstra algorithm */
	Dijkstra();

	/* Output the result */
	WriteResults(argv[2]);

	MPI_Finalize();
}

void ReadGraph(char *file){

	/* open the file */
	fstream in_file;
	in_file.open(file, fstream::in);
	if(!in_file.is_open()){
		cout << "The file " << file << " could not be open." << endl;
		exit(-1);
	}

	/* read the file */
	int from, to;
	double weight;
	string buffer;
	istringstream is;
	vertex = 0;
	while(!in_file.eof()){
		// Read one line from the file
		getline(in_file, buffer);
		if(buffer.size() < 1) continue;
		is.str(buffer);
		is >> to >> from >> weight;
		is.clear();
//		cout << to << " " << from << " " << weight << endl;

		// manage the Adj vector
		if(vertex < to + 1){
			vertex = to + 1;
			graph.resize(vertex);
		}
		if(vertex < from + 1){
			vertex = from + 1;
			graph.resize(vertex);
		}

		// manage one edge
		graph[from].neigh.push_back(to);
		graph[from].weight.push_back(1.0);
	}

	// manage the graph vector
	maxneigh = 0;
	for(int i=0; i<vertex; i++){
		graph[i].total = graph[i].neigh.size();
		graph[i].prev = UNREACH;
		graph[i].length = FAR;
		maxneigh = maxneigh < graph[i].total ? graph[i].total : maxneigh;
	}

	// manage the address vector
	addr = (int*)calloc(vertex, sizeof(int));

//	Test_graph();
}

void DataTransmit(){

	MPI_Bcast(&vertex, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&maxneigh, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(commRank != 0){
		graph.resize(vertex);
		for(int i=0; i<vertex; i++)
			graph[i].total = 0;
		addr = (int*)malloc(sizeof(int) * vertex);
	}
//	printf("[proc %d] vertex %d\n", commRank, vertex);
	if(commRank == 0){
		int *sendbuf = (int*)malloc(sizeof(int) * MAXMSG);
//		double *senddbuf = (double*)malloc(sizeof(double) * MAXMSG);
		int bufptr;
		int tag;

		for(int i=1; i<commSize; i++){
			// Assign vertices according to randomly, uniform block
			int begin = vertex / commSize * i;
			int end = begin + vertex / commSize;
			if(i == commSize - 1) end = vertex;

			// Initialize the buffers
			bufptr = 0;
			tag = 0;
			for(int j=begin; j<end; j++){
				// One integer section consists of three kinds of values
				sendbuf[bufptr++] = j;					// One is the vertex ID
				sendbuf[bufptr++] = graph[j].total;		// Two is the # of edges
				for(int k=0; k<graph[j].total; k++)		// Three is the neighbor ID
					sendbuf[bufptr++] = graph[j].neigh[k];
				addr[j] = i;							// update the addr of vertex

				// If the buffer is almost full, we have to send data to corresponding machine
				if(bufptr > MAXMSG - maxneigh * 2){
					MPI_Send(sendbuf, bufptr, MPI_INT, i, tag++, MPI_COMM_WORLD);
					bufptr = 0;
//					printf("sending tag %d done\n", tag);
				}
			}
			sendbuf[bufptr++] = STOP;
			MPI_Send(sendbuf, bufptr, MPI_INT, i, tag++, MPI_COMM_WORLD);
//			printf("sending tag %d done\n", tag);
		}
	}
	else{
		int *recvbuf = (int*)malloc(sizeof(int) * MAXMSG);
		int bufptr;
		int tag = 0;
		MPI_Status status;
		int length;
		int terminal = 0;

		while(!terminal){	/* until receiving terminal bit */
			MPI_Recv(recvbuf, MAXMSG, MPI_INT, 0, tag++, MPI_COMM_WORLD, &status);

			MPI_Get_count(&status, MPI_INT, &length);
//			printf("recv tag %d done, length %d\n", tag, length);
			bufptr = 0;
			while(bufptr < length){		// analyse the buffer
				int vert = recvbuf[bufptr++];
				if(vert == STOP){
					terminal = 1;
					break;
				}
				int total = recvbuf[bufptr++];
				graph[vert].total = total;
				for(int i=0; i<total; i++)
					graph[vert].neigh.push_back(recvbuf[bufptr++]);
			}
		}
	}
	MPI_Bcast(addr, vertex, MPI_INT, 0, MPI_COMM_WORLD);
}


void Dijkstra(){

	int i, j;

	/* MPI_Alltoallv buffers */
	int *sendbuf = (int*) malloc(sizeof(int) * vertex * commSize);
	int *recvbuf = (int*) malloc(sizeof(int) * vertex * commSize);
	int *sendcnts = (int*) malloc(sizeof(int) * commSize);
	int *recvcnts = (int*) malloc(sizeof(int) * commSize);
	int *sendoff = (int*) malloc(sizeof(int) * commSize);
	int *recvoff = (int*) malloc(sizeof(int) * commSize);

	/* Elements in those array should be fixed */
	for(int i=0; i<commSize; i++){
		recvcnts[i] = vertex;
		recvoff[i] = vertex * i;
		sendoff[i] = vertex * i;
	}

	/* Elements in these array may change during process */
	for(int i=0; i<commSize; i++){
		sendbuf[sendoff[i]] = 0;
		sendcnts[i] = 1;
	}

	int cont, cont_all;					// the signal for checking whether it should terminate
	int position, num;
	int round = 0;
	int start = 0;
	cont = 1;
	while(cont && round < 453526){
	
		// monitor the situation of other processors
		MPI_Alltoallv(sendbuf, sendcnts, sendoff, MPI_INT, recvbuf, recvcnts, recvoff, MPI_INT, MPI_COMM_WORLD);

		// according result from each processor to update the heap
		for(int i=0; i<commSize; i++){
			if(i == commRank) continue;
			position = recvoff[i];
			num = recvbuf[position];
			int prev = recvbuf[position+1];
			position += 2;
			for(int j=0; j<num; j++){
				queue->Update(graph, recvbuf[position + j*2], prev, (double)recvbuf[position + j*2 + 1]);
			}
		}

		// After finishing the update from other ones, each machine should find its own shortest neighbor so far
		for(int i=0; i<commSize; i++){
			sendcnts[i] = 2;
			sendbuf[sendoff[i]] = 0;
		}

		// Dequeue to get the cloest neighbor
		double distance;
		int neigh_node, update_res, district;
		int close_node = queue->DeQueue(graph, &distance);
//		if(round >= 450000) printf("extract %d ", close_node);
//		printf("[proc %d], close one %d\n", commRank, close_node);
		if(close_node == QUEUE_FRESH){
//			printf("fresh here?");
			if(start){
//				printf("here?");
				cont = 0;
			}	
			else cont = 1;
		}
		else if(close_node == QUEUE_EMPTY)
			cont = 0;
		else if(close_node >= 0){
			num = graph[close_node].total;
			start = 1;
			for(int i=0; i<num; i++){
				neigh_node = graph[close_node].neigh[i];
				update_res = queue->Update(graph, neigh_node, close_node, distance+1);

				if(update_res){
					cont = 1;
					district = addr[neigh_node];
					if(district != commRank){
						sendbuf[sendoff[district] + sendcnts[district]] = neigh_node;
						sendcnts[district]++;
						sendbuf[sendoff[district] + sendcnts[district]] = (int) (distance+1);
						sendcnts[district]++;
						sendbuf[sendoff[district]] += 2;
					}
				}
			}

			for(int i=0; i<commSize; i++)
				sendbuf[sendoff[i]+1] = close_node;
		}
//		Test_queue();
		round++;
//		if(commRank == 0) printf("after round %d: %d\n", round, cont);
		MPI_Allreduce(&cont, &cont_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		cont = cont_all;
		
	}
}

void WriteResults(char *file){

	if(commRank == 0){
		int *recvbuf = (int*) malloc (sizeof(int) * vertex * 2);
		MPI_Status status;
		int length;
		int cur_node;

		for(int i=1; i<commSize; i++){
			MPI_Recv(recvbuf, vertex*2, MPI_INT, i, i, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_INT, &length);

			if(length % 3 != 0){
				cout << "ERROR in result" << endl;
				exit(-1);
			}

			for(int j=0; j<length; j+=3){
				cur_node = recvbuf[j];
				graph[cur_node].prev = recvbuf[j+1];
				graph[cur_node].length = (double)recvbuf[j+2];
			}
		}
	}
	else{
		int *sendbuf = (int*)malloc(sizeof(int) * vertex * 2);
		int bufptr = 0;

		for(int i=0; i<vertex; i++){
			if(addr[i] != commRank) continue;
			sendbuf[bufptr++] = i;
			sendbuf[bufptr++] = graph[i].prev;
			sendbuf[bufptr++] = (int)graph[i].length;
		}

		MPI_Send(sendbuf, bufptr, MPI_INT, 0, commRank, MPI_COMM_WORLD);
	}

	if(commRank == 0){
		fstream out_file;
		out_file.open(file, fstream::out);
		if(!out_file.is_open()) cout << file << " could not open" << endl;
		for(int i=0; i<vertex; i++){
			if(graph[i].length == UNDISCOVER)
				out_file << i << " inf" << endl;
			else
				out_file << i << " " << (int)graph[i].length << endl;
		}
		out_file.close();
	}
}

void Test_graph(){

	for(int i=0; i<graph.size(); i++){
		cout << "node " << i << " has " << graph[i].total << " neighbors: ";
		for(int j=0; j<graph[i].total; j++)
			cout << graph[i].neigh[j] << " ";
		cout << endl;
	}
}

void Test_queue(){
	printf("[Proc %d], queue: ", commRank);
	for(int i=0; i<queue->total; i++){
		printf("[%d %.1f] ", queue->nodes[i], queue->length[i]);
	}
	printf("\n");
}
