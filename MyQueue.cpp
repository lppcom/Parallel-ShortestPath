/* MyQueue.cpp -- the MyQueue class */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "MyQueue.h"

MyQueue::MyQueue(vector<Adj> &neigh){
	
	// In this assignment, we assign those vertices uniformly
	int begin = vertex / commSize * commRank;
	int end = begin + vertex / commSize;
	if(commRank == commSize - 1) end = vertex;
	
	total = end - begin;
	nodes.resize(vertex);
	length.resize(vertex);
	index.resize(vertex);
	
	for(int i=0; i<vertex; i++)
		index[i] = UNSIGN;
	
	// Initialize the path heap
	for(int i=0; i<total; i++){
		nodes[i] = begin + i;
		index[begin + i] = i;
		length[i] = UNDISCOVER;
	}
	
	// Initialize the Adjacent vector
	for(int i=0; i<vertex; i++){
		neigh[i].prev = UNSIGN;
		neigh[i].length = UNDISCOVER;
	}
	
	// Assign the source vertex, update the heap
	if(addr[source] == commRank)
		Update(neigh, source, SOURCE, 0);
/*	
	if(commRank == 2){
		for(int i=0; i<vertex; i++){
			printf("[vertex %d] has %d neighbors: ", i, neigh[i].total);
			for(int j=0; j<neigh[i].total; j++){
				printf("%d ", neigh[i].neigh[j]);
			}
			printf("\n");
		}
	}
	if(commRank == 0){
		printf("begin %d end %d\n", begin, end);
		for(int i=begin; i<end; i++){
			printf("%d ", neigh[i].total);
			neigh[i].total = 0;
		}
	}
*/	
	return;
}

MyQueue::~MyQueue(){

	return;
}


int MyQueue::DeQueue(vector<Adj> &neigh, double *dist_v){

	if(total <= 0) return QUEUE_EMPTY;
	
	*dist_v = length[0];
	if(*dist_v >= UNDISCOVER - 1) return QUEUE_FRESH;
	
	int res = nodes[0];
	neigh[res].length = *dist_v;
	
	// change the last item to the top of the heap
	length[0] = length[total-1];		// length value to the top
	nodes[0] = nodes[total-1];			// last vertex to the top
	index[nodes[0]] = 0;				// change the index
	index[res] = -1;					// this vertex doesn't exist in heap
	total--;							// one vertex less
	
	LowerMove(0);
	return res;
}

int MyQueue::Update(vector<Adj> &neigh, int node, int prev, double dist_v){

	int position = index[node];			// the index of the node
	
	if(position < 0){
		if(neigh[node].length <= dist_v)
			return 0;
		Insert(neigh, node, dist_v);
		return 1;
	}
	
	if(length[position] <= dist_v)
		return 0;
		
	length[position] = dist_v;
	neigh[node].length = dist_v;
	neigh[node].prev = prev;
	
	int cont = 1;
	while(position > 0 && cont){
		cont = UpperMove(position);
		position /= 2;
	}
	return 1;
}

void MyQueue::Insert(vector<Adj> &neigh, int node, double dist_v){

	int position = total;
	nodes[position] = node;
	length[position] = dist_v;
	index[node] = position;
	total++;
	
	// Maintain the heap
	int cont = 1;
	while(position > 0 && cont){
		cont = UpperMove(position);
		position /= 2;
	}
	return;
}

int MyQueue::UpperMove(int loc){

	if(loc == 0) return 0;
	int upper = loc / 2;
	
	if(length[upper] <= length[loc]) return 0;
	
	double temp1 = length[loc];
	length[loc] = length[upper];
	length[upper] = temp1;
	
	int node1 = nodes[loc];
	int node2 = nodes[upper];
	int temp2 = nodes[loc];
	nodes[loc] = nodes[upper];
	nodes[upper] = temp2;
	
	index[node1] = upper;
	index[node2] = loc;
	
	return 1;
}

void MyQueue::LowerMove(int loc){

	int left_child = loc*2;
	int right_child = loc*2+1;
	
	int least_node;
	if(left_child < total && length[left_child] < length[loc])
		least_node = left_child;
	else
		least_node = loc;
		
	if(right_child < total && length[right_child] < length[least_node])
		least_node = right_child;
		
	if(least_node != loc){
		double temp1 = length[loc];
		length[loc] = length[least_node];
		length[least_node] = temp1;
		
		int node1 = nodes[loc];
		int node2 = nodes[least_node];
		int temp2 = nodes[loc];
		nodes[loc] = nodes[least_node];
		nodes[least_node] = temp2;
		
		index[node1] = least_node;
		index[node2] = loc;
		
		LowerMove(least_node);
	}
	return;
}
