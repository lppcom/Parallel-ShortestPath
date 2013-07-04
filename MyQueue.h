/* MyQueue.h -- header file for Class Queue and other data structure */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <vector>

using namespace std;

typedef struct Adjacent{
	int total;			// the total number of the neighbors
	int prev;			// the prev node along the path
	double length;		// the shortest path length from the source
	vector<int> neigh;		// the neighbor list
	vector<double> weight;	// the weight list
} Adj;

extern int source;
extern int vertex;
extern int commRank;
extern int commSize;
extern int *addr;

#ifndef QUEUE_H
#define QUEUE_H

#define UNSIGN  	-1		// unsigned vertex in the queue
#define UNDISCOVER	10000		// the infinite number in this one
#define SOURCE	-1		// the source of the path tree

#define QUEUE_EMPTY 	-1		// dequeue return -1 if the queue is empty
#define QUEUE_FRESH	-2		// dequeue return -2 if the queue is out of use

class MyQueue
{
public:
	int total;			// the total elements in the heap
	vector<int> nodes;		// the nodes in the heap
	vector<double> length;	// the shortest path from source
	vector<int> index;		// the position of one node in this heap

	MyQueue(vector<Adj>&);
	~MyQueue();
	int DeQueue(vector<Adj>&, double*);
	int Update(vector<Adj>&, int, int, double);
	void Insert(vector<Adj>&, int, double);
	int UpperMove(int);
	void LowerMove(int);
};

#endif
