#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <map>
#include <limits>

#include "parser.h"
#include "config.h"


extern const double EPSILON;

using parser::MassOfAcid;

class Graph
{

	std::vector< std::vector< int > > adjList; // the vertices achievable from each vertex (achievable -> dif of masses = 1 amino acid)
	std::vector< std::vector< int > > pathsInMasses;
	std::vector< std::vector< double > > foundSequences;
	std::vector < std::vector <int> > pathSymbols;

	std::vector <bool> visited;

	double ParentMass;
	int N;

	bool in(const std::vector< MassOfAcid > &amino_acids, double k, double eps);	
	bool in(const std::vector <int> &path, int k);

	double sum(const int& i, const std::vector <double> &masses);
	void getPath(const int& v, const std::vector <double> &masses, std::vector <int>& path, const int& lenBoundDown, const int& lenBoundUp);
	void chooseSets(const std::vector <double> &masses);
	public:

	Graph(double ParentMass=0) : ParentMass(ParentMass), pathsInMasses(0) // N - the number of peaks
	{
	}

	const std::vector < int > &operator[] (int i) const;
	std::vector < int > &operator[] (int i);



	void initialize( std::vector <double> &masses, const std::vector< MassOfAcid > &amino_acids);

	void getAminoAcidSequences(std::vector< MassOfAcid > &amino_acids, const std::vector <double> &masses);


	double pathWeight(const std::vector <int>& path, const std::vector <double> &masses);
	void pathStarter(const std::vector <double> &masses, const int& lenBoundDown = 0, const int& lenBoundUp = std::numeric_limits<double>::max());

	void initialize(const std::vector <double> &masses);


};

#endif /* GRAPH_H */

