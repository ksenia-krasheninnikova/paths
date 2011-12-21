#include <iostream>
#include <queue>
#include <cmath>
#include <algorithm>
#include "graph.h"

	using std::cin;
	using std::cout;
	using std::endl;
	using std::vector;
	using std::pair;


const double EPSILON = 0.02; 
//const double INF = 100000;

vector < int > &Graph::operator[] (int i)
{
	return adjList[i];
}

const vector < int > &Graph::operator[] (int i) const
{
	return adjList[i];
}
void Graph::initialize( vector <double> &masses, const vector< pair <double,char> > &amino_acids)
{
	N = masses.size();
	cout << N << endl;
	
	adjList.resize(N);
	sort(masses.begin(),masses.end());
    	ParentMass = masses[masses.size() - 1];
	
	visited.resize(N);
	
	visited.assign(N, false);

	for (int i = 0; i < masses.size(); i++)
		cout << masses[i] << " ";
	cout << endl;

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			if (in(amino_acids, fabs(masses[i] - masses[j])))
			{
				adjList[i].push_back(j);				
			}
			
		}
	}

	for (int i = 0; i < N; i++)
	{
		cout << i << ": ";
		for (int j = 0; j < adjList[i].size(); j++)
		{
			cout << adjList[i][j] << " ";
		}
		cout << endl;
	}	

}


bool Graph::in(const vector< pair <double,char> > &amino_acids, double k)
{
	for (int i = 0; i < amino_acids.size(); i++)
	{
		if (fabs(amino_acids[i].first - k) < EPSILON) return true;	
	}

	return false;
}

bool Graph::in(const vector <int> &path, int k)
{
	for (int i = 0; i < path.size(); i++)
	{
		if (path[i] == k) return true;
	}

	return false;
}


void Graph::pathStarter(const vector <double> &masses, const int& lenBoundDown, const int& lenBoundUp)
{
	vector <int> path;
	
	for (int i = 0; i < N; i++)
	{ 
		if (visited[i])
		{
			continue;
		}
		path.push_back(i);
		getPath(i, masses, path, lenBoundDown, lenBoundUp);
		path.pop_back();
	}

	
}

double Graph::pathWeight(const vector <int>& path, const vector <double> &masses)
{
	double res = 0.0;
	for (int i = 1; i < path.size(); i++)
		res += masses[path[i]] - masses[path[i - 1]];

	return res;

}

void Graph::getAminoAcidSequences(vector< pair <double,char> > &amino_acids, const vector <double> &masses)
// pathSymbols[i][j] соответствует номеру аминокислоты или какой-то группе аминокислотных остатков в списке аминокислот
{
 
	pathSymbols.resize(pathsInMasses.size());

	for (int i = 0; i < pathsInMasses.size(); i++)
	{
		for (int j = 0; j < pathsInMasses[i].size(); j++)
		{
			for (int k = 0; k < amino_acids.size(); k++)
			{

				if (j == 0 && fabs(amino_acids[k].first - masses[pathsInMasses[i][j]]) < EPSILON)
				{
					pathSymbols[i].push_back(k);
				}
				else if ((j == 0)&&(k == amino_acids.size() - 1))
				{
					//double val = fabs(masses[pathsInMasses[i][j - 1]] - masses[pathsInMasses[i][j]]);
					double val = masses[pathsInMasses[i][j]];
					amino_acids.push_back(pair <double, char> (val,'?'));
					pathSymbols[i].push_back(k);
					cout << k << " - " << amino_acids[k].second << " ";								
				}
				else if (fabs(amino_acids[k].first - fabs(masses[pathsInMasses[i][j - 1]] - masses[pathsInMasses[i][j]])) < EPSILON) 
				{
					pathSymbols[i].push_back(k);
					cout << k << " - " << amino_acids[k].second << " ";
					
					break;
				}
				
			}
		}
		cout << endl;		
	}
}


void Graph::getPath(const int& v, const vector <double> &masses, vector <int>& path, const int& lenBoundDown, const int& lenBoundUp)
{
	visited[v] = true;

	if (path.size() > lenBoundUp)
	{
		return;
	}

	if ((path.size() >= lenBoundDown) )//&& fabs(masses[v] - ParentMass) < EPSILON) 
	{
		for (int q = 0; q < path.size(); q++)
	
			cout << path[q] << " ";
		
		cout << endl;
		
		pathsInMasses.push_back(path);
		
	}
	
		
	for (int i = 0; i < adjList[v].size(); i++)
	{
			path.push_back(adjList[v][i]);
			getPath(adjList[v][i], masses, path, lenBoundDown, lenBoundUp);
			path.pop_back();
	}
}

