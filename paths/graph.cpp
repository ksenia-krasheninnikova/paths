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


const double EPSILON = 0.07; 
//const double EPS_PPM = 0.000001 * 2.5;
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
	
	adjList.resize(N);
	sort(masses.begin(),masses.end());
    	//ParentMass = masses[masses.size() - 1];
	
	visited.resize(N);
	
	std::fill(visited.begin(), visited.end(), false);

	#ifdef DEBUG
		for (int i = 0; i < masses.size(); i++)
			cout << masses[i] << " ";
		cout << endl;
	#endif

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			if (in(amino_acids, fabs(masses[i] - masses[j]),/*EPSILON*/(EPS_PPM) * (masses[i] + masses[j]) ))
			{
				adjList[i].push_back(j);				
			}
			
		}
	}

	#ifdef DEBUG
	for (int i = 0; i < N; i++)
		{
			cout << i << ": ";
			for (int j = 0; j < adjList[i].size(); j++)
			{
				cout << adjList[i][j] << " ";
			}
			cout << endl;
		}
	#endif	

}


bool Graph::in(const vector< pair <double,char> > &amino_acids, double k, double eps)
{
	for (int i = 0; i < amino_acids.size(); i++)
	{		
		//double eps = EPS_PPM * (amino_acids[i].first + k);
		//cout << eps << endl;
		
		if (fabs(amino_acids[i].first - k) < /*EPSILON*/eps) return true;	
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

	for (int i = 0; i < masses.size(); i++)
	{ 
		if (visited[i])
		{
			continue;
		}
		//visited[v] = true;
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
				//double eps = EPS_PPM*(amino_acids[k].first + masses[pathsInMasses[i][j]]);
				//cout << eps << endl;
				// a begining of a sequence that corresponds to a certain element in amino_acids
				if (j == 0 && fabs(amino_acids[k].first - masses[pathsInMasses[i][j]]) < EPSILON)
				{
					pathSymbols[i].push_back(k);

					if (k == 9) // I-acid
					{
						vector <int> a;
						a.push_back(k + 1);
						pathSymbols.push_back(a);
					}
					cout << amino_acids[k].second << " " << "(" << amino_acids[k].first << ") ";
					//break;
				}
				// a begining of a sequence that does not correspond to any of the elements in amino_acids
				else if ((j == 0)&&(k == amino_acids.size() - 1))
				{
					//double val = fabs(masses[pathsInMasses[i][j - 1]] - masses[pathsInMasses[i][j]]);
					double val = masses[pathsInMasses[i][j]];
					//amino_acids.push_back(pair <double, char> (val,'?'));
					//pathSymbols[i].push_back(k);
					cout << " ? (" << val << ") ";
				}
				// not a begining of a sequence and found in the amino_acids
				else 
				{
					double a = amino_acids[k].first;
					double b = fabs( masses[pathsInMasses[i][j - 1]] - masses[pathsInMasses[i][j]]);
					double eps = EPS_PPM*( masses[pathsInMasses[i][j - 1]] + masses[pathsInMasses[i][j]]);
					//eps = EPSILON;
					if (fabs(a - b) < /*EPSILON*/eps)
					{
						pathSymbols[i].push_back(k);
					
						if (k == 9) // I-acid
						{
							pathSymbols.push_back(pathSymbols[i]);
							pathSymbols[pathSymbols.size() - 1].push_back(k + 1);
						}

						cout << amino_acids[k].second << " " << "(" << amino_acids[k].first << ") ";
						break;
					}
				}
				
			}

			if (j == pathsInMasses[i].size() - 1) cout << "last mass: " << masses[pathsInMasses[i][j]];
		}
		cout << endl;		
	}
}

void Graph::getPath(const int& v, const vector <double> &masses, vector <int>& path, const int& lenBoundDown, const int& lenBoundUp)
{
	//if (visited[v]) return;

	visited[v] = true;

	if (path.size() > lenBoundUp)
	{
		return;
	}

	if ( path.size() >= lenBoundDown  )//&& fabs(masses[v] - ParentMass) < EPSILON) 
	{
		bool finish = true;
		for (int i = 0; i < adjList[v].size(); i++) if (!visited[adjList[v][i]]) finish = false;
		
		if (finish)
		{
			#ifdef DEBUG
				for (int q = 0; q < path.size(); q++)	
					cout << path[q] << "-" << masses [path[q] ] << " ";		
				cout << endl;		
			#endif
			pathsInMasses.push_back(path);
			return;
		}
		
	}
	
		
	for (int i = 0; i < adjList[v].size(); i++)
	{
			path.push_back(adjList[v][i]);
			getPath(adjList[v][i], masses, path, lenBoundDown, lenBoundUp);
			path.pop_back();
	}
}

