#include <iostream>

#include "config.h"

#include "graph.h"
#include "parser.h"



using std::cin;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::pair;

int main()
{
	//Parser P;
	string name = "input.txt";
	vector< parser::MassOfAcid > aminoTable(AMINO_MAX);
	parser::getAminoAcids(name,  aminoTable);

/*	for (int i = 0; i < AMINO_MAX; i++)
	{
		cout << aminoTable[i].second << " " << aminoTable[i].first << endl;
	}
*/
	//string nameMasses = "input_masses.txt";
	string nameMasses = "scan_592.env";
	vector < double > masses;
	vector < double > intensities;
	double ParentMass;
	parser::getExperimentalMassesSpectra(nameMasses, masses, ParentMass, intensities);
		

	/* for (int i = 0; i < (int)masses.size(); i++)
	{
		cout << masses[i] << endl;
	}
	cout << endl;*/

	Graph G(ParentMass);
	G.initialize(masses, aminoTable);

	int lenBoundDown(5), lenBoundUp(masses.size());
	G.pathStarter(masses, lenBoundDown, lenBoundUp);
	G.getAminoAcidSequences(aminoTable, masses);
	
	
	return 0;
}
