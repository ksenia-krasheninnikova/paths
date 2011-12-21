#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include "parser.h"

const double EPS_PPM = 0.000001 * 2.5;

using std::cin;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::pair;

const double  AMINO_MAX = 20;

void parser::getAminoAcids(const string& fileName, vector< MassOfAcid >& aminoTable)
{
	FILE *pinput = fopen (fileName.c_str() , "r");
   	if (pinput == NULL) perror ("Error opening file");
	else
	{
		for (int i = 0; i < AMINO_MAX && ! feof (pinput); i++)
		{
			char buffer[100];
			aminoTable[i].second = fgetc(pinput);
			if (fgets (buffer , 100 , pinput) != NULL)
				aminoTable[i].first = strtod(buffer, 0);
			else
				break;
		
		}
		fclose (pinput);
		return;
	}
}

void parser::getExperimentalMasses(const string& fileName, vector< double >& masses, vector <double>& intensity)
{
	int d;
	std::ifstream inFile;
	inFile.open(fileName.c_str());
	while (!inFile.eof())
	{
		string s;
		inFile >> s;
		
		masses.push_back(strtod(s.c_str(), 0));
		s.clear();
		inFile >> s;
		intensity.push_back(strtod(s.c_str(), 0));
		s.clear();
	}
	inFile.close();
	
}

void parser::getExperimentalMassesSpectra(const string& fileName, vector< double >& masses, double &ParentMass, vector <double>& intensity)
{
	int d;
	std::ifstream inFile;
	inFile.open(fileName.c_str());
	string s;
	size_t found;

	//while (!inFile.eof())
	//{
	do
	{
		inFile >> s;
		found = s.find("MONOISOTOPIC_MASS");

		if (found!=string::npos)
		{
			s.clear();
			inFile >> s;
			ParentMass = strtod(s.c_str(), 0);
			s.clear();		
			do
			{
				inFile >> s;
				found = s.find("CHARGE");	
			}while(found==string::npos);
			s.clear();
			inFile >> s;
			ParentMass = ParentMass * strtod(s.c_str(), 0);
			s.clear();
			break;

		}
	}while(found==string::npos);	
	
	while (!inFile.eof())
	{	
		inFile >> s;	
		found = s.find("THEO_MONO_MASS");	
		if (found!=string::npos)
		{	
			s.clear();
			inFile >> s;
			//double m =

			double val = strtod(s.c_str(), 0);
			double val_rev = ParentMass - strtod(s.c_str(), 0);
			bool exists = false, exists_rev = false;
			for (int i = 0; i < (int)masses.size(); i++)
			{
				if (fabs(masses[i] - val) < EPS_PPM*(masses[i] + val)) exists = true;
				if (fabs(masses[i] - val_rev) < EPS_PPM*(masses[i] + val_rev)) exists_rev = true;
			}
 
			if (!exists) masses.push_back(val);
			if (!exists_rev) masses.push_back(val_rev);
			s.clear();
			//cout << masses[masses.size() - 1] << " ";
			do
			{
				inFile >> s;
				found = s.find("THEO_INTE_SUM");	
			}while(found==string::npos);
			///// что делать с интенсивностью, если мы выяснили, что такая масса уже существует в векторе?
			s.clear();
			inFile >> s;
			intensity.push_back(strtod(s.c_str(), 0));
			intensity.push_back(strtod(s.c_str(), 0));
			s.clear();
			//cout << intensity[intensity.size() - 1] << endl;

		}
		
		
	}
	inFile.close();
}
