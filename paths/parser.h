#ifndef PARSER_H
#define PARSER_H

#include <iostream>
#include <map>
#include <vector>
#include <string>	 
	 
#include "config.h"

extern const double EPS_PPM;
extern const double AMINO_MAX;

namespace parser
{
	//public:
	typedef std::pair <double,char> MassOfAcid;
	void getAminoAcids(const std::string& fileName, std::vector< MassOfAcid >& aminoTable);
	void getExperimentalMasses(const std::string& fileName, std::vector< double >& masses, std::vector < double > &intensity);
	void getExperimentalMassesSpectra(const std::string& fileName, std::vector< double >& masses, double &ParentMass, std::vector < double > &intensity);
};

#endif /* PARSER_H */
