//========================================================================
// Project     : DiversityCaller
// Name        : EnumParser.hpp
// Author      : Xiao Yang
// Created on  : Apr 11, 2012
// Version     : 1.0
// Copyright   : The Broad Institute
//  				 SOFTWARE COPYRIGHT NOTICE AGREEMENT
// 				 This software and its documentation are copyright (2012)
//				 by the Broad Institute. All rights are reserved.
//
// 				 This software is supplied without any warranty or 
//				 guaranteed support whatsoever. The Broad Institute cannot 
//				 be responsible for its use,	misuse, or functionality.
// Description :
//========================================================================


#ifndef ENUMPARSER_HPP_
#define ENUMPARSER_HPP_

/*
 * Example usage:
 *		EnumParser<DiNt> parser;
 *  		std::cout << parser.ParseSomeEnum("TA");
 *		Can define new Enum as needed.
 */

template <typename T>
class EnumParser
{
private:
    std::map <std::string, T> enummap;
public:
    EnumParser();

    T parseEnum(const std::string& value)
    {
    		typename std::map <std::string, T>::iterator it = enummap.find(value);
        if (it  == enummap.end()) {
        		std::cout << "[Err] parseEnum: can't find value: " + value << "\n";
        		exit(1);
        }
        return it->second;
    }
};

/* di-nucleotide */
enum DiNt{ AA, AC, AG, AT,
		   CA, CC, CG, CT,
		   GA, GC, GG, GT,
		   TA, TC, TG, TT,
		   // insertion
		   /*IA, IC, IG, IT,
		   AI, CI, GI, TI,
		   // deletion
		   DA, DC, DG, DT,
		   AD, CD, GD, TD,*/
		   // first base
		   A, C, G, T
} ;

#define DiNtSize 20

template<> inline EnumParser<DiNt>::EnumParser()
{
	// Dinucleotide
    enummap["AA"] = AA; enummap["AC"] = AC; enummap["AG"] = AG; enummap["AT"] = AT;
    enummap["CA"] = CA; enummap["CC"] = CC; enummap["CG"] = CG; enummap["CT"] = CT;
    enummap["GA"] = GA; enummap["GC"] = GC; enummap["GG"] = GG; enummap["GT"] = GT;
    enummap["TA"] = TA; enummap["TC"] = TC; enummap["TG"] = TG; enummap["TT"] = TT;

    /*
    // Insertion
    enummap["IA"] = IA; enummap["IC"] = IC; enummap["IG"] = IG; enummap["IT"] = IT;
    enummap["AI"] = AI; enummap["CI"] = CI; enummap["GI"] = GI; enummap["TI"] = TI;

    // deletion
    enummap["DA"] = DA; enummap["DC"] = DC; enummap["DG"] = DG; enummap["DT"] = DT;
    enummap["AD"] = AD; enummap["CD"] = CD; enummap["GD"] = GD; enummap["TD"] = TD;
	*/
    // first base
    enummap["A"] = A; enummap["C"] = C; enummap["G"] = G; enummap["T"] = T;

}

#endif /* ENUMPARSER_HPP_ */
