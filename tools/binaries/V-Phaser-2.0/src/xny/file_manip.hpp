//============================================================================
// Project     : Diversifier
// Name        : file_manip.hpp
// Author      : Xiao Yang
// Created on  : Aug 22, 2011
// Version     :
// Copyright   :
// Description :
//============================================================================


#ifndef FILE_MANIP_HPP_
#define FILE_MANIP_HPP_

#include <iostream>
#include <string>
#include <dirent.h>
#include <fstream>

namespace xny{

	template <typename Output> bool dir_list(const std::string& name, Output out) {
		DIR* dir = opendir(name.c_str());
		if (dir == 0) return false;

		dirent* entry = readdir(dir);
		while (entry != 0) {
		  *(out++) = entry->d_name;
		  entry = readdir(dir);
		}
		closedir(dir);
	  return true;
	} // dir_list

	/*
	 * find the last/first [delim] in [str] if [dir] is set to be true/false
	 */
	inline std::string get_suffix (const std::string& str, bool last,
			const std::string& delim) {
		std::string suf = str;
		int pos = suf.find_last_of(delim);
		if (!last) pos = suf.find_first_of (delim);
		if (pos != (int) std::string::npos)
			return suf.substr(pos + 1, (int) suf.length() - pos);
		else return "";
	}

	/*
	 * find the first/last [delim] in [str] if [dir] is set to be true/false
	 */
	inline std::string get_prefix (const std::string& str, bool first,
			const std::string& delim) {
		std::string prfx = str;
		int pos = prfx.find_last_of(delim);
		if (first) pos = prfx.find_first_of(delim);
		if (pos!= (int) std::string::npos)
			return prfx.substr(0, pos);
		else return "";
	}

	inline void openfile (std::ofstream& outfile, const std::string& fileName){
		if (!outfile.is_open()) {
			outfile.open(fileName.c_str());
			//outfile.open(fileName.c_str(), std::fstream::app);
			if (! outfile.good() ){
				std::cout << "[WARNING] opening file: " << fileName << " failed!\n"
						  << "Will write to standard output ... \n";
			}
		}
	}

	inline void closefile(std::ofstream& outfile){
		if (outfile.good())	outfile.close();
	}


	/* @brief	Retrieve files [filelist] from directory [iDirNm] such
	 * 			that each file has suffix [suffix]
	 *
	 */
	inline void getInputFilelist(std::vector<std::string>& filelist,
			const std::string& file_suffix, const std::string& iDirNm) {

		if (iDirNm.empty()) {
			std::cout << "\n\nError: no input directory specified !!!\n";
			exit(1);
		} else {
			std::vector<std::string> fileNms;
			if (xny::dir_list(iDirNm, std::back_inserter(fileNms)) == false) {
					std::cout << "can't read from dir: " << iDirNm << "\n";
					exit(1);
			} else {
				for (int i = 0; i < (int) fileNms.size(); ++i) {

					if (fileNms[i].compare(".") == 0 ||
							fileNms[i].compare("..") == 0) continue;

					std::string suf = xny::get_suffix(fileNms[i], true, ".");

					if (suf.length() != file_suffix.length()) continue;

					/* convert to lower case and compare */
					bool isEqual = true;
					for (int j = 0; j < (int) suf.length(); ++ j) {
						if (std::tolower(suf.at(j)) !=
								std::tolower(file_suffix.at(j))){
							isEqual = false;
							break;
						}
					}
					if (isEqual) {
						filelist.push_back(iDirNm + "/" + fileNms[i]);
					}
				}
			} // if
		} // else

	} // getInputFilelist

} // namespace xny

#endif /* FILE_MANIP_HPP_ */
