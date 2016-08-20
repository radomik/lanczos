#ifndef FILE_UTILS_HPP
#define FILE_UTILS_HPP

#include "utils.hpp"

class file_utils {
public:
	static const std::string EXT_CSV;
	static const std::string EXT_BIN;
	static const std::string PATH_SEP;
};

namespace serialization {
	int saveBin(const std::string& path, double value);
	int readBin(const std::string& path, double& value);
	int saveTxt(const std::string& path, double value);
	int readTxt(const std::string& path, double& value);
};

class CsvBinPaths {
public:
	CsvBinPaths(const std::string& work_dir, const std::string& filename_without_ext);
	
	const std::string& strCsv() const { return m_csv; }
	const std::string& strBin() const { return m_bin; }

	const char* csv() const { return m_csv.c_str(); }
	const char* bin() const { return m_bin.c_str(); }
	
private:
	std::string m_csv;
	std::string m_bin;
};

#endif
