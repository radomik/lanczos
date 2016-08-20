#include "file_utils.hpp"

const std::string file_utils::EXT_CSV(".csv");
const std::string file_utils::EXT_BIN(".bin");
const std::string file_utils::PATH_SEP("/");

CsvBinPaths::CsvBinPaths(const std::string& work_dir, const std::string& filename_without_ext) {
	m_csv = work_dir;
	m_csv += file_utils::PATH_SEP;
	m_csv += filename_without_ext;
	m_bin = m_csv;
	m_csv += file_utils::EXT_CSV;
	m_bin += file_utils::EXT_BIN;
}

int serialization::saveBin(const std::string& path, double value) {
	FILE* f = fopen(path.c_str(), "wb");
	if (! f) {
		fprintf(stderr, "%s: Error opening file '%s'\n", __FUNCTION__, strerror(errno));
		return -1;
	}
	if (fwrite(&value, sizeof(value), 1, f) != 1) {
		fprintf(stderr, "%s: Error writing file '%s'\n", __FUNCTION__, strerror(errno));
		fclose(f);
		return -1;
	}
	fclose(f);
	fprintf(stderr, "%s: Saved file '%s'\n", __FUNCTION__, path.c_str());
	return 0;
}

int serialization::readBin(const std::string& path, double& value) {
	FILE* f = fopen(path.c_str(), "rb");
	if (! f) {
		fprintf(stderr, "%s: Error opening file '%s'\n", __FUNCTION__, strerror(errno));
		return -1;
	}
	if (fread(&value, sizeof(value), 1, f) != 1) {
		fprintf(stderr, "%s: Error reading file '%s'\n", __FUNCTION__, strerror(errno));
		fclose(f);
		return -1;
	}
	fclose(f);
	return 0;
}

int serialization::saveTxt(const std::string& path, double value) {
	FILE* f = fopen(path.c_str(), "wb");
	if (! f) {
		fprintf(stderr, "%s: Error opening file '%s'\n", __FUNCTION__, strerror(errno));
		return -1;
	}
	if (fprintf(f, "%lf", value) < 0) {
		fprintf(stderr, "%s: Error writing file '%s'\n", __FUNCTION__, strerror(errno));
		fclose(f);
		return -1;
	}
	fclose(f);
	fprintf(stderr, "%s: Saved file '%s'\n", __FUNCTION__, path.c_str());
	return 0;
}

int serialization::readTxt(const std::string& path, double& value) {
	FILE* f = fopen(path.c_str(), "rb");
	if (! f) {
		fprintf(stderr, "%s: Error opening file '%s'\n", __FUNCTION__, strerror(errno));
		return -1;
	}
	if (fscanf(f, "%lf", &value) < 0) {
		fprintf(stderr, "%s: Error reading file '%s'\n", __FUNCTION__, strerror(errno));
		fclose(f);
		return -1;
	}
	fclose(f);
	return 0;
}
