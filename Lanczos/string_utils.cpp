#include "string_utils.hpp"
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

// trim from start (in place)
void string_utils::ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
void string_utils::rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}
