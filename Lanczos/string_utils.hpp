#ifndef STRING_UTILS_HPP
#define STRING_UTILS_HPP

#include <string>

/// http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
namespace string_utils {
	// trim from start (in place)
	void ltrim(std::string &s);

	// trim from end (in place)
	void rtrim(std::string &s);

	// trim from both ends (in place)
	inline void trim(std::string &s) {
		ltrim(s);
		rtrim(s);
	}

	// trim from start (copying)
	inline std::string ltrimmed(std::string s) {
		ltrim(s);
		return s;
	}

	// trim from end (copying)
	inline std::string rtrimmed(std::string s) {
		rtrim(s);
		return s;
	}

	// trim from both ends (copying)
	inline std::string trimmed(std::string s) {
		trim(s);
		return s;
	}
};

#endif
