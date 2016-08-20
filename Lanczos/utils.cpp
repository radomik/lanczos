#include "utils.hpp"

const bool_index BOOLS[] = {
	bool_index("true", true), bool_index("false", false),
	bool_index("yes", true), bool_index("no", false)
};

int arg2bool(bool& retval, int argc, char** argv, int argIndex) {
	if (argIndex >= argc) {
		fprintf(stderr, "%s: Invalid argument index %d (argc=%d)\n", __FUNCTION__, argIndex, argc);
		return -1;
	}
	
	for (size_t i = 0; i < array_len(BOOLS); i++) {
		if (! strcmp(argv[argIndex], BOOLS[i].str)) {
			retval = BOOLS[i].value; return 0;
		}
	}
	
	fprintf(stderr, "Invalid argument '%s' at index %d (expected one of: {", argv[argIndex], argIndex);
	for (size_t i = 0; i < array_len(BOOLS); i++) {
		fprintf(stderr, "%s'%s'", i?", ":"", BOOLS[i].str);
	}
	fprintf(stderr, "})\n");
	return -1;
}
