#include "Lanczno.hpp"
#include "TrMatrixd.hpp"
#include "Matrixd.hpp"
#include "Vecd.hpp"
#include "utils.hpp"
#include "file_utils.hpp"
#include "mrrr.hpp"

struct Flags {
	bool save_bin;
	bool save_csv;
	bool dbg;
	
	Flags() : save_bin(false), save_csv(false), dbg(false) { }
	
	int read(int argc, char** argv) {
		if (arg2bool(save_bin, argc, argv, argc-3)) return 1;
		if (arg2bool(save_csv, argc, argv, argc-2)) return 1;
		if (arg2bool(dbg, argc, argv, argc-1)) return 1;
		return 0;
	}
};

static void main_lanczno(const double* A, size_t m, size_t n, const Vecd* startvec) {
	double* ritz = new double[m];
	double* S = new double[m*m];
	int ret = lanczno(ritz, S, A, n, m, startvec);
	fprintf(stderr, "%s: lanczno() exited with status %d\n", __FUNCTION__, ret);

	delete[] ritz;
	delete[] S;
}

static int mat_csv2bin(int argc, char** argv) {
	if (argc != 4) {
		printf("Usage: %s %s <input CSV file> <output binary file>\n", argv[0], argv[1]);
		return 1;
	}

	size_t n;
	double* A = TrMatrixd::readCsv(argv[2], ",", &n);
	if (! A) {
		fprintf(stderr, "%s: Error reading input csv file\n", __FUNCTION__);
		return 1;
	}
	
	if (TrMatrixd::saveBin(argv[3], A, n) != 0) {
		fprintf(stderr, "%s: Error saving matrix to bin file\n", __FUNCTION__);
		delete[] A;
		return 1;
	}

	delete[] A;
	return 0;
}

static int vec_csv2bin(int argc, char** argv) {
	if (argc != 4) {
		printf("Usage: %s %s <input CSV file> <output binary file>\n", argv[0], argv[1]);
		return 1;
	}

	Vecd v;
	if (v.readCsv(argv[2], true) != 0) {
		fprintf(stderr, "%s: Error reading input csv file\n", __FUNCTION__);
		return 1;
	}
	
	if (v.saveBin(argv[3]) != 0) {
		fprintf(stderr, "%s: Error saving vector to bin file\n", __FUNCTION__);
		return 1;
	}

	return 0;
}

static int mat_rdbin(int argc, char** argv) {
	if (argc != 3) {
		printf("Usage: %s %s <input binary file>\n", argv[0], argv[1]);
		return 1;
	}
	
	size_t n;
	double* A = TrMatrixd::readBin(argv[2], &n);
	
	if (! A) {
		fprintf(stderr, "%s: Error reading matrix from bin file\n", __FUNCTION__);
		return 1;
	}
	
	TrMatrixd::print(A, n, stdout, "A");
	delete[] A;
	return 0;
}

static int vec_rdbin(int argc, char** argv) {
	if (argc != 3) {
		printf("Usage: %s %s <input binary file>\n", argv[0], argv[1]);
		return 1;
	}
	
	Vecd v;
	if (v.readBin(argv[2]) != 0) {
		fprintf(stderr, "%s: Error reading vector from bin file\n", __FUNCTION__);
		return 1;
	}
	
	v.print(stdout, "v");
	return 0;
}

void lancno_init_work(
	const double* A, const Vecd& startvec, size_t n, size_t m, Vecd& a, Vecd& b, double& anorm, bool dbg) {
	
	Vecd v(n);
	Vecd v2(n);
	Vecd vt(n);
	Vecd r(n);
	Vecd rt(n);

	v.set(startvec);
	a.zero();
	b.zero();

	if (dbg) {
		v.print(stderr, "v");
		a.print(stderr, "a");
		b.print(stderr, "b");
	}

	// initial Lanczos iteration, m steps
	for (uint32_t k = 0; ; ) {
		if (dbg) {
			fprintf(stderr, "for: k == %u\n", k);
		}
		if (k == 0) {
			TrMatrixd::mulByVector(r, A, v);  // r{n} = a{n}{n} * v{n}
		} else {
			TrMatrixd::mulByVector(rt, A, v);     // rt{n} = a{n}{n} * v{n}
			Vecd::mulByScalar(vt, v2, b[k-1]); // vt{n} = b[k-1]{1} * v2{n}
			Vecd::sub(r, rt, vt);              // r{n}  = rt{n} - vt{n}
		}

		if (dbg) {
			r.print(stderr, "r");
		}
		a[k] = Vecd::dotProduct(v, r);  // a[k]{1} = SUM v[i]*r[i]
		if (dbg) {
			a.print(stderr, "a");
		}

		Vecd::mulByScalar(vt, v, a[k]); // vt{n}   = a[k]{1} * v{n}
		Vecd::sub(r, r, vt);            // r{n}    = r{n} - vt{n}
		if (dbg) {
			r.print(stderr, "r");
		}
		b[k] = r.getNorm();             // b[k]{1} = |r{n}|
		if (dbg) {
			b.print(stderr, "b");
		}

		// estimate |A|_2 by |T|_1
		if (k == 0) {
			anorm = fabs(a[0] + b[0]);
		} else {
			anorm = std::max(anorm, b[k-1]+fabs(a[k])+b[k]);
		}

		if (dbg) {
			fprintf(stderr, "anorm = %f\n\n", anorm);
		}

		if (++k == m) {
			break;
		}

		// prepare next step, k = {1 ... m-1}
		v2.set(v);                          // v2{n} = v{n}
		if (dbg) {
			v2.print(stderr, "v2");
		}
		Vecd::divByScalar(v, r, b[k-1]); // EACH_i v[i] = r[i]/b[k-1]
		if (dbg) {
			v.print(stderr, "v");
		}
	}
}

int lancno_init(int argc, char** argv) {
	if (argc != 12) {
		printf(
	"Usage: %s %s work_dir m a_name startvec_name a_vec_prx b_vec_prx anorm_prx save_bin save_csv dbg\n\n\
For example:\n\
	%s %s \"../data\" 200 \"A_n-48_g-10\" \"startvec_n-48\" \"a\" \"b\" \"anorm\" true true false\n",
		argv[0], argv[1], argv[0], argv[1]);
		return 1;
	}
	
	// read flags
	Flags f;
	if (f.read(argc, argv) != 0) return 1;
	
	// read Lanczos iteration count (m)
	uint32_t m;
	if ((sscanf(argv[3], "%u", &m) != 1) || (m <= 1)) {
		fprintf(stderr, "%s: m must be greater than 1\n", __FUNCTION__);
		return 1;
	}
	
	// working directory path
	const std::string work_dir = argv[2];
	
	// starting vector
	Vecd startvec;
	if (startvec.readCsvOrBin(work_dir, argv[5])) {
		fprintf(stderr, "%s: Error reading startvec vector\n", __FUNCTION__);
		return 1;
	}
	
	// A matrix
	double* A;
	size_t n;
	if (! (A=TrMatrixd::readCsvOrBin(work_dir, argv[4], &n))) {
		fprintf(stderr, "%s: Error reading A matrix\n", __FUNCTION__);
		return 1;
	}
	
	if (startvec.size() != n) {
		fprintf(stderr, "%s: A matrix (%zu) and startvec (%zu) size mismatch\n", __FUNCTION__, n, startvec.size());
		delete[] A;
		return 1;
	}
	
	// output vectors
	Vecd a(m), b(m);
	
	// output anorm
	double anorm;
	
	lancno_init_work(A, startvec, n, m, a, b, anorm, f.dbg);
	delete[] A;
	
	char suffix[128];
	snprintf(suffix, sizeof(suffix), "_n-%zu_m-%u", n, m);
	
	if (f.save_bin) {
		const std::string out_a_path = work_dir + file_utils::PATH_SEP + argv[6] + suffix + file_utils::EXT_BIN;
		const std::string out_b_path = work_dir + file_utils::PATH_SEP + argv[7] + suffix + file_utils::EXT_BIN;
		const std::string out_anorm_path = work_dir + file_utils::PATH_SEP + argv[8] + suffix + file_utils::EXT_BIN;
		
		a.saveBin(out_a_path.c_str());
		b.saveBin(out_b_path.c_str());
		serialization::saveBin(out_anorm_path, anorm);
	}
	if (f.save_csv) {
		const std::string out_a_path = work_dir + file_utils::PATH_SEP + argv[6] + suffix + file_utils::EXT_CSV;
		const std::string out_b_path = work_dir + file_utils::PATH_SEP + argv[7] + suffix + file_utils::EXT_CSV;
		const std::string out_anorm_path = work_dir + file_utils::PATH_SEP + argv[8] + suffix + file_utils::EXT_CSV;
		
		a.saveCsv(out_a_path.c_str());
		b.saveCsv(out_b_path.c_str());
		serialization::saveTxt(out_anorm_path, anorm);
	}
	
	return 0;
}

int main_mrrr(int argc, char** argv) {
	if (argc != 10) {
		printf("Usage: %s %s work_dir a_vec_name b_vec_name s_prx ritz_prx save_bin save_csv dbg\n\n\
For example:\n\
	%s %s \"../data\" \"a_n-48_m-100\" \"b_n-48_m-100\" \"S\" \"ritz\" true true false\n",
		argv[0], argv[1], argv[0], argv[1]);
		return 1;
	}
	
	// read flags
	Flags f;
	if (f.read(argc, argv) != 0) return 1;
	
	const std::string work_dir = argv[2];
	
	Vecd a;
	if (a.readCsvOrBin(work_dir, argv[3])) {
		fprintf(stderr, "%s: Error reading a vector\n", __FUNCTION__);
		return 1;
	}
	
	Vecd b;
	if (b.readCsvOrBin(work_dir, argv[4])) {
		fprintf(stderr, "%s: Error reading b vector\n", __FUNCTION__);
		return 1;
	}
	
	if (a.size() != b.size()) {
		fprintf(stderr, "%s: a vector (%zu) and b vector (%zu) size mismatch\n", __FUNCTION__, a.size(), b.size());
		return 1;
	}
	
	size_t m = a.size();
	double* S = new double[m*m];
	Vecd ritz(m);
	
	int ret = mrrr(ritz.data(), S, a, b, f.dbg);
	if (ret) {
		fprintf(stderr, "%s: mrrr() failed with status %d\n", __FUNCTION__, ret);
		delete[] S;
		return 1;
	}
	
	char suffix[128];
	snprintf(suffix, sizeof(suffix), "_m-%zu", m);
	
	if (f.save_bin) {
		const std::string out_s_path = work_dir + file_utils::PATH_SEP + argv[5] + suffix + file_utils::EXT_BIN;
		const std::string out_ritz_path = work_dir + file_utils::PATH_SEP + argv[6] + suffix + file_utils::EXT_BIN;
		
		Matrixd::saveBin(out_s_path.c_str(), S, m);
		ritz.saveBin(out_ritz_path.c_str());
	}
	if (f.save_csv) {
		const std::string out_s_path = work_dir + file_utils::PATH_SEP + argv[5] + suffix + file_utils::EXT_CSV;
		const std::string out_ritz_path = work_dir + file_utils::PATH_SEP + argv[6] + suffix + file_utils::EXT_CSV;
		
		Matrixd::saveCsv(out_s_path.c_str(), S, m);
		ritz.saveCsv(out_ritz_path.c_str());
	}
	
	delete[] S;
	return 0;
}

int cmp_vec(int argc, char** argv) {
	if (argc != 6) {
		printf("Usage: %s %s work_dir x_vec_name y_vec_name eps\n\n\
For example:\n\
	%s %s \"../data\" \"a_n-48_m-100\" \"b_n-48_m-100\" 0.0001\n",
		argv[0], argv[1], argv[0], argv[1]);
		return 1;
	}
	
	double eps;
	if (sscanf(argv[5], "%lf", &eps) != 1) {
		fprintf(stderr, "%s: Error reading eps argument\n", __FUNCTION__);
		return 1;
	}
	
	float eps_f = (float)eps;
	fprintf(stderr, "%s: eps: %.20f\n", __FUNCTION__, eps_f);
	
	const std::string work_dir = argv[2];
	
	Vecd x;
	if (x.readCsvOrBin(work_dir, argv[3])) {
		fprintf(stderr, "%s: Error reading x vector\n", __FUNCTION__);
		return 1;
	}
	
	Vecd y;
	if (y.readCsvOrBin(work_dir, argv[4])) {
		fprintf(stderr, "%s: Error reading y vector\n", __FUNCTION__);
		return 1;
	}
	
	size_t first_diff_index;
	float first_diff;
	bool equal = x.equals(y, eps_f, &first_diff_index, &first_diff);
	
	fprintf(stderr, "%s: a.size: %zu, b.size: %zu, equal: %s, first_diff_index = %zu, first_diff = %.10f\n", 
		__FUNCTION__, x.size(), y.size(), sbool(equal), first_diff_index, first_diff);
	return 0;
}

int cmp_mat(int argc, char** argv) {
	if (argc != 6) {
		printf("Usage: %s %s work_dir x_mat_name y_mat_name eps\n\n\
For example:\n\
	%s %s \"../data\" \"a_n-48_m-100\" \"b_n-48_m-100\" 0.0001\n",
		argv[0], argv[1], argv[0], argv[1]);
		return 1;
	}
	
	double eps;
	if (sscanf(argv[5], "%lf", &eps) != 1) {
		fprintf(stderr, "%s: Error reading eps argument\n", __FUNCTION__);
		return 1;
	}
	
	float eps_f = (float)eps;
	fprintf(stderr, "%s: eps: %.20f\n", __FUNCTION__, eps_f);
	
	const std::string work_dir = argv[2];
	
	double *X, *Y;
	size_t Xn, Yn;
	if (! (X=TrMatrixd::readCsvOrBin(work_dir, argv[3], &Xn))) {
		fprintf(stderr, "%s: Error reading X matrix\n", __FUNCTION__);
		return 1;
	}
	
	if (! (Y=TrMatrixd::readCsvOrBin(work_dir, argv[4], &Yn))) {
		fprintf(stderr, "%s: Error reading Y matrix\n", __FUNCTION__);
		delete[] X;
		return 1;
	}
	
	size_t first_diff_index;
	float first_diff;
	bool equal = TrMatrixd::equals(X, Xn, Y, Yn, eps_f, &first_diff_index, &first_diff);
	
	fprintf(stderr, "%s: X.size: %zu, Y.size: %zu, equal: %s, first_diff_index = %zu, first_diff = %.10f\n", 
		__FUNCTION__, Xn, Yn, sbool(equal), first_diff_index, first_diff);
		
	delete[] X;
	delete[] Y;
	return 0;
}

int main_rescon(int argc, char** argv) {
	if (argc != 17) {
		printf(
	"Usage: %s \
%s \
work_dir \
s_name \
ritz_name \
b_name \
anorm_name \
eps_name \
s_prx \
ritz_prx \
lres_prx \
cul_prx \
idx_prx \
k_prx \
save_bin \
save_csv \
dbg\
\n\n\
For example:\n\
	%s %s \"../data\" \"c-S_m-200\" \"c-ritz_m-200\" \"c-b_n-48_m-200\" \"c-anorm_n-48_m-200\" \
\"c-rescon-eps\" \"c-rescon-S\" \"c-rescon-ritz\" \"c-rescon-lres\" \"c-rescon-cul\" \"c-rescon-idx\" \
\"c-rescon-k\" true true false\n",
		argv[0], argv[1], argv[0], argv[1]);
		return 1;
	}
	
	// read flags
	Flags f;
	if (f.read(argc, argv) != 0) return 1;
	
	const std::string work_dir = argv[2];
	
	double* S;
	size_t Sm
	if (! (S=Matrixd::readCsvOrBin(work_dir, argv[3], &Sm))) {
		fprintf(stderr, "%s: Error reading S matrix\n", __FUNCTION__);
		return 1;
	}
	
	Vecd a;
	if (a.readCsvOrBin(work_dir, argv[3])) {
		fprintf(stderr, "%s: Error reading a vector\n", __FUNCTION__);
		return 1;
	}
	
	Vecd b;
	if (b.readCsvOrBin(work_dir, argv[4])) {
		fprintf(stderr, "%s: Error reading b vector\n", __FUNCTION__);
		return 1;
	}
	
	if (a.size() != b.size()) {
		fprintf(stderr, "%s: a vector (%zu) and b vector (%zu) size mismatch\n", __FUNCTION__, a.size(), b.size());
		return 1;
	}
	
	size_t m = a.size();
	
	double* S = new double[m*m];
	Vecd ritz(m);
	
	int ret = mrrr(ritz.data(), S, a, b, f.dbg);
	if (ret) {
		fprintf(stderr, "%s: mrrr() failed with status %d\n", __FUNCTION__, ret);
		delete[] S;
		return 1;
	}
	
	char suffix[128];
	snprintf(suffix, sizeof(suffix), "_m-%zu", m);
	
	if (f.save_bin) {
		const std::string out_s_path = work_dir + file_utils::PATH_SEP + argv[5] + suffix + file_utils::EXT_BIN;
		const std::string out_ritz_path = work_dir + file_utils::PATH_SEP + argv[6] + suffix + file_utils::EXT_BIN;
		
		Matrixd::saveBin(out_s_path.c_str(), S, m);
		ritz.saveBin(out_ritz_path.c_str());
	}
	if (f.save_csv) {
		const std::string out_s_path = work_dir + file_utils::PATH_SEP + argv[5] + suffix + file_utils::EXT_CSV;
		const std::string out_ritz_path = work_dir + file_utils::PATH_SEP + argv[6] + suffix + file_utils::EXT_CSV;
		
		Matrixd::saveCsv(out_s_path.c_str(), S, m);
		ritz.saveCsv(out_ritz_path.c_str());
	}
	
	delete[] S;
	return 0;
}

int main(int argc, char** argv) {
	if (argc >= 2) {
		if (! strcmp(argv[1], "mat_csv2bin")) {
			return mat_csv2bin(argc, argv);
		}
		if (! strcmp(argv[1], "vec_csv2bin")) {
			return vec_csv2bin(argc, argv);
		}
		if (! strcmp(argv[1], "mat_rdbin")) {
			return mat_rdbin(argc, argv);
		}
		if (! strcmp(argv[1], "vec_rdbin")) {
			return vec_rdbin(argc, argv);
		}
		if (! strcmp(argv[1], "lancno_init")) {
			return lancno_init(argc, argv);
		}
		if (! strcmp(argv[1], "mrrr")) {
			return main_mrrr(argc, argv);
		}
		if (! strcmp(argv[1], "cmp_vec")) {
			return cmp_vec(argc, argv);
		}
		if (! strcmp(argv[1], "cmp_mat")) {
			return cmp_mat(argc, argv);
		}
		if (! strcmp(argv[1], "rescon")) {
			return main_rescon(argc, argv);
		}
	}

	if (argc < 3) {
		printf("Usage: %s <m> <input weight matrix CSV file> [<initial vector>]\n", argv[0]);
		return 0;
	}

	uint32_t m;

	if ((sscanf(argv[1], "%u", &m) != 1) || (m <= 1)) {
		fprintf(stderr, "%s: First parameter must be integer above or equal 2\n", __FUNCTION__);
		return 1;
	}

	fprintf(stderr, "%s: m: %u\n", __FUNCTION__, m);

	size_t  n;
	double* A = TrMatrixd::readCsv(argv[2], ",", &n);

	if (! A) {
		fprintf(stderr, "%s: Error reading weight matrix from '%s'\n", __FUNCTION__, argv[2]);
		return 1;
	}

	fprintf(stderr, "%s: Weight matrix size (n): %zu\n", __FUNCTION__, n);

	TrMatrixd::print(A, n, stderr, "A");

	if (argc > 3) {
		Vecd startvec(n);
		if (startvec.readCsv(argv[3]) != 0) {
			fprintf(stderr, "%s: Error reading start vector from file '%s'\n", __FUNCTION__, argv[3]);
			delete[] A;
			return 1;
		}

		fprintf(stderr, "%s: Using predefinied start vector\n", __FUNCTION__);
		startvec.print(stderr, "startvec");
		main_lanczno(A, m, n, &startvec);
	}
	else {
		fprintf(stderr, "%s: Using random start vector\n", __FUNCTION__);
		main_lanczno(A, m, n, NULL);
	}

	delete[] A;

	return 0;
}
