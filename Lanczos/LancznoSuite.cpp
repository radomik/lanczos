#include "utils.hpp"
#include "string_utils.hpp"
#include "file_utils.hpp"
#include <vector>
#include <fstream>
#include <algorithm>

typedef int32_t  s32;
typedef uint32_t u32;
typedef uint64_t u64;
typedef double flt;

using namespace std;

class Vecd : public vector<flt> {
public:
	Vecd() : vector<flt>() { }
	Vecd(size_t size) : vector<flt>(size) { }
	Vecd(const Vecd& src) : vector<flt>(src) { }

	void random() {
		random(size(), &front());
	}
	
	flt norm() {
		return norm(size(), &front());
	}
	
	void zero() {
		zero(size(), &front());
	}
	
	flt getNorm() const {
		return getNorm(size(), &front());
	}
	
	flt getNormPow() const {
		return getNormPow(size(), &front());
	}
	
	void print(FILE* file, const char* name) {
		print(size(), &front(), file, name);
	}
	
	void mulByScalar(flt s) {
		mulByScalar(size(), &front(), &front(), s);
	}
	
	void mulByScalar(const Vecd& v, flt s) {
		assert(size() == v.size());
		mulByScalar(size(), &front(), &v.front(), s);
	}
	
	void divByScalar(flt s) {
		divByScalar(size(), &front(), &front(), s);
	}
	
	void divByScalar(const Vecd& v, flt s) {
		assert(size() == v.size());
		divByScalar(size(), &front(), &v.front(), s);
	}
	
	void mulInvByScalar(flt s) {
		mulInvByScalar(size(), &front(), &front(), s);
	}
	
	void mulInvByScalar(const Vecd& v, flt s) {
		assert(size() == v.size());
		mulInvByScalar(size(), &front(), &v.front(), s);
	}
	
	void add(const Vecd& b) {
		assert(size() == b.size());
		add(size(), &front(), &front(), &b.front());
	}
	
	void add(const Vecd& a, const Vecd& b) {
		assert(a.size() == b.size());
		assert(size() == b.size());
		add(size(), &front(), &a.front(), &b.front());
	}
	
	void sub(const Vecd& b) {
		assert(size() == b.size());
		sub(size(), &front(), &front(), &b.front());
	}
	
	void sub(const Vecd& a, const Vecd& b) {
		assert(a.size() == b.size());
		assert(size() == b.size());
		sub(size(), &front(), &a.front(), &b.front());
	}
	
	flt dotProduct(const Vecd& b) const {
		assert(size() == b.size());
		return dotProduct(size(), &front(), &b.front());
	}
	
	static void zero(size_t size, flt* data) {
		for (size_t i = 0; i < size; i++) {
			data[i] = 0.0;
		}
	}
	
	static flt getNorm(size_t size, const flt* data) {
		return sqrt(getNormPow(size, data));
	}
	
	static void random(size_t size, flt* data) {
		assert(sizeof(data[0]) % sizeof(s32) == 0);
		srand(time(NULL));

		const s32* vi_end = (s32*)(data+size);
		for (s32* vi = (s32*)data; vi != vi_end; vi++) {
			*vi = rand();
		}
	}
	
	static flt norm(size_t size, flt* data) {
		const flt s = getNorm(size, data);
		for (size_t i = 0; i < size; i++) {
			data[i] /= s;
		}
		return s;
	}
	
	static flt getNormPow(size_t size, const flt* data) {
		flt s = 0.0;
		for (size_t i = 0; i < size; i++) {
			s += data[i]*data[i];
		}
		return s;
	}
	
	static void print(size_t size, const flt* data, FILE* file, const char* name) {
		fprintf(file, "%s = [\n", name);
		for (size_t i = 0; i < size; i++) {
			fprintf(file, "%s%.4f ", i ? " ;\n" : "", data[i]);
		}
		fprintf(file, " ]\n\n");
	}
	
	static void mulByScalar(size_t size, flt* data, const flt* src, flt s) {
		for (size_t i = 0; i < size; i++) {
			data[i] = src[i] * s;
		}
	}
	
	static void divByScalar(size_t size, flt* data, const flt* src, flt s) {
		for (size_t i = 0; i < size; i++) {
			data[i] = src[i] / s;
		}
	}
	
	static void mulInvByScalar(size_t size, flt* data, const flt* src, flt s) {
		for (size_t i = 0; i < size; i++) {
			data[i] *= 1.0 / (src[i] * s);
		}
	}
	
	static void add(size_t size, flt* data, const flt* a, const flt* b) {
		for (size_t i = 0; i < size; i++) {
			data[i] = a[i] + b[i];
		}
	}
	
	static void sub(size_t size, flt* data, const flt* a, const flt* b) {
		for (size_t i = 0; i < size; i++) {
			data[i] = a[i] - b[i];
		}
	}
	
	static flt dotProduct(size_t size, const flt* a, const flt* b) {
		flt s = 0.0;
		for (size_t i = 0; i < size; i++) {
			s += a[i] * b[i];
		}
		return s;
	}
};

class TrMatrixd : public Vecd {
public:
	TrMatrixd() : Vecd() { }
	TrMatrixd(size_t size) : Vecd(size) { }
	TrMatrixd(const TrMatrixd& src) : Vecd(src) { }
	
	//NOTE: size() returns nelem
	//NOTE: dim() returns n
	
	size_t dim() const { return m_dim; }
	
	int readCsv(const char* filename, const char* separator) {
		std::ifstream file(filename);

		if (! file) {
			fprintf(stderr, "%s: Failed to open '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
			return -1;
		}

		std::string delim = separator;
		std::string line;
		size_t  i = 0;
		size_t  di = 0;
		flt*    data = NULL;
		bool    error = false;
		m_dim = 0;

		while (std::getline(file, line)) {
			if (i == 0) {
				m_dim = std::count(line.begin(), line.end(), delim[0]);
				if (m_dim == 0) {
					break;
				}
				m_dim++;
				size_t datasize = sizeFromDim(m_dim);
				resize(datasize);
				fprintf(stderr, "%s: dim=%zu, size=%zu\n", __FUNCTION__, m_dim, size());
				data = &front();
			}

			size_t j = 0;
			size_t pos;
			while ((pos = line.find(delim)) != std::string::npos) {
				if (j >= i) {
					if (j >= m_dim) {
						fprintf(stderr, "%s: [1] Invalid column index %zu at row %zu\n", __FUNCTION__, j, i);
						error = true;
						break;
					}

					std::string token = line.substr(0, pos);
					string_utils::trim(token);

					if (sscanf(token.c_str(), "%lf", &data[di++]) != 1) {
						fprintf(stderr, "%s: [1] Cannot read value at row index: %zu and column index: %zu\n", __FUNCTION__, i, j);
						error = true;
						break;
					}
				}

				line.erase(0, pos + delim.length());
				j++;
			}

			if (error) {
				break;
			}

			string_utils::trim(line);
			if (! line.empty()) {
				if (j >= m_dim) {
					fprintf(stderr, "%s: [2] Invalid column index %zu at row %zu\n", __FUNCTION__, j, i);
					error = true;
					break;
				}

				if (sscanf(line.c_str(), "%lf", &data[di++]) != 1) {
					fprintf(stderr, "%s: [2] Cannot read value at row index: %zu and column index: %zu\n", __FUNCTION__, i, j);
					error = true;
					break;
				}
			}

			i++;
		}
		
		file.close();

		if (m_dim == 0) {
			fprintf(stderr, "%s: No data in file '%s'\n", __FUNCTION__, filename);
			return -1;
		}
		if (error) {
			return -1;
		}
		return 0;
	}
	
	int readBin(const char* filename) {
		FILE* bin = fopen(filename, "rb");
		if (! bin) {
			fprintf(stderr, "%s: Error opening input file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
			return -1;
		}
		
		uint32_t nu;
		size_t ret;
		if (((ret=fread(&nu, sizeof(nu), 1, bin)) != 1) || (nu == 0)) {
			fprintf(stderr, "%s: Error reading matrix size (dim=%u, ret=%zu) [%s]\n", __FUNCTION__, nu, ret, strerror(errno));
			fclose(bin);
			return -1;
		}
		
		size_t nelem = sizeFromDim(nu);
		
		fprintf(stderr, "%s: Reading matrix of dim=%u, nelem=%zu, sizeof(A[0])=%zu\n", __FUNCTION__, nu, nelem, sizeof((&front())[0]));
		
		resize(nelem);
		fprintf(stderr, "%s: Size of data is %zu\n", __FUNCTION__, size());
		
		flt* data = &front();
		
		if ((ret=fread(data, sizeof(data[0]), nelem, bin)) != nelem) {
			fprintf(stderr, "%s: Error reading matrix content (nelem=%zu, ret=%zu) [%s]\n", __FUNCTION__, nelem, ret, strerror(errno));
			fclose(bin);
			return -1;
		}
		
		fclose(bin);
		return 0;
	}
	
	int readCsvOrBin(const std::string& work_dir, const std::string& filename_without_ext) {
		CsvBinPaths p(work_dir, filename_without_ext);
		if (readBin(p.bin()) != 0) {
			if (readCsv(p.csv(), ",") != 0) {
				fprintf(stderr, "%s: Error reading matrix [csv=%s, bin=%s]\n", __FUNCTION__, p.csv(), p.bin());
				return -1;
			}
			saveBin(p.bin());
		}
		return 0;
	}
	
	void print(FILE* file, const char* name) const {
		print(&front(), dim(), file, name);
	}
	
	int saveCsv(const char* filename) const {
		return saveCsv(&front(), dim(), filename);
	}
	
	int saveBin(const char* filename) const {
		return saveBin(&front(), dim(), filename);
	}
	
	bool equals(const TrMatrixd& other, flt eps, size_t* first_diff_index, flt* first_diff) const {
		return equals(&front(), dim(), &other.front(), other.dim(), eps, first_diff_index, first_diff);
	}
	
	static size_t sizeFromDim(size_t datadim) {
		return (datadim*(datadim+1))>>1;
	}
	
	static void mulByVector(Vecd& res, const TrMatrixd& A, const Vecd& b) {
		assert(A.size() == b.size());
		assert(res.size() == b.size());
		mulByVector(&res.front(), &A.front(), &b.front(), b.size());
	}
	
	static void mulByVector(flt* res, const flt* A, const flt* b, size_t n) {
		const size_t t = (n<<1) + 1;
		const size_t u = n - 1;

		// loop is optimized so each iteration is independent and may be parallelized
		for (size_t i = 0; i < n; i++) {
			size_t ai = ((t - i) * i)>>1, bi = 0;
			flt res_i = 0.0;

			if (i != 0) { // elements left from diagonal
				size_t step = u;
				size_t sta_curr = i;
				for ( ; sta_curr < ai; ) {
					res_i += b[bi++] * A[sta_curr];
					sta_curr += step;
					step--;
				}
			}

			// elements right of and including diagonal
			const size_t ai_end = ai + (n-i);
			for ( ; ai < ai_end; ) {
				res_i += b[bi++] * A[ai++];
			}

			res[i] = res_i;
		}
	}
	
	static void print(const flt* A, size_t dim, FILE* file, const char* name) {
		const size_t t = (dim<<1) + 1;
		const size_t u = dim - 1;

		fprintf(file, "%s = [\n", name);
		for (size_t i = 0; i < dim; i++) {
			size_t ai = ((t - i) * i)>>1;

			if (i != 0) { // elements left from diagonal
				size_t step = u;
				size_t sta_curr = i;

				fprintf(file, ";\n");

				for ( ; sta_curr < ai; ) {
					fprintf(file, "%15.10f ", A[sta_curr]);
					sta_curr += step;
					step--;
				}
			}

			// elements right of and including diagonal
			const size_t ai_end = ai + (dim-i);
			for ( ; ai < ai_end; ) {
				fprintf(file, "%15.10f ", A[ai++]);
			}
		}
		fprintf(file, "]\n");
	}
	
	static int saveCsv(const flt* A, size_t dim, const char* filename) {
		FILE* f = fopen(filename, "wb");
		if (! f) {
			fprintf(stderr, "%s: Error opening file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
			return -1;
		}
		
		const size_t t = (dim<<1) + 1;
		const size_t u = dim - 1;

		for (size_t i = 0; i < dim; i++) {
			size_t ai = ((t - i) * i)>>1;

			if (i != 0) { // elements left from diagonal
				size_t step = u;
				size_t sta_curr = i;

				fprintf(f, "\n");

				for ( ; sta_curr < ai; ) {
					fprintf(f, "%lf,", A[sta_curr]);
					sta_curr += step;
					step--;
				}
			}

			// elements right of and including diagonal
			const size_t ai_end = ai + (dim-i);
			for ( ; ; ) {
				fprintf(f, "%lf", A[ai]);
				if (++ai == ai_end) break;
				fprintf(f, ",");
			}
		}
		
		fclose(f);
		fprintf(stderr, "%s: Saved file '%s'\n", __FUNCTION__, filename);
		return 0;
	}
	
	static int saveBin(const flt* A, size_t dim, const char* filename) {
		if (dim > UINT_MAX) {
			fprintf(stderr, "%s: Error. Matrix dimension %zu is above UINT_MAX=%u\n", __FUNCTION__, dim, UINT_MAX);
			return -1;
		}

		uint32_t nu = (uint32_t)dim;
		size_t nelem = sizeFromDim(dim);

		fprintf(stderr, "%s: Saving matrix of dim=%u, nelem=%zu, sizeof(A[0])=%zu\n", __FUNCTION__, nu, nelem, sizeof(A[0]));

		FILE* bin = fopen(filename, "wb");
		if (! bin) {
			fprintf(stderr, "%s: Error opening output file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
			return -1;
		}
		
		int ret = (fwrite(&nu, sizeof(nu), 1, bin) != 1) || (fwrite(A, sizeof(A[0]), nelem, bin) != nelem);
		
		if (ret) {
			fprintf(stderr, "%s: Error writing file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
			ret = -1;
		}
		else {
			fprintf(stderr, "%s: Saved file '%s'\n", __FUNCTION__, filename);
		}
		
		fclose(bin);
		return ret;
	}
	
	static bool equals(const flt* X, size_t Xn, const flt* Y, size_t Yn, flt eps, size_t* first_diff_index, flt* first_diff) {
		if (!(X && Y)) {
			fprintf(stderr, "%s: Warning at least one matrix is null (X: %p, Y: %p)\n", __FUNCTION__, X, Y);
			return false;
		}
		if (Xn != Yn) {
			return false;
		}
		if (X == Y) {
			return true;
		}
		
		size_t nelem = sizeFromDim(Xn);
		
		for (size_t i = 0; i < nelem; i++) {
			const flt diff = fabsf(X[i]-Y[i]);
			if (diff > eps) {
				if (first_diff_index) {
					*first_diff_index = i;
				}
				if (first_diff) {
					*first_diff = diff;
				}
				return false;
			}
		}
		
		if (first_diff_index) {
			*first_diff_index = 0;
		}
		if (first_diff) {
			*first_diff = 0.0;
		}
		return true;
	}

	
	size_t m_dim;
};

class LancznoInit {
public:
	LancznoInit(u32 m) : m_avec(m), m_bvec(m) { }
	
	/**
	 * @param A{n.n}      linearized reduced symmetric square matrix used to compute eigen values/vectors for
	 * @param startvec{n} optional starting vector 
	 * @param dbg         whether print debug information
	 * 
	 * Produces results in m_avec, m_bvec, m_anorm
	 * 
	 */
	void run(const TrMatrixd& A, const Vecd& startvec, bool dbg) {
		u32 n = startvec.size();
		u32 m = m_avec.size();
		Vecd v(startvec);
		Vecd v2(n);
		Vecd vt(n);
		Vecd r(n);
		Vecd rt(n);
		
		m_avec.zero();
		m_bvec.zero();

		if (dbg) {
			v.print(stderr, "v");
			m_avec.print(stderr, "a");
			m_bvec.print(stderr, "b");
		}

		// initial Lanczos iteration, m steps
		for (uint32_t k = 0; ; ) {
			if (dbg) {
				fprintf(stderr, "for: k == %u\n", k);
			}
			if (k == 0) {
				TrMatrixd::mulByVector(r, A, v);  // r{n} = a{n}{n} * v{n}
			} else {
				TrMatrixd::mulByVector(rt, A, v); // rt{n} = a{n}{n} * v{n}
				vt.mulByScalar(v2, m_bvec[k-1]); // vt{n} = b[k-1]{1} * v2{n}
				r.sub(rt, vt);              // r{n}  = rt{n} - vt{n}
			}

			if (dbg) {
				r.print(stderr, "r");
			}
			m_avec[k] = v.dotProduct(r);  // a[k]{1} = SUM v[i]*r[i]
			if (dbg) {
				m_avec.print(stderr, "a");
			}

			vt.mulByScalar(v, m_avec[k]); // vt{n}   = a[k]{1} * v{n}
			r.sub(vt);               // r{n}    = r{n} - vt{n}
			if (dbg) {
				r.print(stderr, "r");
			}
			m_bvec[k] = r.getNorm();             // b[k]{1} = |r{n}|
			if (dbg) {
				m_bvec.print(stderr, "b");
			}

			// estimate |A|_2 by |T|_1
			if (k == 0) {
				m_anorm = fabs(m_avec[0] + m_bvec[0]);
			} else {
				m_anorm = std::max(m_anorm, m_bvec[k-1]+fabs(m_avec[k])+m_bvec[k]);
			}

			if (dbg) {
				fprintf(stderr, "anorm = %f\n\n", m_anorm);
			}

			if (++k == m) {
				break;
			}

			// prepare next step, k = {1 ... m-1}
			v2 = v;                          // v2{n} = v{n}
			if (dbg) {
				v2.print(stderr, "v2");
			}
			v.divByScalar(r, m_bvec[k-1]); // EACH_i v[i] = r[i]/b[k-1]
			if (dbg) {
				v.print(stderr, "v");
			}
		}
	}
	
	Vecd m_avec;
	Vecd m_bvec;
	flt  m_anorm;
};

class Mrrr {
public:
	Mrrr() { }
};

class LancznoSuite {
public:
	/**
	 * @param n input matrix size
	 * @param m Lanczos iteration count
	*/
	LancznoSuite(u32 n, u32 m) : m_n(n), m_m(m), m_startvec(n), m_init(m) { }
	
	/**
	 * @param A{n.n}      linearized reduced symmetric square matrix used to compute eigen values/vectors for
	 * @param startvec{n} optional starting vector 
	 * @param dbg         whether print debug information
	 */
	void run(const TrMatrixd& A, const Vecd* startvec, bool dbg) {
		if (startvec) {
			m_startvec = *startvec;
		} else {
			m_startvec.random();
			m_startvec.norm();
		}
		
		m_init.run(A, m_startvec, dbg);
		
		
	}
	
	u32 m_n;
	u32 m_m;
	Vecd m_startvec;
	LancznoInit m_init;
};
