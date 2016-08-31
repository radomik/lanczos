#include "utils.hpp"
#include "string_utils.hpp"
#include "file_utils.hpp"
#include <vector>
#include <fstream>
#include <algorithm>
#include <clapack/clapack.h>

typedef uint8_t  u8;
typedef int32_t  s32;
typedef uint32_t u32;
typedef uint64_t u64;

// must use double because of MRRR step
typedef double flt;

class Vecd : public std::vector<flt> {
public:
	Vecd() : std::vector<flt>() { }
	Vecd(size_t size) : std::vector<flt>(size) { }
	Vecd(size_t size, flt val) : std::vector<flt>(size, val) { }
	Vecd(const Vecd& src) : std::vector<flt>(src) { }
	
	flt* data() { return &front(); }
	const flt* data() const { return &front(); }

	Vecd& random() {
		random(size(), data());
		return *this;
	}
	
	flt norm() {
		return norm(size(), data());
	}
	
	Vecd& zero() {
		zero(size(), data());
		return *this;
	}
	
	flt getNorm() const {
		return getNorm(size(), data());
	}
	
	flt getNormPow() const {
		return getNormPow(size(), data());
	}
	
	void print(bool dbg, const char* name) const {
		if (dbg) print(size(), data(), stderr, name);
	}
	
	void print(FILE* file, const char* name) const {
		print(size(), data(), file, name);
	}
	
	Vecd& mulByScalar(flt s) {
		mulByScalar(size(), data(), data(), s);
		return *this;
	}
	
	Vecd& mulByScalar(const Vecd& v, flt s) {
		assert(size() == v.size());
		mulByScalar(size(), data(), v.data(), s);
		return *this;
	}
	
	Vecd& divByScalar(flt s) {
		divByScalar(size(), data(), data(), s);
		return *this;
	}
	
	Vecd& divByScalar(const Vecd& v, flt s) {
		assert(size() == v.size());
		divByScalar(size(), data(), v.data(), s);
		return *this;
	}
	
	Vecd& mulInvByScalar(flt s) {
		mulInvByScalar(size(), data(), data(), s);
		return *this;
	}
	
	Vecd& mulInvByScalar(const Vecd& v, flt s) {
		assert(size() == v.size());
		mulInvByScalar(size(), data(), v.data(), s);
		return *this;
	}
	
	Vecd& add(const Vecd& b) {
		assert(size() == b.size());
		add(size(), data(), data(), b.data());
		return *this;
	}
	
	Vecd& add(const Vecd& a, const Vecd& b) {
		assert(a.size() == b.size());
		assert(size() == b.size());
		add(size(), data(), a.data(), b.data());
		return *this;
	}
	
	Vecd& sub(const Vecd& b) {
		assert(size() == b.size());
		sub(size(), data(), data(), b.data());
		return *this;
	}
	
	Vecd& sub(const Vecd& a, const Vecd& b) {
		assert(a.size() == b.size());
		assert(size() == b.size());
		sub(size(), data(), a.data(), b.data());
		return *this;
	}
	
	flt dotProduct(const Vecd& b) const {
		assert(size() == b.size());
		return dotProduct(size(), data(), b.data());
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

class Veci : public std::vector<s32> {
public:
	Veci() : std::vector<s32>() { }
	Veci(size_t size) : std::vector<s32>(size) { }
	Veci(size_t size, s32 val) : std::vector<s32>(size, val) { }
	Veci(const Veci& src) : std::vector<s32>(src) { }
	
	s32* data() { return &front(); }
	const s32* data() const { return &front(); }
	
	Veci& zero() {
		zero(size(), data());
		return *this;
	}
	
	void print(FILE* file, const char* name) const {
		print(size(), data(), file, name);
	}
	
	static void zero(size_t size, s32* data) {
		memset(data, 0, size*sizeof(data[0]));
	}
	
	static void print(size_t size, const s32* data, FILE* file, const char* name) {
		fprintf(file, "%s = [\n", name);
		for (size_t i = 0; i < size; i++) {
			fprintf(file, "%s%d ", i ? " ;\n" : "", data[i]);
		}
		fprintf(file, " ]\n\n");
	}
};

class Vecb : public std::vector<u8> {
public:
	Vecb() : std::vector<u8>() { }
	Vecb(size_t size) : std::vector<u8>(size) { }
	Vecb(size_t size, u8 val) : std::vector<u8>(size, val) { }
	Vecb(const Vecb& src) : std::vector<u8>(src) { }
	
	u8* data() { return &front(); }
	const u8* data() const { return &front(); }
	
	Vecb& zero() {
		zero(size(), data());
		return *this;
	}
	
	void print(FILE* file, const char* name) const {
		print(size(), data(), file, name);
	}
	
	static void zero(size_t size, u8* data) {
		memset(data, 0, size*sizeof(data[0]));
	}
	
	static void print(size_t size, const u8* data, FILE* file, const char* name) {
		fprintf(file, "%s = [\n", name);
		for (size_t i = 0; i < size; i++) {
			fprintf(file, "%s%u ", i ? " ;\n" : "", data[i]);
		}
		fprintf(file, " ]\n\n");
	}
};

class Matrixd : public Vecd {
public:
	Matrixd() : Vecd() { }
	Matrixd(size_t size) : Vecd(size) { }
	Matrixd(size_t size, flt val) : Vecd(size, val) { }
	Matrixd(const Matrixd& src) : Vecd(src) { }
};

class TrMatrixd : public Vecd {
public:
	TrMatrixd() : Vecd() { }
	TrMatrixd(size_t size) : Vecd(size) { }
	TrMatrixd(size_t size, flt val) : Vecd(size, val) { }
	TrMatrixd(const TrMatrixd& src) : Vecd(src) { }
	
	//NOTE: size() returns nelem
	//NOTE: dim() returns n
	
	size_t dim() const { return m_dim; }
	
	int readCsv(const char* filename, const char* separator) {
		std::ifstream file(filename);

		if (!file) {
			fprintf(stderr, "%s: Failed to open '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
			return -1;
		}

		const std::string delim = separator;
		std::string line;
		size_t  i = 0;
		size_t  di = 0;
		flt*    dt = NULL;
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
				dt = data();
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

					if (sscanf(token.c_str(), "%lf", &dt[di++]) != 1) {
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

				if (sscanf(line.c_str(), "%lf", &dt[di++]) != 1) {
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
			fprintf(stderr, "%s: CSV parsing of '%s' error occured\n", __FUNCTION__, filename);
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
		
		const size_t nelem = sizeFromDim(nu);
		
		fprintf(stderr, "%s: Reading matrix of dim=%u, nelem=%zu, sizeof(A[0])=%zu\n", __FUNCTION__, nu, nelem, sizeof(data()[0]));
		
		resize(nelem);
		fprintf(stderr, "%s: Size of data is %zu\n", __FUNCTION__, size());
		
		flt* dt = data();
		if ((ret=fread(dt, sizeof(dt[0]), nelem, bin)) != nelem) {
			fprintf(stderr, "%s: Error reading matrix content [nelem=%zu, ret=%zu, err=%s]\n", __FUNCTION__, nelem, ret, strerror(errno));
			fclose(bin);
			return -1;
		}
		
		fclose(bin);
		return 0;
	}
	
	int readCsvOrBin(const std::string& work_dir, const std::string& filename_without_ext) {
		const CsvBinPaths p(work_dir, filename_without_ext);
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
		print(data(), dim(), file, name);
	}
	
	int saveCsv(const char* filename) const {
		return saveCsv(data(), dim(), filename);
	}
	
	int saveBin(const char* filename) const {
		return saveBin(data(), dim(), filename);
	}
	
	bool equals(const TrMatrixd& other, flt eps, size_t* first_diff_index, flt* first_diff) const {
		return equals(data(), dim(), other.data(), other.dim(), eps, first_diff_index, first_diff);
	}
	
	static size_t sizeFromDim(size_t datadim) {
		return (datadim*(datadim+1))>>1;
	}
	
	static void mulByVector(Vecd& res, const TrMatrixd& A, const Vecd& b) {
		assert(A.size() == b.size());
		assert(res.size() == b.size());
		mulByVector(res.data(), A.data(), b.data(), b.size());
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

		const uint32_t nu = (uint32_t)dim;
		const size_t nelem = sizeFromDim(dim);
		
		assert(nu == dim);

		fprintf(stderr, "%s: Saving matrix of dim=%u, nelem=%zu, sizeof(A[0])=%zu\n", __FUNCTION__, nu, nelem, sizeof(A[0]));

		FILE* bin = fopen(filename, "wb");
		if (! bin) {
			fprintf(stderr, "%s: Error opening output file '%s' [%s]\n", __FUNCTION__, filename, strerror(errno));
			return -1;
		}
		
		int ret = (fwrite(&nu, sizeof(nu), 1, bin) != 1) || (fwrite(A, sizeof(A[0]), nelem, bin) != nelem);
		
		if (ret) {
			fprintf(stderr, "%s: Error writing file '%s' [ret=%d, err=%s]\n", __FUNCTION__, filename, ret, strerror(errno));
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
		
		const size_t nelem = sizeFromDim(Xn);
		
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
	LancznoInit(size_t m) : m_a(m, 0.0), m_b(m, 0.0) { }
	
	/**
	 * @param A{n.n}      linearized reduced symmetric square matrix used to compute eigen values/vectors for
	 * @param stv{n}      starting vector 
	 * @param dbg         whether print debug information
	 * 
	 * Produces results in m_a, m_b, m_anorm
	 * 
	 */
	void run(const TrMatrixd& A, const Vecd& stv, bool dbg) {
		const size_t n = stv.size();
		const size_t m = m_a.size();
		Vecd v(stv);                             // v{n} = stv{n}
		Vecd v2(n);                              // v2{n}
		Vecd r(n);                               // r{n}

		v.print(dbg, "v");
		m_a.print(dbg, "a");
		m_b.print(dbg, "b");

		// initial Lanczos iteration, m steps
		for (size_t k = 0; ; ) {
			if (dbg) fprintf(stderr, "for: k == %zu\n", k);
			
			if (k == 0) {
				TrMatrixd::mulByVector(r, A, v); // r{n}  = A{n}{n} * v{n}
			} else {
				TrMatrixd::mulByVector(r, A, v); // r{n}  = A{n}{n} * v{n}
				v2.mulByScalar(m_b[k-1]);        // v2{n} = b[k-1]{1} * v2{n}
				r.sub(v2);                       // r{n}  = r{n} - v2{n}
			}

			r.print(dbg, "r");
			
			m_a[k] = v.dotProduct(r);            // a[k]{1} = SUM v[i]*r[i]
			m_a.print(dbg, "a");

			v2.mulByScalar(v, m_a[k]);           // v2{n}   = a[k]{1} * v{n}
			r.sub(v2).print(dbg, "r");           // r{n}    = r{n} - v2{n}
			
			m_b[k] = r.getNorm();                // b[k]{1} = (SUM r[i]^2)^0.5
			
			m_b.print(dbg, "b");

			// estimate |A|_2 by |T|_1
			if (k == 0) {
				m_anorm = fabs(m_a[0] + m_b[0]); // anorm = |a[0] + b[0]|
			} else {
				m_anorm = max(m_anorm, m_b[k-1]+fabs(m_a[k])+m_b[k]);
			}                                    // anorm = max(anorm, b[k-1]+|a[k]|+b[k])

			if (dbg) fprintf(stderr, "anorm = %f\n\n", m_anorm);

			if (++k == m) {
				break;
			}

			// prepare next step, k = {1 ... m-1}
			v2 = v;                              // v2{n} = v{n}
			v2.print(dbg, "v2");
			
			v.divByScalar(r, m_b[k-1]).print(dbg, "v"); // EACH v[i] = r[i]/b[k-1]
		}
	}
	
	Vecd m_a;
	Vecd m_b;
	flt  m_anorm;
};

class Mrrr {
public:
	Mrrr(size_t m) : m_S(m), m_ritz(m) {
		 m_eps = 1e-13;
	}
	
	int run(const Vecd& a, const Vecd& b, bool dbg) {
		const size_t m = m_ritz.size();
		assert(a.size() == m);
		assert(b.size() == m);
		assert(m_S.size() == m*m);
		
		// matrix size
		int im = m; // m
		assert((size_t)im == m);

		char jobz[] = "V";
		char range[] = "A";

		// main diagonal
		Vecd d(a);

		// sub-/superdiagonal
		Vecd e(b);

		double vl = 0.0, vu = 0.0;  // unused
		int il = 0, iu = 0;         // unused
		double abstol = m_eps;      // unused?
		int numEval;                // number of eigenvalues found

		// misc arguments
		double v_work1 = 0.0;
		int v_iwork1 = 0;
		int ldz = m;
		int lwork = -1;
		int liwork = -1;
		int info;

		Veci isuppz(m<<1, 0);
		double* work1 = &v_work1;
		int* iwork1 = &v_iwork1;

		// query for workspace size
		dstegr_( jobz, range, &im, d.data(), e.data(), &vl, &vu, &il, &iu, &abstol,
			&numEval, m_ritz.data(), m_S.data(), &ldz, isuppz.data(), work1, &lwork, iwork1, &liwork,
			&info );

		if (info < 0) {
			fprintf(stderr, "%s: [ERROR 1] Input to dstegr had an illegal value [info: %d]\n", __FUNCTION__, info);
			return info;
		}

		if (dbg) {
			fprintf(stderr, "%s: [1] lwork: %d (%f), liwork: %d (%d)\n", __FUNCTION__, lwork, work1[0], liwork, iwork1[0]);
		}

		lwork  = (int)work1[0];
		liwork = iwork1[0];

		if (dbg) {
			fprintf(stderr, "%s: [2] lwork: %d (%f), liwork: %d (%d)\n", __FUNCTION__, lwork, work1[0], liwork, iwork1[0]);
		}
		
		assert(lwork > 0);
		assert(liwork > 0);
		
		Vecd work(lwork, 0.0);
		Veci iwork(liwork, 0);

		// call LAPACK routine
		dstegr_( jobz, range, &im, d.data(), e.data(), &vl, &vu, &il, &iu, &abstol,
			&numEval, m_ritz.data(), m_S.data(), &ldz, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
			&info );

		if (info != 0) {
			fprintf(stderr, "%s: [ERROR 2] Input to dstegr had an illegal value [info: %d]\n", __FUNCTION__, info);
		}
		return info;
	}
	
	Matrixd m_S;
	Vecd    m_ritz;
	flt     m_eps;
};

class Rescon {
public:
	Rescon(size_t m) : m_idx(m) {
		 m_eps = 1e-13;
	}
	
	/**
	 * Residual estimation and removing non-converged and spurious Ritz values.
	 * 
	 * @param S     Input non-symmetric matrix of size m*m obtained in MRRR step
	 * @param ritz  Input ritz vector of size m obtained in MRRR step
	 * @param b     Input b vector of size m obtained in LANCNO_INIT step
	 * @param anorm Input value of anorm obtained in LANCNO_INIT step
	 * @param dbg   Print debug info
	 * 
	 * @return 0 on success, 1 when all RITZ values has been removed
	 * 
	 */
	int run(
		const Matrixd  S,
		const Vecd&    ritz,
		const Vecd&    b,
		const double   anorm,
		bool           dbg
	)
	{
		const size_t m        = b.size();
		size_t mm             = 0;
		const flt* S_last_row = S.data() + m*(m-1);
		const flt  b_last     = b[m-1];
		const flt  tol        = m_eps * anorm;
		Vecd   temp_lres(m);
		Vecd   temp_cul(m);
		
		if (dbg) {
			fprintf(stderr, "%s: eps=%g, anorm=%g, tol=%g\n", __FUNCTION__, m_eps, anorm, tol);
		}
		
		for (size_t i = 0; i < m; i++) {
			const flt lres_i = fabs(b_last * S_last_row[i]); // lres = abs(b(m)*S(m,:))'; % residual estimation
			const flt cul_i  = fabs(S[i]);                    // cul = abs(S(1,:))';
			
			// remove non-converged and spurious Ritz values
			if ((m_idx[i] = (lres_i < tol) && (cul_i > tol))) {
				mm++;
			}
			
			temp_lres[i] = lres_i;
			temp_cul[i] = cul_i;
		}
		
		if (dbg) {
			fprintf(stderr, "%s: m=%zu, mm=%zu\n", __FUNCTION__, m, mm);
		}
		
		if (mm == 0) {
			fprintf(stderr, "%s: [ERROR] All Ritz values are marked to be removed\n", __FUNCTION__);
			return -1;
		}
		
		m_lres.resize(mm);
		m_cul.resize(mm);
		m_ritz.resize(mm);
		m_S.resize(m * mm);
		m_mm = mm;
		for (size_t i = 0, k = 0; i < m; i++) {
			if (m_idx[i]) {
				m_lres[k] = temp_lres[i];
				m_cul[k]  = temp_cul[i];
				m_ritz[k] = ritz[i];
				k++;
			}
		}
		
		// rewrite only columns s[i][k] which has idx[k] = 1
		flt*       dst_data = m_S.data();
		const flt* src_data = S.data();
		for (size_t i = 0; i < m; i++) {
			const size_t sd = i * mm;
			const size_t ss = i * m;
			for (size_t j = 0, k = 0; j < m; j++) {
				if (m_idx[j]) {
					dst_data[sd + k] = src_data[ss + j];
					k++;
				}
			}
		}
		
		return 0;
	}
	
	Vecd    m_lres; ///< output vector lres of size mm, mm <= m
	Vecd    m_cul;  ///< output vector cul of size mm, mm <= m
	Vecd    m_ritz; ///< output vector ritz of size mm, mm <= m
	Vecb    m_idx;  ///< output 0/1 vector idx of size m
	Matrixd m_S;    ///< output s matrix of size m*mm, mm <= m
	flt     m_mm;   ///< output mm value
	flt     m_eps;  ///< parameter, double epsilon value
};

class Dc {
public:
	Dc() {
		m_tol = 1e-8;
	}

	int run(const Vecd& ritz, size_t mm, bool dbg) {
		if (mm > UINT_MAX) {
			fprintf(stderr, "%s: mm=%zu shall not be greater than UINT_MAX=%u\n", __FUNCTION__, mm, UINT_MAX);
			return -1;
		}
		m_ci.clear();
		u32 i = 0;
		for (u32 k = 1; k < mm; k++) {
			const flt d = fabs(ritz[k] - ritz[i]);
			if (d > m_tol) {
				m_ci.push_back(std::make_pair(i, k-1));
				i = k++;
			}
		}
		m_ci.push_back(std::make_pair(i, mm));
		return 0;
	}
	
	std::vector< std::pair<u32, u32> > m_ci;
	flt m_tol;
};

class Householder {
public:
	Householder(size_t m, size_t cis) : m_S2(m*cis), m_e(cis), m_c(cis), m_m(m) { }
	
	void run(
		const std::vector< std::pair<u32, u32> >& ci,
		const Matrixd&  S,
		const Vecd&     ritz,
		const Vecd&     lres,
		const Vecd&     cul,
		Vecb&           idx
	) {
		Vecd x;
		Vecd u;
		const size_t m = m_m;
		const size_t cis = ci.size();
		for (size_t i = 0; i < cis; i++) {
			const std::pair<u32,u32> cii = ci[i];
			u32 max_cul_idx = 0;                // j
			flt max_cul_value = cul[cii.first]; // y
			for (u32 j = cii.first+1, ix = 1; j <= cii.second; j++, ix++) {
				if (cul[j] > max_cul_value) {
					max_cul_value = cul[j];
					max_cul_idx = ix;
				}
			}
			const size_t ji = cii.first + max_cul_idx;
			
			// x and u vectors are of size no more than mm
			// x = {S[0][cii.first], S[0][cii.first+1], ..., S[0][cii.second]}
			x.resize(cii.second - cii.first + 1);
			memcpy(x.data(), S.data()+cii.first, x.size());
			u = x;
			
			const flt sd = x[max_cul_idx];
			if (sd > 0.0) {
				u[max_cul_idx] += x.norm();
			} else {
				if (sd < 0.0) {
					u[max_cul_idx] -= x.norm();
				}
			}
			
			// s = S(:,ji)-2/(u'*u)*(S(:,cii)*u)*u(j);
			// S2(:,i) = s;
			const flt un = u.getNormPow();
			const flt uj = u[max_cul_idx];
			
			for (u32 p = 0; p < m; p++) { //TODO: paralellize this loop with ISPC/OpenMP
				flt st = 0.0;
				const u64 pm = p*m;
				for (u32 q = cii.first, r = 0; q <= cii.second; q++, r++) {
					st += S[pm + q] * u[r];
				}
				m_S2[p*cis+i] = S[pm + ji] - 2.0 / un*st*uj;
			}
			
			for (u32 q = cii.first; q <= cii.second; q++) {
				idx[q] = false; // idx is a 0/1 vector obtained from Rescon step
			}
			idx[ji] = true;
			m_e[i] = ritz[ji];
			m_c[i] = lres[ji];
		}
	}
	
	Matrixd       m_S2;
	Vecd          m_e;
	Vecd          m_c;
	const size_t  m_m;
};

class EigenVec {
public:
	EigenVec(size_t n, size_t m, size_t cis) : m_X(n*cis), m_m(m), m_cis(cis) { }
	
	void run(
		const TrMatrixd&  A, 
		const Vecd&       stv, 
		const Matrixd&    S2,
		bool              dbg) 
	{
		const size_t n   = stv.size();
		const size_t m   = m_m;
		const size_t cis = m_cis;
		Vecd v(stv);                             // v{n} = stv{n}  
		Vecd v2(n), r(n);
		flt* X = m_X.data();
		flt  a, b, b2 = 0.0;

		for (size_t k = 0; ; k++) {
			const flt* S2k = S2.data() + k * cis;
			for (size_t i = 0; i < n; i++) {
				flt* Xi = X + i * cis;
				for (size_t j = 0; j < cis; j++) {
					Xi[j] += v[i] * S2k[j];
				}
			}
			
			if (k == m-1) {
				break;
			}
			
			if (k == 0) {
				TrMatrixd::mulByVector(r, A, v); // r{n} = A{n}{n}  * v{n}
			} else {
				TrMatrixd::mulByVector(r, A, v); // r{n}  = A{n}{n} * v{n}
				v2.mulByScalar(b2);              // v2{n} = b2{1}   * v2{n}
				r.sub(v2);                       // r{n}  = r{n}    - v2{n}
			}
			
			a = v.dotProduct(r);                 // a{1}  = SUM v[i]*r[i]
			v2.mulByScalar(v, a);                // v2{n} = a{1} * v{n}
			r.sub(v2);                           // r{n}  = r{n} - v2{n}
			b = r.getNorm();                     // b{1}  = (SUM r[i]^2)^0.5
			v2 = v;                              // v2{n} = v{n}
			b2 = b;                              // b2{1} = b{1}
			v.divByScalar(r, b);                 // EACH v[i] = r[i]/b
		}

		Vecd xv(cis, 0.0);

		// xv is a vector of sums of squared elements in each column of X
		for (size_t i = 0; i < n; i++) {
			const flt* Xi = &X[i*cis];
			for (size_t j = 0; j < cis; j++) {
				xv[j] += Xi[j]*Xi[j];
			}
		}

		for (size_t i = 0; i < cis; i++) {
			// xv[i] = 1.0 / diag(X)[i] / sqrt(xv[i])
			xv[i] = 1.0 / X[i*cis + i] / sqrt(xv[i]);
		}

		//TODO: Determine results of this and above operations
		// mainly (see lancno.m):
		//		X = X + v*S2(m,:);
		//		X = X*diag(1./normc(X));
		//Matrixd::mulByVector(X, X, xv);
	}
	
	Matrixd       m_X;
	const size_t  m_m;
	const size_t  m_cis;
};

class LancznoSuite {
public:
	/**
	 * @param n input matrix size
	 * @param m Lanczos iteration count
	*/
	LancznoSuite(u32 n, u32 m) : m_n(n), m_m(m) { 
		m_eps_mrrr   = 1e-13;
		m_eps_rescon = 1e-13;
		m_tol_dc     = 1e-8;
	}
	
	/**
	 * @param A{n.n}      linearized reduced symmetric square matrix used to compute eigen values/vectors for
	 * @param startvec{n} optional starting vector 
	 * @param dbg         whether print debug information
	 */
	int run(const TrMatrixd& A, const Vecd* startvec, bool dbg) {
		size_t m = m_m;
		size_t n = m_n;
		Vecd stv(n);
		
		// 1. Initialize starting vector stv{n}
		if (startvec) {
			stv = *startvec;
		} else {
			stv.random().norm();
		}
		
		// 2. Initial Lanczos iteration
		// input:  A{T,n.n}, stv{n}, m, n
		// output: a{m}, b{m}, anorm
		LancznoInit init(m);
		init.run(A, stv, dbg);
		
		// 3. MRRR step
		// input: a{m}, b{m}, m, eps_mrrr
		// output: S{X,m.m}, ritz{m}
		Mrrr mrrr(m);
		mrrr.m_eps = m_eps_mrrr;
		if (mrrr.run(init.m_a, init.m_b, dbg) != 0) {
			fprintf(stderr, "%s: MRRR step failed\n", __FUNCTION__);
			return -1;
		}
		
		// 4. Rescon step
		// input: S{X,m.m}, ritz{m}, b{m}, anorm, m, eps_rescon
		// output: S'{X,m.mm}, idx{m}, ritz'{mm}, lres{mm}, cul{mm}, mm
		Rescon rescon(m);
		rescon.m_eps = m_eps_rescon;
		if (rescon.run(mrrr.m_S, mrrr.m_ritz, init.m_b, init.m_anorm, dbg) != 0) {
			fprintf(stderr, "%s: Rescon step failed\n", __FUNCTION__);
			return -1;
		}
		
		// 5. DC step
		// input: ritz'{mm}, mm, tol_dc
		// output: ci (vector<pair<u32,u32>> of size <= mm)
		Dc dc;
		dc.m_tol = m_tol_dc;
		if (dc.run(rescon.m_ritz, rescon.m_mm, dbg) != 0) {
			fprintf(stderr, "%s: Rescon step failed\n", __FUNCTION__);
			return -1;
		}
		
		// 6. Householder step
		// input: ci, m, S'{X,m.mm}, ritz'{mm}, lres{mm}, cul{mm}, idx{m}
		// output: idx'{m}, S2{X,m.ci.size()}, m_e{ci.size()}, m_c{cis.size()}
		Householder h(m, dc.m_ci.size());
		h.run(dc.m_ci, rescon.m_S, rescon.m_ritz, rescon.m_lres, rescon.m_cul, rescon.m_idx);
	}
	
	size_t m_n;
	size_t m_m;
	flt    m_eps_mrrr;
	flt    m_eps_rescon;
	flt    m_tol_dc;
};
