module math.matrix;

private import std.traits;
private import std.conv : to;
private import std.exception : assumeUnique;
private import std.typecons;
private import std.math;
private import std.complex;
private import std.functional : binaryFun;
private import std.algorithm : min, max;
private import std.random : uniform;
private import std.format;

private import mir.ndslice;
private import math.linear;
private import math.numerics;
import mir.random;

/**
 * Dense matrix stored in a (row-major) contigous chunk of memory.
 * Copy-Construction is disabled. Use .dup or .move accordingly.
 */
struct Matrix(T)
{
	ContiguousMatrix!T data;

	/** allocate new (uninitialized) matrix */
	this(size_t height, size_t width)
	{
		data = stdcUninitSlice!T(height, width);
	}

	/** copy data from arbitrary matrix-like object */
	this(M)(auto ref M m)
	{
		this(m.length!0, m.length!1);
		for(size_t i = 0; i < height; ++i)
			for(size_t j = 0; j < width; ++j)
				this[i,j] = m[i,j];
	}

	/** create zero-matrix */
	static Matrix zero(size_t n, size_t m)
	{
		auto r = Matrix(n, m);
		r[] = T(0);
		return r;
	}

	/** create (square) identity matrix */
	static Matrix identity(size_t n)
	{
		auto r = Matrix(n, n);
		r[] = T(0);
		r[].diagonal[] = T(1);
		return r;
	}

	/** free memory */
	~this()
	{
		stdcFreeSlice(data);
	}

	/** disable copy-construction. use `.dup` or `.move` instead */
	@disable this(this);

	/** height of the matrix */
	size_t height() const @property
	{
		return data.length!0;
	}

	/** width of the matrix */
	size_t width() const @property
	{
		return data.length!1;
	}

	ContiguousMatrix!T opSlice()
	{
		return data;
	}

	ContiguousMatrix!(const(T)) opSlice() const
	{
		return data.toConst;
	}

	ContiguousMatrix!(immutable(T)) opSlice() immutable
	{
		return data.toImmutable;
	}

	ref T opIndex(size_t i, size_t j)
	{
		return this[][i, j];
	}

	ref const(T) opIndex(size_t i, size_t j) const
	{
		return this[][i, j];
	}

	ref immutable(T) opIndex(size_t i, size_t j) immutable
	{
		return this[][i, j];
	}

	void opIndexAssign(B)(auto ref B b)
	{
		this[][] = b;
	}

	void opIndexAssign(B)(auto ref B b, size_t i, size_t j)
	{
		this[][i,j] = b;
	}

	/** human-readable repesenation */
	string toString() const @property
	{
		string s;
		auto strings = slice!string(height, width);
		auto pitch = new size_t[width];

		for(size_t i = 0; i < height; ++i)
			for(size_t j = 0; j < width; ++j)
			{
				static if(isFloatingPoint!T || is(T : Complex!R, R))
					strings[i,j] = format("%.3g", this[i,j]);
				else
					strings[i,j] = to!string(this[i,j]);
				pitch[j] = max(pitch[j], strings[i,j].length);
			}

		for(size_t i = 0; i < height; ++i)
		{
			if(i == 0)
				s ~= "⎛";
			else if(i == height-1)
				s ~= "⎝";
			else
				s ~= "⎜";


			for(size_t j = 0; j < width; ++j)
			{
				for(int k = 0; k < pitch[j]+1-strings[i,j].length; ++k)
					s ~= " ";
				s ~= strings[i,j];
			}

			if(i == 0)
				s ~= " ⎞\n";
			else if(i == height-1)
				s ~= " ⎠";
			else
				s ~= " ⎟\n";
		}
		return s;
	}

	/** make a new copy of the matrix */
	Matrix dup() const
	{
		auto r = Matrix(height, width);
		r[] = this[];
		return r;
	}

	/** dup with row permutation */
	Matrix dup(const(int)[] p) const
	{
		auto r = Matrix(height, width);
		for(int i = 0; i < height; ++i)
			r[][i, 0..$] = data[][p[i], 0..$];
		return r;
	}

	/** dup with row  and column permutation */
	Matrix dup(const(int)[] p, const(int)[] q) const
	{
		auto r = Matrix(height, width);
		for(size_t i = 0; i < height; ++i)
			for(size_t j = 0; j < width; ++j)
				r[i,j] = this[p[i], q[j]];
		return r;
	}

	/** Matrix * / scalar */
	Matrix opBinary(string op, S)(auto ref S b) const
		if((op == "*" || op == "/") && is(typeof(T(b))))
	{
		auto a = Matrix(height, width);
		a[] = this[].opBinary!op(b);
		return a;
	}

	/** scalar * Matrix */
	Matrix opBinaryRight(string op, S)(auto ref S b) const
		if(op == "*" && is(typeof(T(b))))
	{
		Matrix!T a = Matrix(height, width);
		a[] = this[].opBinaryRight!op(b);
		return a;
	}

	/** Matrix +- Matrix */
	Matrix opBinary(string op)(auto ref const Matrix b) const
		if(op == "+" || op == "-")
	{
		if(width != b.width || height != b.height)
			throw new Exception("matrix dimension mismatch");

		auto a = Matrix(height, width);
		a[] = this[].opBinary!op(b[]);
		return a;
	}

	/** Matrix * Matrix */
	Matrix opBinary(string op)(auto ref const Matrix b) const
		if(op == "*")
	{
		if(width != b.height)
			throw new Exception("matrix dimension mismatch");

		auto a = Matrix(height, b.width);
		for(size_t i = 0; i < height; ++i)
			for(size_t j = 0; j < b.width; ++j)
				a[i,j] = T(0).reduce!"a+b"(this[][i,0..$]*b[][0..$,j]);
		return a;
	}

	auto adjoint() const
	{
		static if(is(T : Complex!R, R))
			return this[].transposed.map!conj;
		else
			return this[].transposed;
	}

	/** L1 norm */
	RealTypeOf!T norm1() const @property
	{
		return RealTypeOf!T(0).reduce!"a+b"(this[].map!abs);
	}

	/** squared L2 norm */
	RealTypeOf!T sqNorm2() const @property
	{
		return RealTypeOf!T(0).reduce!"a+b"(this[].map!sqAbs);
	}

	/** L2 norm */
	RealTypeOf!T norm2() const @property
	{
		return sqNorm2.sqrt;
	}

	/** L_infinity norm */
	RealTypeOf!T normInf() const @property
	{
		return RealTypeOf!T(0).reduce!max(this[].map!abs);
	}

	/** most common choice for numerics */
	alias norm = normInf;

	/** matrix ^^ int */
	Matrix!T opBinary(string op)(long exp) const
		if(op == "^^")
	{
		assert(width > 0);
		if(width != height)
			throw new Exception("can not take a power of a non-square matrix");
		if(exp <= 0)
			throw new Exception("non-positive powers not implemented (yet?)");

		Matrix!T r;
		auto base = this.dup;

		while(exp)
		{
			if(exp & 1)
			{
				if(r.width == 0)
					r = base.dup;
				else
					r = r * base;
			}
			exp >>= 1;
			base = base * base;
		}
		return r;
	}

	/** compute LU decomposition */
	DenseLU!T lu()()
	{
		return DenseLU!T(this);
	}

	/** compute LDL decomposition */
	DenseLDL!T ldl()()
	{
		return DenseLDL!T(this);
	}

	/** compute QR decomposition */
	DenseQR!T qr()()
	{
		return DenseQR!T(this);
	}

	DenseHessenberg!T hessenberg()()
	{
		return DenseHessenberg!T(this);
	}

	DenseSchur!T schur()()
	{
		return DenseSchur!T(this);
	}

	/** generate a matrix with random elements */
	static Matrix!T random(Rng)(ref Rng rng, size_t height, size_t width)
		if (isSaturatedRandomEngine!Rng)
	{
		auto r = Matrix!T(height, width);
		foreach(ref x; r[].flattened)
			static if(isFloatingPoint!T)
				x = rng.rand!T;
			else static if(is(T : Complex!R, R))
				x = T(rng.rand!R, rng.rand!R);
			else
				x = T.random(rng);
		return r;
	}

	/** generate a hermitian matrix with random elements */
	static Matrix!T randomHermitian(Rng)(ref Rng rng, size_t n)
		if (isSaturatedRandomEngine!Rng)
	{
		auto r = Matrix!T(n, n);
		static if(isFloatingPoint!T)
		{
			r[].eachUploPair!((ref u, ref l){ u = rng.rand!T; l = u; }, false);
			r[].diagonal.each!((ref x){ x = rng.rand!T; });
		}
		else static if(is(T : Complex!R, R))
		{
			r[].eachUploPair!((ref u, ref l){ u = T(rng.rand!R, rng.rand!R); l = u.conj; }, false);
			r[].diagonal.each!((ref x){ x = T(rng.rand!R); });
		}
		else
		{
			static assert(false);
		}
		return r;
	}
}

/+
/** efficient storage of a square band matrix */
struct BandMatrix(T)
{
	const int kl, ku; // number of sub/super-diagonals
	ContiguousMatrix!(immutable T) data;
	mixin MatrixHelper!T;

	this(int kl, int ku, ContiguousMatrix!(immutable T) data)
	{
		this.data = data;
		this.kl = kl;
		this.ku = ku;

		if(kl+ku+1 != data.length!0 || kl >= height || ku >= height)
			throw new Exception("invalid dimensions of band matrix");
	}

	size_t height() const @property
	{
		return data.length!1; // width == height
	}

	size_t width() const @property
	{
		return data.length!1;
	}

	ref const(T) opIndex(size_t i, size_t j) const
	{
		if(i+ku < j || i > j+kl)
			return zero;

		return data[i-j+ku,j];
	}

	static const(T) zero;

	static this()
	{
		zero = T(0);
	}
}
+/

/**
 * read-only access to lower/upper triangular matrix with/without implicit 1's on the diagonal
 */
struct TriangularView(T, bool lower, bool implicitOne, int offDiag = 0)
{
	ContiguousMatrix!(const(T)) data;

	size_t length(size_t i)() const @property
	{
		return data.length!i;
	}

	this(ContiguousMatrix!(const(T)) data)
	{
		this.data = data;
	}

	size_t height() const @property
	{
		return data.length!0;
	}

	size_t width() const @property
	{
		return data.length!1;
	}

	const(T) opIndex(size_t i, size_t j) const
	{
		if(implicitOne && i == j)
			return T(1);

		if((lower && i+offDiag < j) || (!lower && i > j+offDiag))
			return T(0);

		return data[i,j];
	}
}

struct DenseLU(T)
{
	Matrix!T m; // L and U matrix stored together
	int[] p, pInv; // row permutation

	this(ref const Matrix!T mat)
	{
		m = mat.dup;
		p = new int[m.height];
		denseComputeLU!T(m[], p);
		pInv = new int[m.height];
		for(int i = 0; i < m.height; ++i)
			pInv[p[i]] = i;
	}

	Matrix!T solve(ref const Matrix!T b)
	{
		auto r = b.dup(p);
		denseSolveLU!T(m[], r[]);
		return r;
	}

	/** lower/left part of the decomposition */
	auto l() const
	{
		return TriangularView!(T, true, true)(m[].toConst);
	}

	/** upper/right part of the decomposition */
	auto u() const
	{
		return TriangularView!(T, false, false)(m[].toConst);
	}

	/**
	 * Reconstruct the original matrix from the decomposition.
	 */
	Matrix!T a()
	{
		return (Matrix!T(l)*Matrix!T(u)).dup(pInv);
	}
}

struct DenseLDL(T)
{
	Matrix!T m; // uses only lower triangle

	this(ref const Matrix!T mat)
	{
		m = mat.dup;
		denseComputeLDL!T(m[]);
	}

	Matrix!T solve(ref const Matrix!T b)
	{
		auto r = b.dup;
		denseSolveLDL!T(m[], r[]);
		return r;
	}

	auto l() const
	{
		return TriangularView!(T, true, true)(m[].toConst);
	}

	Matrix!T a()
	{
		auto b = Matrix!T(Matrix!T(l).adjoint);
		for(int i = 0; i < m.width; ++i)
			b[][i, 0..$] *= m[i,i];
		return Matrix!T(l)*b;
	}
}

struct DenseQR(T)
{
	Matrix!T m;
	RealTypeOf!T[] beta;

	this(ref const Matrix!T mat)
	{
		m = mat.dup;
		beta = new RealTypeOf!T[m.height];
		denseComputeQR!T(m[], beta);
	}

	Matrix!T solve(ref const Matrix!T b)
	{
		auto r = b.dup();
		denseSolveQR!T(m[], beta, r[]);
		return r;
	}

	auto q()
	{
		auto q = Matrix!T.identity(beta.length);
		for(int k = cast(int)beta.length-1; k >= 0; --k)
			applyHouseholder!T(q[][k..$,0..$], m[][k+1..$,k], beta[k]);
		return q;
	}

	auto r()
	{
		return TriangularView!(T, false, false)(m[]);
	}

	Matrix!T a()
	{
		return q*Matrix!T(r);
	}
}

struct DenseHessenberg(T)
{
	Matrix!T m;
	RealTypeOf!T[] beta;

	this(ref const Matrix!T mat)
	{
		m = mat.dup;
		beta = new RealTypeOf!T[m.height];
		denseComputeHessenberg!T(m[], beta);
	}

	auto q()
	{
		auto q = Matrix!T.identity(beta.length);
		for(int k = cast(int)beta.length-2; k >= 0; --k)
			applyHouseholder!T(q[][k+1..$,0..$], m[][k+2..$,k], beta[k]);
		return q;
	}

	auto h()
	{
		return TriangularView!(T, false, false, 1)(m[]);
	}

	Matrix!T a()
	{
		return q*Matrix!T(h)*Matrix!T(q.adjoint);
	}
}

struct DenseSchur(T)
{
	Matrix!T u;
	Matrix!T q;

	this(ref const Matrix!T mat)
	{
		u = mat.dup;
		q = Matrix!T.identity(u.height);
		denseComputeSchur!T(u[], q[]);
	}

	Matrix!T a()
	{
		return q*u*Matrix!T(q.adjoint);
	}
}
