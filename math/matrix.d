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

private import jive.array;
private import math.linear;
private import math.numerics;
private import math.permutation : Permutation;

template MatrixHelper(T)
{
	string toString() const @property
	{
		string s;
		auto strings = Array2!string(height, width);
		auto pitch = Array!size_t(width, 0);

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

	Slice2!T dup() const @property
	{
		auto r = Slice2!T(height, width);
		foreach(i,j, ref x; r)
			x = this[i,j];
		return r;
	}

	/** dup with row permutation */
	Slice2!T dup(Permutation p) const @property
	{
		auto r = Slice2!T(height, width);
		foreach(i,j, ref x; r)
			x = this[p(cast(int)i),j];
		return r;
	}

	/** dup with row  and column permutation */
	Slice2!T dup(Permutation p, Permutation q) const @property
	{
		auto r = Slice2!T(height, width);
		foreach(i,j, ref x; r)
			x = this[p(cast(int)i),q(cast(int)j)];
		return r;
	}

    // TODO: efficient + - * for matrices of special structure

	Matrix!T opBinary(string op, M)(M b) const
		if(op == "+" || op == "-")
	{
		if(width != b.width || height != b.height)
			throw new Exception("matrix dimension mismatch");

		auto a = Slice2!T(height, b.width);
		for(size_t j = 0; j < b.width; ++j)
			for(size_t i = 0; i < height; ++i)
				a[i,j] = mixin("this[i,j] "~op~" b[i,j]");
		return Matrix(a.assumeUnique);
	}

	Matrix!T opBinary(string op, M)(M b) const
		if(op == "*")
	{
		if(width != b.height)
			throw new Exception("matrix dimension mismatch");

		auto a = Slice2!T(height, b.width);
		for(size_t i = 0; i < height; ++i)
			for(size_t j = 0; j < b.width; ++j)
			{
				a[i,j] = this[i,0]*b[0,j];
				for(size_t k = 1; k < width; ++k)
					a[i,j] = a[i,j] + this[i,k]*b[k,j];
			}
		return Matrix!T(a.assumeUnique);
	}

	Matrix!T permute(Permutation p) const
	{
		return Matrix!T(this.dup(p).assumeUnique);
	}

	Matrix!T permute(Permutation p, Permutation q) const
	{
		return Matrix!T(this.dup(p, q).assumeUnique);
	}
}

/**
 * dense storage of a (possibly non-square) matrix.
 * column major order in a contigous array.
 */
struct Matrix(T)
{
	Slice2!(immutable(T)) data;
	mixin MatrixHelper!T;

	this(Slice2!(immutable(T)) data)
	{
		this.data = data;
	}

	size_t height() const @property
	{
		return data.size[0];
	}

	size_t width() const @property
	{
		return data.size[1];
	}

	ref const(T) opIndex(size_t i, size_t j) const
	{
		return data[i, j];
	}

	/** transposed matrix (without duplication) */
	Matrix transpose() const
	{
		return Matrix(data.transpose);
	}

	Matrix adjoint() const
	{
		auto r = this.dup.transpose;
		static if(is(T : Complex!R, R))
			foreach(i, j, ref x; r)
				x = conj(x);
		return Matrix(r.assumeUnique);
	}

	/** return squared L2 norm */
	RealTypeOf!T sqNorm() const @property
	{
		RealTypeOf!T sum = 0;
		for(size_t j = 0; j < width; ++j)
			for(size_t i = 0; i < height; ++i)
				sum += this[i,j].sqAbs;
		return sum;
	}

	/** return L2 norm */
	RealTypeOf!T norm() const @property
	{
		return sqNorm.sqrt;
	}

	Matrix!T pow(long exp)
	{
		assert(width > 0);
		if(width != height)
			throw new Exception("can not take a power of a non-square matrix");
		if(exp <= 0)
			throw new Exception("non-positive powers not implemented (yet?)");

		Matrix!T r;
		Matrix!T base = this;

		while(exp)
		{
			if(exp & 1)
			{
				if(r.width == 0)
					r = base;
				else
					r = r * base;
			}
			exp >>= 1;
			base = base * base;
		}
		return r;
	}

	/** compute LU decomposition */
	DenseLU!(T) lu()()
	{
		return DenseLU!(T)(this);
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

	static Matrix!T random(size_t height, size_t width)
	{
		static if(isFloatingPoint!T)
			return build!((i,j)=> cast(T)uniform(-1.0, 1.0))(height, width);
		else static if(is(T : Complex!R, R))
			return build!((i,j)=> T(uniform(-1.0,1.0), uniform(-1.0,1.0)))(height, width);
		else
			return build!((i,j)=> T.random())(height, width);
	}

	static Matrix!T randomHermitian(size_t n)
	{
		static if(isFloatingPoint!T)
			return buildHermitian!((i,j)=> cast(T)uniform(-1.0, 1.0))(n);
		else static if(is(T : Complex!R, R))
			return buildHermitian!((i,j)=> T(uniform(-1.0,1.0), uniform(-1.0,1.0)))(n);
		else
			return buildHermitian!((i,j)=> T.random())(n);
	}

	static Matrix!T random(Ring)(size_t height, size_t width, Ring ring)
	{
		auto data = Slice2!T(height, width);
		for(size_t j = 0; j < width; ++j)
			for(size_t i = 0; i < height; ++i)
				data[i,j] = ring.random();
		return Matrix!T(data.assumeUnique);
	}

	static Matrix!T build(alias fun)(size_t h, size_t w)
	{
		auto data = Slice2!T(h, w);
		for(size_t j = 0; j < w; ++j)
			for(size_t i = 0; i < h; ++i)
				data[i,j] = binaryFun!(fun,"i","j")(i, j);
		return Matrix!T(data.assumeUnique);
	}

	/** only explicitly generates upper/right half */
	static Matrix!T buildHermitian(alias fun)(size_t n)
	{
		auto data = Slice2!T(n, n);
		for(size_t j = 0; j < n; ++j)
			for(size_t i = 0; i <= j; ++i)
			{
				data[i,j] = binaryFun!(fun,"i","j")(i, j);

				static if(is(T : Complex!R, R))
				{
					if(i == j)
						data[i,j] = 0.5*(data[i,j] + conj(data[i,j]));
					else
						data[j,i] = conj(data[i,j]);
				}
				else
				{
					if(i != j)
						data[j,i] = data[i,j];
				}
			}
		return Matrix!T(data.assumeUnique);
	}

	static BandMatrix!T buildBand(alias fun)(size_t n, int kl, int ku)
	{
		assert(kl >= 0 && ku >= 0);
		auto data = Slice2!T(kl+ku+1, n);
		for(size_t j = 0; j < n; ++j)
			for(size_t i = max(ku,j)-ku; i <= min(n-1,j+kl); ++i)
				data[i+ku-j,j] = binaryFun!(fun,"i","j")(i, j);
		return BandMatrix!T(kl, ku, data.assumeUnique);
	}

	/** only explicitly genrates upper/right half */
	static BandMatrix!T buildSymmetricBand(alias fun)(size_t n, int kd)
	{
		assert(kd >= 0);
		auto data = Slice2!T(2*kd+1, n);
		for(size_t j = 0; j < n; ++j)
			for(size_t i = max(kd,j)-kd; i <= j; ++i)
			{
				data[i+kd-j,j] = binaryFun!(fun,"i","j")(i, j);
				if(i != j)
					data[j+kd-i, i] = data[i+kd-j, j];
			}
		return BandMatrix!T(kd, kd, data.assumeUnique);
	}

	static BandMatrix!T buildIdentity(size_t n, T x = T(1))
	{
		auto data = Slice2!T(1, n);
		for(size_t j = 0; j < n; ++j)
			data[j,j] = x;
		return BandMatrix!T(0, 0, data.assumeUnique);
	}

	/**
	 * positive definite matrix (and thus invertible)
	 * but with very bad condition number. Used for numerics testing.
	 */
	static Matrix!T buildHilbert(size_t n)
	{
		return build!((i,j)=>T(1)/(1+cast(int)i+cast(int)j))(n, n);
	}
}

/** efficient storage of a square band matrix */
struct BandMatrix(T)
{
	const int kl, ku; // number of sub/super-diagonals
	Slice2!(immutable(T)) data;
	mixin MatrixHelper!T;

	this(int kl, int ku, Slice2!(immutable(T)) data)
	{
		this.data = data;
		this.kl = kl;
		this.ku = ku;

		if(kl+ku+1 != data.size[0] || kl >= height || ku >= height)
			throw new Exception("invalid dimensions of band matrix");
	}

	size_t height() const @property
	{
		return data.size[1]; // width == height
	}

	size_t width() const @property
	{
		return data.size[1];
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

/**
 * lower/upper triangular matrix with/without implicit 1's on the diagonal
 */
struct TriangularMatrix(T, bool lower, bool implicitOne, int offDiag = 0)
{
	Slice2!(immutable(T)) data;
	mixin MatrixHelper!T;

	this(Slice2!(immutable(T)) data)
	{
		this.data = data;
	}

	size_t height() const @property
	{
		return data.size[0];
	}

	size_t width() const @property
	{
		return data.size[1];
	}

	ref const(T) opIndex(size_t i, size_t j) const
	{
		if(implicitOne && i == j)
			return one;

		if((lower && i+offDiag < j) || (!lower && i > j+offDiag))
			return zero;

		return data[i,j];
	}

	static const(T) zero, one;

	static this()
	{
		zero = T(0);
		one = T(1);
	}
}

struct DenseLU(T)
{
	Slice2!(immutable(T)) m;
	Permutation p;	// row permutation

	this(Matrix!T _m)
	{
		auto p = new int[_m.height];
		auto m = _m.dup;
		denseComputeLU!T(m, p);

		this.m = m.assumeUnique;
		this.p = Permutation(p.assumeUnique);
	}

	Matrix!T solve(Matrix!T b)
	{
		auto rhs = b.dup(p);
		denseSolveLU!T(m, rhs);
		return Matrix!T(rhs.assumeUnique);
	}

	/** lower/left part of the decomposition */
	auto l()
	{
		return TriangularMatrix!(T, true, true)(m);
	}

	/** upper/right part of the decomposition */
	auto u()
	{
		return TriangularMatrix!(T, false, false)(m);
	}

	auto a()
	{
		return (l*u).permute(p.inverse);
	}
}

struct DenseQR(T)
{
	Slice2!(immutable(T)) m;
	immutable(RealTypeOf!T)[] beta;

	this(Matrix!T _m)
	{
		auto beta = new RealTypeOf!T[_m.height];
		auto m = _m.dup;

		denseComputeQR!T(m, beta);

		this.beta = beta.assumeUnique;
		this.m = m.assumeUnique;
	}

	Matrix!T solve(Matrix!T b)
	{
		auto rhs = b.dup();
		denseSolveQR!T(m, beta, rhs);
		return Matrix!T(rhs.assumeUnique);
	}

	auto q()
	{
		auto q = Slice2!T(beta.length, beta.length);
		foreach(i, j, ref x; q)
			x = i==j ? 1 : 0;

		for(int k = cast(int)beta.length-1; k >= 0; --k)
			applyHouseholder!T(q[k..$,0..$], m[k+1..$,k], beta[k]);

		return Matrix!T(q.assumeUnique);
	}

	auto r()
	{
		return TriangularMatrix!(T, false, false)(m);
	}

	auto a()
	{
		return q*r;
	}
}

struct DenseHessenberg(T)
{
	Slice2!(immutable(T)) m;
	immutable(RealTypeOf!T)[] beta;

	this(Matrix!T _m)
	{
		auto beta = new RealTypeOf!T[_m.height];
		auto m = _m.dup;

		denseComputeHessenberg!T(m, beta);

		this.beta = beta.assumeUnique;
		this.m = m.assumeUnique;
	}

	auto q()
	{
		auto q = Slice2!T(beta.length, beta.length);
		foreach(i, j, ref x; q)
			x = i==j ? 1 : 0;

		for(int k = cast(int)beta.length-2; k >= 0; --k)
			applyHouseholder!T(q[k+1..$,0..$], m[k+2..$,k], beta[k]);

		return Matrix!T(q.assumeUnique);
	}

	auto h()
	{
		return TriangularMatrix!(T, false, false, 1)(m);
	}

	auto a()
	{
		return q*h*q.adjoint;
	}
}

struct DenseSchur(T)
{
	Matrix!T u;
	Matrix!T q;

	this(Matrix!T _m)
	{
		auto m = _m.dup;
		auto q = Slice2!T(m.size[0], m.size[1]);

		denseComputeSchur!T(m, q);

		this.u = Matrix!T(m.assumeUnique);
		this.q = Matrix!T(q.assumeUnique);
	}

	auto a()
	{
		return q*u*q.adjoint;
	}
}
