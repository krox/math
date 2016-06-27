module math.matrix;

private import std.traits;
private import std.conv : to;
private import std.exception : assumeUnique;
private import std.typecons;
private import std.complex;
private import std.functional : binaryFun;
private import std.algorithm : min, max;
private import std.random : uniform;

private import jive.array;
private import math.linear;
private import math.complex;

class Matrix(T)
{
	abstract size_t height() const @property;

	abstract size_t width() const @property;

	abstract ref const(T) opIndex(size_t i, size_t j) const;

	final ref const(T) opIndex(size_t i) const
	{
		if(width == 1)
			return this[i,0];
		if(height == 1)
			return this[0,i];

		throw new Exception("not a vector");
	}

	final override string toString() const @property
	{
		string s;
		auto strings = Array2!string(height, width);
		auto pitch = Array!size_t(width, 0);

		for(size_t i = 0; i < height; ++i)
			for(size_t j = 0; j < width; ++j)
			{
				strings[i,j] = to!string(this[i,j]);
				pitch[j] = max(pitch[j], strings[i,j].length);
			}

		for(size_t i = 0; i < height; ++i)
		{
			if(i == 0)
				s ~= "⎛ ";
			else if(i == height-1)
				s ~= "⎝ ";
			else
				s ~= "⎜ ";


			for(size_t j = 0; j < width; ++j)
			{
				s ~= strings[i,j];
				for(int k = 0; k < pitch[j]+1-strings[i,j].length; ++k)
					s ~= " ";
			}

			if(i == 0)
				s ~= "⎞\n";
			else if(i == height-1)
				s ~= "⎠";
			else
				s ~= "⎟\n";
		}
		return s;
	}

	/** return squared L2 norm */
	final RealTypeOf!T sqNorm() const @property
	{
		RealTypeOf!T sum = 0;
		for(size_t j = 0; j < width; ++j)
			for(size_t i = 0; i < height; ++i)
				sum += sqAbs(this[i,j]);
		return sum;
	}

	/** return L2 norm */
	final RealTypeOf!T norm() const @property
	{
		return std.math.sqrt(sqNorm);
	}

	final Matrix opBinary(string op)(Matrix b) const
		if(op == "+" || op == "-")
	{
		if(width != b.width || height != b.height)
			throw new Exception("matrix dimension mismatch");

		auto a = Slice2!T(height, b.width);
		for(size_t j = 0; j < b.width; ++j)
			for(size_t i = 0; i < height; ++i)
			{
				a[i,j] = mixin("this[i,j] "~op~" b[i,j]");
			}
		return Matrix(a.assumeUnique);
	}

	final Matrix opBinary(string op)(Matrix b) const
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
		return Matrix(a.assumeUnique);
	}

	final Matrix pow(long exp)
	{
		assert(width > 0);
		if(width != height)
			throw new Exception("can not take a power of a non-square matrix");
		if(exp <= 0)
			throw new Exception("non-positive powers not implemented (yet?)");

		Matrix r;
		Matrix base = this;

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
	DenseLU!T lu()
	{
		return new DenseLU!T(this);
	}

	/** compute QR decomposition */
	DenseQR!T qr()
	{
		return new DenseQR!T(this);
	}

	DenseHessenberg!T hessenberg()
	{
		return new DenseHessenberg!T(this);
	}

	static auto opCall(Slice!(immutable(T), 2) data)
	{
		return new DenseMatrix!T(data);
	}

	static auto opCall(size_t height, size_t width, immutable(T)[] data)
	{
		return new DenseMatrix!T(Slice!(immutable(T), 2)(height, width, data));
	}

	static DenseMatrix!T random(size_t height, size_t width)
	{
		static if(isFloatingPoint!T)
			return build!((i,j)=> cast(T)uniform(-1.0, 1.0))(height, width);
		else static if(isComplex!T)
			return build!((i,j)=> T(uniform(-1.0,1.0), uniform(-1.0,1.0)))(height, width);
		else
			return build!((i,j)=> T.random())(height, width);
	}

	static DenseMatrix!T random(Ring)(size_t height, size_t width, Ring ring)
	{
		auto data = Slice2!T(height, width);
		for(size_t j = 0; j < width; ++j)
			for(size_t i = 0; i < height; ++i)
				data[i,j] = ring.random();
		return new DenseMatrix!T(data.assumeUnique);
	}

	static DenseMatrix!T build(alias fun)(size_t h, size_t w)
	{
		auto data = Slice2!T(h, w);
		for(size_t j = 0; j < w; ++j)
			for(size_t i = 0; i < h; ++i)
				data[i,j] = binaryFun!(fun,"i","j")(i, j);
		return new DenseMatrix!T(data.assumeUnique);
	}

	/** only explicitly genrates upper/right half */
	static DenseMatrix!T buildSymmetric(alias fun)(size_t n)
	{
		auto data = Slice2!T(n, n);
		for(size_t j = 0; j < n; ++j)
			for(size_t i = 0; i <= j; ++i)
			{
				data[i,j] = binaryFun!(fun,"i","j")(i, j);
				if(i != j)
					data[j,i] = data[i,j];
			}
		return new DenseMatrix!T(data.assumeUnique);
	}

	static BandMatrix!T buildBand(alias fun)(size_t n, int kl, int ku)
	{
		assert(kl >= 0 && ku >= 0);
		auto data = Slice2!T(kl+ku+1, n);
		for(size_t j = 0; j < n; ++j)
			for(size_t i = max(ku,j)-ku; i <= min(n-1,j+kl); ++i)
				data[i+ku-j,j] = binaryFun!(fun,"i","j")(i, j);
		return new BandMatrix!T(kl, ku, data.assumeUnique);
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
		return new BandMatrix!T(kd, kd, data.assumeUnique);
	}

	static BandMatrix!T buildIdentity(size_t n)
	{
		return  buildBand!"i==j?1:0"(n, 0, 0);
	}

	Slice2!T dup() const @property
	{
		auto r = Slice2!T(height, width);
		foreach(i,j, ref x; r)
			x = this[i,j];
		return r;
	}

	/** dup with row permutation */
	Slice2!T dup(const int[] p) const @property
	{
		auto r = Slice2!T(p.length, width);
		foreach(i,j, ref x; r)
			x = this[p[i],j];
		return r;
	}

	static T zero = 0;
	static T one = 1;
}

/**
 * simple storage of a (possibly non-square) matrix.
 * column major order in a contigous array.
 */
final class DenseMatrix(T) : Matrix!T
{
	// TODO: decide if this should be Array2 instead of Slice2
	private Slice2!(immutable(T)) data;

	this(Slice2!(immutable(T)) data)
	{
		this.data = data;
	}

	override size_t height() const @property
	{
		return data.size[0];
	}

	override size_t width() const @property
	{
		return data.size[1];
	}

	override ref const(T) opIndex(size_t i, size_t j) const
	{
		return data[i,j];
	}

	DenseMatrix transpose() const
	{
		auto r = new DenseMatrix(data);
		swap(r.data.size[0], r.data.size[1]);
		swap(r.data.pitch[0], r.data.pitch[1]);
		return r;
	}

	DenseMatrix adjoint() const
	{
		auto r = this.dup;
		swap(r.size[0], r.size[1]);
		swap(r.pitch[0], r.pitch[1]);
		static if(isComplex!T)
			foreach(i, j, ref x; r)
				x = conj(x);
		return new DenseMatrix(r.assumeUnique);
	}
}

/** efficient storage of a square band matrix */
final class BandMatrix(T) : Matrix!T
{
	const int kl, ku; // number of sub/super-diagonals

	private Slice2!(immutable(T)) data;

	this(int kl, int ku, Slice2!(immutable(T)) data)
	{
		this.data = data;
		this.kl = kl;
		this.ku = ku;

		if(kl+ku+1 != data.size[0] || kl >= height || ku >= height)
			throw new Exception("invalid dimensions of band matrix");
	}

	override size_t height() const @property
	{
		return data.size[1]; // width == height
	}

	override size_t width() const @property
	{
		return data.size[1];
	}

	override ref const(T) opIndex(size_t i, size_t j) const
	{
		if(i+ku < j || i > j+kl)
			return zero;

		return data[i-j+ku,j];
	}
}

/**
 * lower/upper triangular matrix with/without implicit 1's on the diagonal
 */
final class TriangularMatrix(T, bool lower, bool implicitOne, int offDiag = 0) : Matrix!T
{
	private Slice2!(const(T)) data;

	this(Slice2!(const(T)) data)
	{
		this.data = data;
	}

	override size_t height() const @property
	{
		return data.size[0];
	}

	override size_t width() const @property
	{
		return data.size[1];
	}

	override ref const(T) opIndex(size_t i, size_t j) const
	{
		if(implicitOne && i == j)
			return one;

		if((lower && i+offDiag < j) || (!lower && i > j+offDiag))
			return zero;

		return data[i,j];
	}
}

final class PermutationMatrix(T) : Matrix!T
{
	immutable int[] p;

	this(immutable(int)[] p)
	{
		this.p = p;
	}

	override size_t height() const @property
	{
		return p.length;
	}

	override size_t width() const @property
	{
		return p.length;
	}

	override ref const(T) opIndex(size_t i, size_t j) const
	{
		if(p[i] == j)
			return one;
		else
			return zero;
	}
}

final class DenseLU(T)
{
	Slice2!T m;
	int[] p;

	this(Matrix!T _m)
	{

		p = new int[_m.height];
		m = _m.dup;
		denseComputeLU!T(m, p);
	}

	Matrix!T solve(Matrix!T b)
	{
		auto rhs = b.dup(p);
		denseSolveLU!T(m, rhs);
		return Matrix!T(rhs.assumeUnique);
	}

	auto L()
	{
		return new TriangularMatrix!(T, true, true)(m);
	}

	auto U()
	{
		return new TriangularMatrix!(T, false, false)(m);
	}

	auto P()
	{
		return new PermutationMatrix!T(cast(immutable(int)[])p);
	}
}

final class DenseQR(T)
{
	Slice2!T m;
	RealTypeOf!T[] beta;

	this(Matrix!T _m)
	{
		beta = new RealTypeOf!T[_m.height];
		m = _m.dup;
		denseComputeQR!T(m, beta);
	}

	Matrix!T solve(Matrix!T b)
	{
		auto rhs = b.dup();
		denseSolveQR!T(m, beta, rhs);
		return Matrix!T(rhs.assumeUnique);
	}

	auto Q()
	{
		auto q = Slice2!T(beta.length, beta.length);
		foreach(i, j, ref x; q)
			x = i==j ? 1 : 0;

		for(int k = cast(int)beta.length-1; k >= 0; --k)
			householderReflection!T(q[k..$,0..$], m[k+1..$,k], beta[k]);

		return Matrix!T(q.assumeUnique);
	}

	auto R()
	{
		return new TriangularMatrix!(T, false, false)(m);
	}
}

final class DenseHessenberg(T)
{
	Slice2!T m;
	RealTypeOf!T[] beta;

	this(Matrix!T _m)
	{
		beta = new RealTypeOf!T[_m.height];
		m = _m.dup;
		denseComputeHessenberg!T(m, beta);
	}

	auto Q()
	{
		auto q = Slice2!T(beta.length, beta.length);
		foreach(i, j, ref x; q)
			x = i==j ? 1 : 0;

		for(int k = cast(int)beta.length-2; k >= 0; --k)
			householderReflection!T(q[k+1..$,0..$], m[k+2..$,k], beta[k]);

		return Matrix!T(q.assumeUnique);
	}

	auto H()
	{
		return new TriangularMatrix!(T, false, false, 1)(m);
	}
}
