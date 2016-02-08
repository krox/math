module math.matrix;

private import std.traits;
private import std.conv : to;
private import std.exception : assumeUnique;
private import std.typecons;
private import std.complex;
private import std.functional : binaryFun;
private import std.algorithm : min, max;
private import std.random : uniform;

import jive.array;
private import math.lapacke, math.linear;

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
	final T norm2() const @property
	{
		T sum = 0;
		for(size_t j = 0; j < width; ++j)
			for(size_t i = 0; i < height; ++i)
				sum = sum + this[i,j]*this[i,j];
		return sum;
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

	/** solve the linear equations this * x = b */
	Matrix!T solve(Matrix!T b)
	{
		auto n = cast(int)this.width;
		if(this.height != n || b.height != n)
			throw new Exception("invalid matrix dimensions");

		// compute (P)LU decomposition
		auto m = this.dup;
		auto p = new int[this.height];
		denseComputeLU!T(m, p);

		// solve it
		auto rhs = b.dup(p);
		denseSolveLU!T(m, rhs);

		return Matrix!T(rhs.assumeUnique);
	}

	static if(is(T == float) || is(T == double) || is(T == Complex!float) || is(T == Complex!double))
	{

	/** compute (complex) eigenvalues of this */
	Array!(ComplexType!T) eigenvalues()
	{
		auto n = cast(int)this.width;
		if(this.height != n)
			throw new Exception("invalid matrix dimensions");

		auto a = this.dup;
		auto w = Array!(ComplexType!T)(n);
		int r;

		     static if(is(T == float))       r = my_LAPACKE_sgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a.ptr, n, w.ptr, null, n, null, n);
		else static if(is(T == double))	     r = my_LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a.ptr, n, w.ptr, null, n, null, n);
		else static if(is(T == Complex!float))  r = LAPACKE_cgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a.ptr, n, w.ptr, null, n, null, n);
		else static if(is(T == Complex!double)) r = LAPACKE_zgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a.ptr, n, w.ptr, null, n, null, n);
		else throw new Exception("unsupported matrix type");

		if(r != 0)
			throw new Exception("lapack error");

		return w;
	}

	/** compute (real) eigenvalues of this assuming matrix is hermitian */
	Array!(RealType!T) hermitianEigenvalues()
	{
		auto n = cast(int)this.width;
		if(this.height != n)
			throw new Exception("invalid matrix dimensions");

		auto a = this.dup;
		auto w = Array!(RealType!T)(n);
		int r;

		     static if(is(T == float))          r = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'N', 'U', n, a.ptr, n, w.ptr);
		else static if(is(T == double))         r = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'N', 'U', n, a.ptr, n, w.ptr);
		else static if(is(T == Complex!float))  r = LAPACKE_cheev(LAPACK_COL_MAJOR, 'N', 'U', n, a.ptr, n, w.ptr);
		else static if(is(T == Complex!double)) r = LAPACKE_zheev(LAPACK_COL_MAJOR, 'N', 'U', n, a.ptr, n, w.ptr);
		else throw new Exception("unsupported matrix type");

		if(r != 0)
			throw new Exception("lapack error");

		return w;
	}

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
			return T(0);

		return data[i-j+ku,j];
	}

	/** solve the linear equations this * x = b */
	override Matrix!T solve(Matrix!T b)
	{
		auto n = cast(int)this.width;
		auto nrhs = cast(int)b.width;
		if(this.height != n || b.height != n)
			throw new Exception("invalid matrix dimensions");

		auto a = this.data.data.dup;
		auto rhs = b.dup;
		auto p = new int[n];
		int r;

		     static if(is(T == float))          r = LAPACKE_sgbsv(LAPACK_COL_MAJOR, n, kl, ku, nrhs, a.ptr, 1+kl+ku, p.ptr, rhs.ptr, n);
		else static if(is(T == double))         r = LAPACKE_dgbsv(LAPACK_COL_MAJOR, n, kl, ku, nrhs, a.ptr, 1+kl+ku, p.ptr, rhs.ptr, n);
		else static if(is(T == Complex!float))  r = LAPACKE_cgbsv(LAPACK_COL_MAJOR, n, kl, ku, nrhs, a.ptr, 1+kl+ku, p.ptr, rhs.ptr, n);
		else static if(is(T == Complex!double)) r = LAPACKE_zgbsv(LAPACK_COL_MAJOR, n, kl, ku, nrhs, a.ptr, 1+kl+ku, p.ptr, rhs.ptr, n);
		else throw new Exception("unsupported matrix type");

		if(r != 0)
			throw new Exception("lapack error");

		return Matrix!T(rhs.assumeUnique);
	}

	/** compute (complex) eigenvalues of this */
	override Array!(ComplexType!T) eigenvalues()
	{
		throw new Exception("non-symmetric eigenvalues on band matrices not supported");
	}

	/** compute (real) eigenvalues of this, assuming matrix is hermitian */
	override Array!(RealType!T) hermitianEigenvalues()
	{
		auto n = cast(int)this.width;
		if(this.height != n)
			throw new Exception("invalid matrix dimensions");

		auto a = this.data.data.dup;
		auto w = Array!(RealType!T)(n);
		int r;

	         static if(is(T == float))          r = LAPACKE_ssbev(LAPACK_COL_MAJOR, 'N', 'U', n, ku, a.ptr, 1+kl+ku, w.ptr, null, n);
		else static if(is(T == double))         r = LAPACKE_dsbev(LAPACK_COL_MAJOR, 'N', 'U', n, ku, a.ptr, 1+kl+ku, w.ptr, null, n);
		else static if(is(T == Complex!float))  r = LAPACKE_chbev(LAPACK_COL_MAJOR, 'N', 'U', n, ku, a.ptr, 1+kl+ku, w.ptr, null, n);
		else static if(is(T == Complex!double)) r = LAPACKE_zhbev(LAPACK_COL_MAJOR, 'N', 'U', n, ku, a.ptr, 1+kl+ku, w.ptr, null, n);
		else throw new Exception("unsupported matrix type");

		if(r != 0)
			throw new Exception("lapack error");

		return w;
	}
}
