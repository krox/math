module math.matrix;

private import std.traits;
private import std.conv : to;
private import std.exception : assumeUnique;
private import std.typecons;
private import std.complex;

import jive.slice;
import jive.array;
private import math.lapacke;

struct Matrix(T)
{
	private Slice!(immutable(T), 2) data;

	this(Slice!(immutable(T), 2) data)
	{
		this.data = data;
	}

	this(size_t height, size_t width, immutable(T)[] data)
	{
		this(Slice!(immutable(T),2)(height,width,data));
	}

	size_t height() const @property
	{
		return data.size[0];
	}

	size_t width() const @property
	{
		return data.size[1];
	}

	T opIndex(size_t i, size_t j) const
	{
		return data[i,j];
	}

	T opIndex(size_t i) const
	{
		if(width == 1)
			return data[i,0];
		if(height == 1)
			return data[0,i];

		throw new Exception("not a vector");
	}

	string toString() const @property
	{
		string s;
		for(size_t i = 0; i < height; ++i)
		{
			if(i == 0)
				s ~= "⎛ ";
			else if(i == height-1)
				s ~= "⎝ ";
			else
				s ~= "⎜ ";


			for(size_t j = 0; j < width; ++j)
				s ~= to!string(this[i,j]) ~ " ";

			if(i == 0)
				s ~= "⎞\n";
			else if(i == height-1)
				s ~= "⎠";
			else
				s ~= "⎟\n";
		}
		return s;
	}

	Matrix opBinary(string op)(Matrix b) const
		if(op == "+" || op == "-")
	{
		if(width != b.width || height != b.height)
			throw new Exception("matrix dimension mismatch");
		T[] a = new T[width*height];
		mixin("a[] = this.data[]"~op~"b.data[];");
		return Matrix(height, width, a.assumeUnique);
	}

	Matrix opBinary(string op)(Matrix b) const
		if(op == "*")
	{
		if(width != b.height)
			throw new Exception("matrix dimension mismatch");

		auto a = Slice!(T,2)(height, b.width);
		for(size_t i = 0; i < height; ++i)
			for(size_t j = 0; j < b.width; ++j)
			{
				a[i,j] = this[i,0]*b[0,j];
				for(size_t k = 1; k < width; ++k)
					a[i,j] = a[i,j] + this[i,k]*b[k,j];
			}
		return Matrix(a.assumeUnique);
	}

	Matrix pow(long exp) const
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
	Matrix solve(Matrix!T b)
	{
		auto n = cast(int)this.width;
		auto nrhs = cast(int)b.width;
		if(this.height != n || b.height != n)
			throw new Exception("invalid matrix dimensions");

		auto a = this.data.dup;
		auto rhs = b.data.dup;
		auto p = new int[n];

		int r;

		static if(is(T == float))
			r = LAPACKE_sgesv(LAPACK_COL_MAJOR, n, nrhs, a.data.ptr, n, p.ptr, rhs.data.ptr, n);
		else static if(is(T == double))
			r = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, a.data.ptr, n, p.ptr, rhs.data.ptr, n);
		else static if(is(T == Complex!float))
			r = LAPACKE_cgesv(LAPACK_COL_MAJOR, n, nrhs, a.data.ptr, n, p.ptr, rhs.data.ptr, n);
		else static if(is(T == Complex!double))
			r = LAPACKE_zgesv(LAPACK_COL_MAJOR, n, nrhs, a.data.ptr, n, p.ptr, rhs.data.ptr, n);
		else throw new Exception("unsupported matrix type");

		if(r != 0)
			throw new Exception("lapack error");

		return Matrix(rhs.assumeUnique);
	}

	/** compute eigenvalues of this */
	auto eigenvalues()
	{
		auto n = cast(int)this.width;
		if(this.height != n)
			throw new Exception("invalid matrix dimensions");

		auto a = this.data.dup;

		int r;

		static if(is(T == float) || is(T == double))
		{
			auto wr = new T[n];
			auto wi = new T[n];

			static if(is(T == float))
				r = LAPACKE_sgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a.ptr, n, wr.ptr, wi.ptr, null, n, null, n);
			static if(is(T == double))
				r = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a.ptr, n, wr.ptr, wi.ptr, null, n, null, n);

			auto w = Array!(Complex!T)(n);
			for(int i = 0; i < n; ++i)
				w[i] = Complex!T(wr[i], wi[i]);
		}
		else static if(is(T == Complex!float) || is(T == Complex!double))
		{
			auto w = Array!T(n);
			static if(is(T == Complex!float))
				r = LAPACKE_cgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a.ptr, n, w.ptr, null, n, null, n);
			static if(is(T == Complex!double))
				r = LAPACKE_zgeev(LAPACK_COL_MAJOR, 'N', 'N', n, a.ptr, n, w.ptr, null, n, null, n);
		}
		else throw new Exception("unsupported matrix type");

		if(r != 0)
			throw new Exception("lapack error");

		return w;
	}
}


struct Mat(T, size_t N, size_t M)
{
	T[N][M] c;

	enum Mat zero = constant(0);
	enum Mat identity = _identity();

	static Mat constant(T v) pure
	{
		Mat r;
		foreach(n; 0..N)
			foreach(m; 0..M)
				r[n,m] = v;
		return r;
	}

	static Mat _identity() pure
	{
		Mat r;
		foreach(n; 0..N)
			foreach(m; 0..M)
				r[n,m] = n==m ? 1 : 0;
		return r;
	}

	ref T opIndex(size_t n, size_t m)
	{
		return c[m][n];
	}

	ref const(T) opIndex(size_t n, size_t m) const
	{
		return c[m][n];
	}

	static if(M == 1)
	{
		ref inout(T) opIndex(size_t n) inout
		{
			return c[0][n];
		}
	}

	static if(N == 1)
	{
		ref inout(T) opIndex(size_t m) inout
		{
			return c[m][0];
		}
	}

	void opAddAssign(const Mat b)
	{
		foreach(n; 0..N)
			foreach(m; 0..M)
				 this[n,m] += b[n,m];
	}

	void opSubAssign(const Mat b)
	{
		foreach(n; 0..N)
			foreach(m; 0..M)
				 this[n,m] -= b[n,m];
	}

	void opMulAssign(T s)
	{
		foreach(n; 0..N)
			foreach(m; 0..M)
				 this[n,m] *= s;
	}

	Mat opMul(T s) const
	{
		Mat r = this;
		r *= s;
		return r;
	}

	static if(isFloatingPoint!T)
	{
		void opDivAssign(T s)
		{
			opMulAssign(1.0/s);
		}

		Mat opDiv(T s) const
		{
			return opMul(1.0/s);
		}
	}

	Mat opAdd(const Mat b) const
	{
		Mat r = this;
		r += b;
		return r;
	}

	Mat!(T,N,K) opMul(size_t K)(ref const Mat!(T,M,K) rhs) const
	{
		auto r = Mat!(T,N,K).zero;
		foreach(n; 0..N)
			foreach(k; 0..K)
				foreach(m; 0..M)
					r[n,k] += this[n,m]*rhs[m,k];
		return r;
	}

	static if(N==M) Mat pow(long exp) const
	{
		assert(exp >= 0);
		Mat r = identity;
		Mat base = this;

		while(exp)
		{
			if(exp & 1)
				r = r * base;
			exp >>= 1;
			base = base * base;
		}
		return r;
	}

	string toString() const @property
	{
		return to!string(c);	// TODO: transpose this
	}
}

alias Vec(T, size_t N)  = Mat!(T, N, 1);

T cross(T)(Vec!(T,2) a, Vec!(T,2) b)
{
    return a[0] * b[1] - a[1] * b[0];
}

bool ccw(T)(Vec!(T,2) a, Vec!(T,2) b)
{
	return cross!T(a, b) > 0;
}

bool ccw(T)(Vec!(T,2) a, Vec!(T,2) b, Vec!(T,2) c)
{
	return cross(b-a, c-a) > 0;
}
