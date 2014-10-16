module math.matrix;

private import std.traits;
private import std.conv : to;
private import std.exception : assumeUnique;

import jive.slice;

struct Matrix(T)
{
	private Slice!(const(T),2) data; // TODO: replace "const" with "immutable"

	this(Slice!(const(T),2) data)
	{
		this.data = data;
	}

	this(size_t height, size_t width, const(T)[] data)
	{
		this(Slice!(const(T),2)(height,width,data));
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
		return Matrix(height, width, a);
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
		return Matrix(a.toConst);
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
