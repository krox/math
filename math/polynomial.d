module math.polynomial;

private import std.algorithm : map, min, max;
private import std.conv : to;
private import std.array : join;
private import std.exception : assumeUnique;

private import math.integer;

struct Polynomial(T, string _x = "x")
{
	immutable(T)[] coeffs; // highest one is always != 0

	this(immutable(T)[] coeffs)
	{
		while(coeffs.length && coeffs[$-1] == 0)
			coeffs = coeffs[0..$-1];
		this.coeffs = coeffs;
	}

	/** power of highest non-zero term. -1 for the 0-polynomial */
	int degree() const @property
	{
		return cast(int)coeffs.length-1;
	}

	string toString() const @property
	{
		if(coeffs.length == 0)
			return "0";

		string s;
		foreach_reverse(i, c; coeffs)
			if(c != 0)
			{
				if(s)
					s ~= " + ";
				if(i == 0 || c != 1)
					s ~= to!string(c);
				if(i >= 1)
					s ~= _x;
				if(i >= 2)
					s ~= "^" ~ to!string(i);
			}
		return s;
	}

	Polynomial opBinary(string op)(T b) const
		if(op == "*")
	{
		if(b == 0 || degree < 0)
			return Polynomial(null);
		if(b == 1)
			return this;

		auto r = new T[degree + 1];
		for(int i = 0; i <= degree; ++i)
			r[i] = b*coeffs[i];
		return Polynomial(assumeUnique(r));
	}

	Polynomial opBinary(string op)(T b) const
		if(op == "/")
	{
		return this*(b.inverse);
	}

	Polynomial opBinary(string op)(Polynomial b) const
		if(op == "+" || op == "-")
	{
		auto r = new T[max(coeffs.length, b.coeffs.length)];
		for(size_t i = 0; i < r.length; ++i)
		{
			if(i < coeffs.length && i < b.coeffs.length)
				r[i] = mixin("this.coeffs[i]"~op~"b.coeffs[i]");
			else if(i < coeffs.length)
				r[i] = this.coeffs[i];
			else
				static if(op == "+")
					r[i] = b.coeffs[i];
				else
					r[i] = -b.coeffs[i];
		}
		return Polynomial(assumeUnique(r));
	}

	Polynomial opBinary(string op)(Polynomial b) const
		if(op == "*")
	{
		if(degree < 0 || b.degree < 0)
			return Polynomial(null);

		auto r = new T[degree + b.degree + 1];

		for(int k = 0; k < r.length; ++k)
		{
			r[k] = coeffs[max(0, k-b.degree)] * b.coeffs[k-max(0, k-b.degree)];
			for(int i = max(0, k-b.degree)+1; i <= min(degree, k); ++i)
				r[k] = r[k] + coeffs[i] * b.coeffs[k-i];
		}

		return Polynomial(assumeUnique(r));
	}

	Polynomial opBinary(string op)(int e) const
		if(op == "^^")
	{
		assert(e >= 1);
		Polynomial r = this;

		for(int i = bsr(e)-1; i >= 0; --i)
		{
			r = r*r;
			if(e & (1<<i))
				r = r * this;
		}
		return r;
	}

	/** returns leading coefficient */
	T leading() const @property
	{
		return coeffs[$-1];
	}

	void divRem(Polynomial b, ref Polynomial quot, ref Polynomial rem)
	{
		if(b.degree < 0)
			throw new Exception("tried to divide by 0-polynomial");

		if(this.degree < b.degree)
		{
			quot = Polynomial(null);
			rem = this;
			return;
		}

		auto q = new T[this.degree - b.degree + 1];
		auto r = coeffs.dup;

		for(int i = cast(int)r.length-1; i >= b.degree; --i)
		{
			auto c = q[i-b.degree] = r[i]/b.leading;
			for(int j = 0; j <= b.degree; ++j)
				r[i-b.degree+j] = r[i-b.degree+j] - c*b.coeffs[j];
		}
		rem = Polynomial(assumeUnique(r));
		quot = Polynomial(assumeUnique(q));
	}

	Polynomial opBinary(string op)(Polynomial b)
		if(op == "/")
	{
		if(b.degree < 0)
			throw new Exception("tried to divide by 0-polynomial");

		if(this.degree < b.degree)
			return Polynomial(null);

		auto q = new T[this.degree - b.degree + 1];
		auto r = coeffs.dup;

		for(int i = cast(int)r.length-1; i >= b.degree; --i)
		{
			auto c = q[i-b.degree] = r[i]/b.leading;
			for(int j = 0; j <= b.degree; ++j)
				r[i-b.degree+j] = r[i-b.degree+j] - c*b.coeffs[j];
		}

		return Polynomial(assumeUnique(q));
	}

	Polynomial opBinary(string op)(Polynomial b) const
		if(op == "%")
	{
		if(b.degree < 0)
			throw new Exception("tried to divide by 0-polynomial");

		if(this.degree < b.degree)
			return this;

		auto r = coeffs.dup;
		for(int i = cast(int)r.length-1; i >= b.degree; --i)
		{
			auto c = r[i]/b.leading;
			for(int j = 0; j <= b.degree; ++j)
				r[i-b.degree+j] = r[i-b.degree+j] - c*b.coeffs[j];
		}

		return Polynomial(assumeUnique(r));
	}

	/** divide polynomial by leading coefficient */
	Polynomial normalize() const
	{
		if(degree < 0) // in this case, leading() is not defined
			return this;
		return this / leading;
	}

	/** formal derivative. I.e. (a*x^n)' = a*n*x^(n-1) */
	Polynomial derivative() const @property
	{
		if(degree < 1)
			return Polynomial(null);

		auto r = new T[degree];
		for(int i = 1; i <= degree; ++i)
			r[i-1] = coeffs[i]*i;

		return Polynomial(assumeUnique(r));
	}

	bool opEquals(T r) const
	{
		if(r == 0)
			return degree < 0;
		return degree == 0 && coeffs[0] == r;
	}

	bool opEquals(Polynomial r) const
	{
		return coeffs[] == r.coeffs[];
	}

	bool isSquareFree() const @property
	{
		// TODO: check for which fields T this is actually true
		return gcd(this, this.derivative).degree == 0;
	}

	Polynomial powmod(Integer e, Polynomial mod) const @property
	{
		if(e < 1)
			throw new Exception("polynomial powers < 1 not implemented yet"); // TODO

		auto r = this % mod;
		for(size_t i = e.length-1; i > 0; --i)
		{
			r = (r*r)%mod;
			if(e[i-1])
				r = (r*this)%mod;
		}
		return r;
	}
}

Polynomial!(T, _x) gcd(T, string _x)(Polynomial!(T,_x) a, Polynomial!(T,_x) b)
{
	// TODO: make it non-recursive and with fewer allocations
	if(a.degree<0)
		return b;
	else
		return gcd(b%a, a);
}
