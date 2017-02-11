module math.polynomial;

private import std.algorithm : map, min, max;
private import std.conv : to;
private import std.exception : assumeUnique;
private import std.functional : unaryFun;
private import core.bitop : bsr;
private import math.integer;
private import math.numtheory : binomial;

/**
 * univariate polynomials
 *  - T is the type of coefficients
 *  - x is the same of the variable (currently only used for toString)
 *  - delta can be set for skew-polynomial rings, i.e.
 *    x*a = a*x + delta(a) for any coefficient a.
 *
 * major usecases
 *  - Polynomials over finite fields: Polynomial!FFE
 *    (see modules math.finitefields and math.polyfactor)
 *  - Differential operators with polynomial functions:
 *    Polynomial!(Polynomial!double, "D", "a.derivative")
 *
 */
struct Polynomial(T, string _x = "x", alias _delta = "0")
{
	//////////////////////////////////////////////////////////////////////
	/// constructors and internals
	//////////////////////////////////////////////////////////////////////

	enum bool skew = _delta != "0";
	alias delta = unaryFun!_delta;

	immutable(T)[] coeffs; // highest one is always != 0
	alias coeffs this;

	this(int c) pure
	{
		this(T(c));
	}

	this(immutable(T)[] coeffs...) pure
	{
		while(coeffs.length && coeffs[$-1] == 0)
			coeffs = coeffs[0..$-1];
		this.coeffs = coeffs;
	}

	static Polynomial random(Field)(int d, Field field) @property
	{
		if(d < 0)
			return Polynomial(null);
		auto coeffs = new T[d+1];
		for(size_t i = 0; i < coeffs.length; ++i)
			coeffs[i] = field.random;
		return Polynomial(assumeUnique(coeffs));
	}

	//////////////////////////////////////////////////////////////////////
	/// basic properties
	//////////////////////////////////////////////////////////////////////

	/** power of highest non-zero term. -1 for the 0-polynomial */
	int degree() const pure nothrow @property
	{
		return cast(int)coeffs.length-1;
	}

	/** coefficient of x^i */
	T opIndex(size_t i) const pure nothrow @property
	{
		if(i < 0 || i >= coeffs.length)
			return T(0);
		return coeffs[i];
	}

	/** returns leading coefficient */
	T leading() const pure nothrow @property
	{
		return coeffs[$-1];
	}

	static if(is(typeof(1/T.init)))
	{
		Polynomial normalize() const pure @property
		{
			return (1/leading)*this;
		}
	}

	/** human readable string */
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

	bool opEquals(int r) const pure
	{
		if(r == 0)
			return degree < 0;
		return degree == 0 && coeffs[0] == r;
	}

	bool opEquals(T r) const pure
	{
		if(r == 0)
			return degree < 0;
		return degree == 0 && coeffs[0] == r;
	}

	bool opEquals(Polynomial r) const pure
	{
		return this.coeffs[] == r.coeffs[];
	}

	//////////////////////////////////////////////////////////////////////
	/// additive arithmetic
	//////////////////////////////////////////////////////////////////////

	/** - polynomial */
	Polynomial opUnary(string op)() const pure
		if(op == "-")
	{
		auto r = new T[degree+1];
		for(int n = 0; n <= degree; ++n)
			r[n] = -coeffs[n];
		return Polynomial(assumeUnique(r));
	}

	/** polynomial +- polynomial */
	Polynomial opBinary(string op)(Polynomial b) const pure
		if(op == "+" || op == "-")
	{
		auto r = new T[max(degree, b.degree)+1];
		for(size_t i = 0; i < r.length; ++i)
			r[i] = mixin("this[i]"~op~"b[i]");
		return Polynomial(assumeUnique(r));
	}

	//////////////////////////////////////////////////////////////////////
	/// multiplicative arithmetic
	//////////////////////////////////////////////////////////////////////

	/** scalar * polynomial */
	Polynomial opBinaryRight(string op, S)(S b) const pure
		if(op == "*" && (is(S == T) || is(S == int)))
	{
		if(b == 0 || degree < 0)
			return Polynomial(null);
		if(b == 1)
			return this;

		auto r = new T[degree + 1];
		for(int n = 0; n <= degree; ++n)
			r[n] = b*coeffs[n];
		return Polynomial(assumeUnique(r));
	}

	/** polynomial * scalar */
	Polynomial opBinary(string op, S)(S b) const pure
		if(op == "*" && (is(S == T) || is(S == int)))
	{
		if(b == 0 || degree < 0)
			return Polynomial(null);
		if(b == 1)
			return this;

		auto r = new T[degree + 1];

		static if(skew)
		{
			r[] = T(0);
			for(int k = 0; k <= degree; ++k)
				for(int n = k; n <= degree; ++n)
				{
					r[n-k] = r[n-k] + T(cast(int)binomial(n,k))*coeffs[n]*b;
					b = delta(b);
				}
		}
		else
		{
			for(int n = 0; n <= degree; ++n)
				r[n] = coeffs[n]*b;
		}

		return Polynomial(assumeUnique(r));
	}

	/** polynomial / scalar */
	Polynomial opBinary(string op, S)(S b) const pure
		if(op == "/" && (is(S == T) || is(S == int)))
	{
		return this*(1/b);
	}

	/** polynomial * polynomial */
	Polynomial opBinary(string op)(Polynomial b) const pure
		if(op == "*")
	{
		if(degree < 0 || b.degree < 0)
			return Polynomial(null);

		auto r = new T[degree + b.degree + 1];
		r[] = T(0);

		static if(skew)
		{
			for(int m = 0; m <= b.degree; ++m)
			{
				T tb = b.coeffs[m];
				for(int k = 0; k <= degree; ++k)
				{
					for(int n = k; n <= degree; ++n)
						r[n+m-k] = r[n+m-k] + T(cast(int)binomial(n,k))*coeffs[n]*tb;
					tb = delta(tb);
				}
			}
		}
		else
		{
			for(int n = 0; n <= degree; ++n)
				for(int m = 0; m <= b.degree; ++m)
					r[n+m] += coeffs[n]*b.coeffs[m];
		}

		return Polynomial(assumeUnique(r));
	}

	/** polynomial ^ integer */
	Polynomial opBinary(string op)(int e) const pure
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


	//////////////////////////////////////////////////////////////////////
	/// modular arithmetic
	//////////////////////////////////////////////////////////////////////

	// actually, modular arithmetic for skew polynomials is possible if you
	// distinguish between left- and right- gcd and so on.
	static if(!skew)
	{
		void divRem(Polynomial b, ref Polynomial quot, ref Polynomial rem) pure
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

		/* don't use with numerical types */
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
}

Polynomial!(T, _x) gcd(T, string _x)(Polynomial!(T,_x) a, Polynomial!(T,_x) b)
{
	// TODO: make it non-recursive and with fewer allocations
	if(a.degree<0)
		return b;
	else
		return gcd(b%a, a);
}
