module math.numberfield;

// NOTE: try not to include any GMP dependent module here

import std.format;

/**
 * Numbers of the form (a + b * sqrt(d)), with a, b and d a (fixed) non-square.
 * should (at least) work with T == Rational and T == IntMod
 */
struct Quadratic(T)
{
	T a, b, d;

	this(T a, T b, T d)
	{
		// would be nice to check d for non-square here

		this.a = a;
		this.b = b;
		this.d = d;
	}

	string toString() const pure @property
	{
		return format("%s+%s√%s", a, b, d);
	}

	/** (additive) inverse) */
	Quadratic opUnary(string op)() const pure nothrow @property
		if(op == "-")
	{
		return Quadratic(-a, -b, d);
	}

	/** negate the 'imaginary' part */
	Quadratic conjugate() const pure nothrow
	{
		return Quadratic(a, -b, d);
	}

	/** returns the norm N(x) = x*conj(x) */
	T norm() const pure nothrow @property
	{
		return a*a-b*b*d;
	}

	/** returns 1/this */
	Quadratic inverse() const pure nothrow
	{
		return conjugate / norm;
	}

	Quadratic opBinary(string op)(T rhs) const pure nothrow
	{
		static if(op == "+")
			return Quadratic(a+rhs, b, d);
		else static if(op == "-")
			return Quadratic(a-rhs, b, d);
		else static if(op == "*")
			return Quadratic(a*rhs, b*rhs, d);
		else static if(op == "/")
			return Quadratic(a/rhs , b/rhs, d);
		else static assert(false, "binary assign '"~op~"' is not defined");
	}

	Quadratic opBinary(string op)(Quadratic rhs) const pure nothrow
	{
		assert(this.d == rhs.d);

		     static if(op == "+") return Quadratic(a+rhs.a, b+rhs.b, d);
		else static if(op == "-") return Quadratic(a-rhs.a, b-rhs.b, d);
		else static if(op == "*") return Quadratic(a*rhs.a + b*rhs.b*d, a*rhs.b + b*rhs.a, d);
		else static if(op == "/") return Quadratic(a*rhs.a - b*rhs.b*d, b*rhs.a - a*rhs.b, d) / rhs.norm;
		else static assert(false, "binary assign '"~op~"' is not defined");
	}

	Quadratic opBinary(string op)(long e) const pure nothrow
		if(op == "^^")
	{
		Quadratic b = this;
		if(e < 0)
		{
			b = b.inverse;
			e = -e;
		}

		auto r = this/this; // TODO: remove hack
		for(; e != 0; e >>= 1, b = b*b)
			if(e&1)
				r = r * b;
		return r;
	}

	bool opEquals(Quadratic r) const pure nothrow
	{
		assert(d == r.d);
		return a == r.a && b == r.b;
	}

	bool opEquals(int r) const pure nothrow
	{
		return a == r && b == 0;
	}

	/** returns largest integer <= this */
	/+static if(is(T == Rational))
	{
		Integer floor() const pure nothrow @property
		{
			if(d.sign != 1)
				throw new Exception("floor is not defined in imaginary quadratic fields");

			Integer x = a.denom*b.num;
			x = isqrt(x*x*d); // this root is never exact...

			if(a.denom.sign * b.num.sign == -1)
				x = (-1)-x; // ... therefore the "-1" is always necessary

			return (x + a.num * b.denom) / (a.denom*b.denom);
		}
	}+/
}

unittest
{
	import math.numtheory;
	alias Q = Quadratic!IntMod;
	long p = 17;
	auto d = IntMod(3, p); // thats a non-square so the quadratic field is GF(17^^2)
	auto x = Q(IntMod(1,p), IntMod(1,p), d);
	assert(x^^(p*p-2) == x.inverse);
}
