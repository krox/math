module math.rational;

private import std.string : toStringz;
private import std.algorithm : move, swap;

import math.integer;


/**
 * Rational numbers. Based on math.integer.
 */
struct Rational
{
	public Integer num, denom;
	// invariants:
	// denom > 0
	// gcd(num, denom) == 1

	/** constructor for given value */
	this(int v) pure nothrow
	{
		num = Integer(v);
		denom = Integer(1);
	}

	/** ditto */
	this(int n, int d) pure
	{
		this(Integer(n), Integer(d));
		// TODO: cancel factors before creating Integer
	}

	/** ditto */
	this(Integer n, Integer d) pure
	{
		if(d == 0)
			throw new Exception("rational with denominator = 0");
		auto g = gcd(n,d);
		if(g != 1)
		{
			n = n.divExact(g);
			d = d.divExact(g);
		}
		if(d < 0)
		{
			n = -n;
			d = -d;
		}
		num = n;
		denom = d;
	}

	/** ditto */
	this(Integer v) pure nothrow
	{
		num = v;
		denom = Integer(1);
	}

	/** ditto */
	this(string v) pure
	{
		assert(false, "FIXME");
	}

	string toString() const pure nothrow @property
	{
		if(denom == 1)
			return num.toString;
		else
			return num.toString ~ "/" ~ denom.toString;
	}

	double opCast(T)() const pure nothrow
		if(is(T == double))
	{
		return cast(double)num / cast(double)denom;
	}

	/** return -1 / 0 / +1, faster than actual compare */
	int sign() const pure nothrow
	{
		return num.sign;
	}

	Rational opUnary(string op)() const pure nothrow
		if(op == "-")
	{
		return Rational(-num, denom);
	}

	/** returns 1/this */
	Rational inverse() const pure
	{
		if(num == 0)
			throw new Exception("tried to invert rational 0");

		if(num < 0)
			return Rational(-denom, -num);
		else
			return Rational(denom, num);
	}

	Rational opBinary(string op)(int b) const pure
	{
		     static if(op == "+") return Rational(num + denom*b, denom);
		else static if(op == "-") return Rational(num - denom*b, denom);
		else static if(op == "*") return Rational(num*b, denom);
		else static if(op == "/") return Rational(num, denom*b);
		else static if(op == "^^") return Rational(num^^b, denom^^b);
		else static assert(false, "binary '"~op~"' is not defined");
	}

	Rational opBinary(string op)(Integer b) const pure
	{
		     static if(op == "+") return Rational(num + denom*b, denom);
		else static if(op == "-") return Rational(num - denom*b, denom);
		else static if(op == "*") return Rational(num*b, denom);
		else static if(op == "/") return Rational(num, denom*b);
		else static if(op == "^^") return Rational(num^^b, denom^^b);
		else static assert(false, "binary '"~op~"' is not defined");
	}

	Rational opBinary(string op)(Rational b) const pure
	{
		static if(op == "+")
			return Rational(num*b.denom + b.num*denom, denom*b.denom);
		else static if(op == "-")
			return Rational(num*b.denom - b.num*denom, denom*b.denom);
		else static if(op == "*")
			return Rational(num*b.num, denom*b.denom);
		else static if(op == "/")
			return Rational(num*b.denom, denom*b.num);
		else static assert(false, "binary '"~op~"' is not defined");
	}

	Rational opBinaryRight(string op)(int a) const pure
	{
		     static if(op == "+") return Rational(denom*a + num, denom);
		else static if(op == "-") return Rational(denom*a - num, denom);
		else static if(op == "*") return Rational(a*num, denom);
		else static if(op == "/") return Rational(a*denom, num);
		else static assert(false, "binary '"~op~"' is not defined");
	}

	Rational opBinaryRight(string op)(Integer a) const pure
	{
		     static if(op == "+") return Rational(denom*a + num, denom);
		else static if(op == "-") return Rational(denom*a - num, denom);
		else static if(op == "*") return Rational(a*num, denom);
		else static if(op == "/") return Rational(a*denom, num);
		else static assert(false, "binary '"~op~"' is not defined");
	}

	Rational opOpAssign(string op, T)(T b) pure
	{
		this = this.opBinary!op(b);
		return this;
	}

	/** round to integer towards -infinity */
	Integer floor() const pure
	{
		return num / denom;
	}

	bool opEquals(int b) const pure nothrow
	{
		return denom == 1 && num == b;
	}

	bool opEquals(Integer b) const pure nothrow
	{
		return denom == 1 && num == b;
	}

	bool opEquals(Rational b) const pure nothrow
	{
		return num == b.num && denom == b.denom; // both numbers need to be normalized for this
	}

	int opCmp(Rational b) const pure nothrow
	{
		return (num*b.denom).opCmp(denom*b.num); // denominators need to be positive for this
	}

	int opCmp(Integer b) const pure nothrow
	{
		return num.opCmp(denom*b);
	}

	int opCmp(int b) const pure nothrow
	{
		return num.opCmp(denom*b);
	}
}
