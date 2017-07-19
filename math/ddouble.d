/**
 * This module is based on work by Yozu Hida, Xiaoye S. Li and David H. Bailey
 * and others. A 2007 paper describing all relevant ideas is available at
 * http://www.jaist.ac.jp/~s1410018/papers/qd.pdf. Their C/C++/Fortran library
 * is available at http://crd-legacy.lbl.gov/~dhbailey/mpdist/ (BSD license).
 *
 * This module implements only a small subset of their algorithms. All code in
 * here was written by myself from scratch and released into public domain.
 */

/**
 * Remarks / Issues / TODO:
 *  - Correct propagation of infinity and NaN is not tested.
 *  - double-double operations during compile-time are broken. CTFE is
 *    explicitly disabled for the relevant functions to avoid problems.
 *  - Algorithms rely on certain setting of the FPU like rounding mode. Should
 *    be correct by default, but no checks are done.
 *  - The error of a sum a+b is of the order O(eps*(|a|+|b|)). For comparison
 *    when using IEEE-floating-point numbers or the MPFR library, the error
 *    would be of the order  O(eps*|a+b|). Similar weaker-than-IEEE error
 *    bounds hold for other operations. This is a design choice, as more
 *    precise algorithms are significantly slower and usually not necessary.
 *  - All operations are implemented as member functions, i.e. x.sqrt()
 *    instead of sqrt(x). This avoids collision with functions from std.math.
 *  - Basic operations are implemented without control-flow. This should make
 *    vectorization straight-forward if the need ever arises.
 *  - toString/fromString are not as fast or accurate as they could be because
 *    they rely on many multiplications/divisions by 10.
 */

module math.ddouble;

private import std.math;
private import std.traits;
private import std.format;
private import std.random;

//debug = ddouble;

/**
 * "double double" type which is a (non-IEEE) floating point type implemented
 * as the sum of two doubles and thus has about 107 bits of effective
 * mantissa. This should be significantly faster than an actual multi-precision
 * library like MPFR. So as long as quadruple precision is not supported by
 * hardware, this is the most efficient alternative in cases where only a
 * little more than double-precision is needed.
 */
struct ddouble
{
	//////////////////////////////////////////////////////////////////////
	/// internals and some constants
	//////////////////////////////////////////////////////////////////////

	private double hi;
	private double lo;

	enum nan = ddouble(double.nan, double.nan);
	enum infinity = ddouble(double.infinity, double.infinity);
	enum epsilon = ddouble(double.epsilon^^2, 0); // this might be over-optimistic

	enum pi = ddouble(3.141592653589793116e+00, 1.224646799147353207e-16);
	enum e = ddouble(2.718281828459045091e+00, 1.445646891729250158e-16);
	enum log2 = ddouble(6.931471805599452862e-01, 2.319046813846299558e-17);
	enum log10 = ddouble(2.302585092994045901e+00, -2.170756223382249351e-16);

	debug(ddouble) invariant
	{
		if(isNaN(hi) || isNaN(lo))
			return;

		// During CTFE, doubles are promoted to real for extra precision.
		// That makes double-double impossible.
		if(__ctfe)
			return;

		assert(hi + lo == hi);
	}


	//////////////////////////////////////////////////////////////////////
	/// constructors and basic casts
	//////////////////////////////////////////////////////////////////////

	private this(double hi, double lo) pure @safe nothrow
	{
		this.hi = hi;
		this.lo = lo;
	}

	/** constructor taking a float/double/int/... */
	this(S)(S x) pure @safe nothrow
	{
		if(__ctfe)
			assert(false, "compile-time double-double computations are not supported");
		this.hi = x;
		this.lo = x-hi;
	}

	/**
	 * Constructor taking a (decimal) string in the format
	 * [0-9]*(.[0-9]*)?((e|E)[0-9]+)?
	 */
	this(string s) pure @safe
	{
		this(0,0);
		size_t k = 0;

		// sign
		bool sign = false;
		if(k < s.length && s[k] == '-')
		{
			++k;
			sign = true;
		}

		// whole part
		for(; k < s.length && '0' <= s[k] && s[k] <= '9'; ++k)
			this = 10 * this + (s[k]-'0');

		// fractional part
		if(k < s.length && s[k] == '.')
		{
			++k;
			ddouble inc = 1;
			for(; k < s.length && '0' <= s[k] && s[k] <= '9'; ++k)
			{
				inc /= 10;
				this += inc * (s[k] - '0');
			}
		}

		// exponent
		if(k < s.length && (s[k] == 'e' || s[k] == 'E'))
		{
			++k;
			bool esign = false;
			if(k < s.length && s[k] == '-')
			{
				++k;
				esign = true;
			}
			int e;
			for(; k < s.length && '0' <= s[k] && s[k] <= '9'; ++k)
				e = 10 * e + (s[k]-'0');
			if(esign)
				this *= ddouble(10)^^-e;
			else
				this *= ddouble(10)^^e;
		}

		if(k != s.length)
			throw new Exception("invalid floating point string: '"~s~"'");
		if(sign)
		{
			hi = -hi;
			lo = -lo;
		}
	}

	/**
	 * generate a random number in the interval [0,1). Note that only 64 bit of
	 * entropy are used, i.e. the trailing ~43 bits of mantissa are always zero.
	 */
	static ddouble random()
	{
		return ddouble(uniform!ulong).ldexp(-64);
	}

	/** only scientific notation is supported right now */
	void toString(scope void delegate(const(char)[]) sink, FormatSpec!char fmt) const
	{
		// special cases and sign
		if(isNaN(hi))
			return sink("nan");
		if(hi < 0)
			sink("-");
		if(hi == double.infinity || hi == -double.infinity)
			sink("inf");
		if(hi == 0)
			return sink("0");
		assert(isFinite(hi) && isFinite(lo));

		// split abs(this) = x * 10^e
		auto e = cast(int)std.math.log10(std.math.abs(hi));
		auto x = this.abs * (ddouble(10) ^^ -e);
		if(x < 1)
		{
			x *= 10;
			e -= 1;
		}
		else if(x >= 10)
		{
			x /= 10;
			e += 1;
		}
		assert(1 <= x && x < 10);

		// write mantissa
		enum digits = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"];
		int prec = fmt.precision;
		if(prec == FormatSpec!char.UNSPECIFIED)
			prec = 7; // same default as for float/double
		assert(prec > 0);
		sink(digits[cast(int)x]);
		if(prec > 1)
			sink(".");
		for(int i = 1; i < prec; ++i)
		{
			x -= cast(int)x;
			x = 10*x;
			sink(digits[cast(int)x]);
		}

		// write exponent
		if(e != 0)
		{
			sink("e");
			sink(format("%s", e));
		}
	}

	/** default formatting */
	string toString() const
	{
		return format("%s");
	}

	/** casts to builtin types */
	T opCast(T)() const pure @safe nothrow
	{
		static if(isFloatingPoint!T)
			return hi;
		else static if(isIntegral!T)
			return cast(T)hi; // TODO: not correct for (u)long
		else static if(is(T == string))
			return toString();
		else static assert(false);
	}

	void opAssign(double x) pure @safe nothrow
	{
		hi = x;
		lo = 0;
	}


	//////////////////////////////////////////////////////////////////////
	/// arithmetic
	//////////////////////////////////////////////////////////////////////

	/** +- ddouble */
	ddouble opUnary(string op)() const pure @safe nothrow
	{
		static if(op == "+")
			return ddouble(hi, lo);
		else static if(op == "-")
			return ddouble(-hi, -lo);
		else static assert(false);
	}

	/** ddouble <-> double */
	ddouble opBinary(string op)(double b) const pure @safe nothrow
	{
		static if(op == "+")
		{
			double rhi, rlo;
			twoSum(hi, b, rhi, rlo);
			twoSumQuick(rhi, rlo + lo, rhi, rlo);
			return ddouble(rhi, rlo);
		}
		else static if(op == "-")
		{
			return this + (-b);
		}
		else static if(op == "*")
		{
			double rhi, rlo;
			twoProd(hi, b, rhi, rlo);
			twoSumQuick(rhi, rlo + lo*b, rhi, rlo);
			return ddouble(rhi, rlo);
		}
		else static if(op == "/")
		{
			double q1, q2;
			double p1, p2;
			double s, e;

			q1 = hi / b;
			twoProd(q1, b, p1, p2);
			twoSum(hi, -p1, s, e);
			q2 = (s + (e + lo - p2)) / b;
			twoSumQuick(q1, q2, q1, q2);
			return ddouble(q1, q2);
		}
		else static assert(false);
	}

	/** double <-> ddouble */
	ddouble opBinaryRight(string op)(double a) const pure @safe nothrow
	{
		static if(op == "+")
			return this + a;
		else static if(op == "-")
			return (-this) + a;
		else static if(op == "*")
			return this * a;
		else static if(op == "/")
			return ddouble(a) / this; // TODO: implement special case
		else static assert(false);
	}

	/** ddouble <-> ddouble */
	ddouble opBinary(string op)(ddouble b) const pure @safe nothrow
	{
		static if(op == "+")
		{
			double s1, s2;
			twoSum(hi, b.hi, s1, s2);
			twoSumQuick(s1, s2 + (lo + b.lo), s1, s2);
			return ddouble(s1, s2);

		}
		else static if(op == "-")
		{
			return this + (-b);
		}
		else static if(op == "*")
		{
			double p1, p2;
			twoProd(hi, b.hi, p1, p2);
			twoSumQuick(p1, p2 + (hi * b.lo + lo * b.hi), p1, p2);
			return ddouble(p1, p2);
		}
		else static if(op == "/")
		{
			double s1, s2;
			double q1, q2;
			q1 = hi / b.hi;
			ddouble r = b * q1;
			twoSum(hi, -r.hi, s1, s2);
			s2 -= r.lo;
			s2 += lo;
			q2 = (s1 + s2) / b.hi;
			twoSumQuick(q1, q2, q1, q2);
			return ddouble(q1, q2);
		}
		else static assert(false);
	}

	/** ddouble ^^ long */
	ddouble opBinary(string op)(long n) const pure @safe nothrow
		if(op == "^^")
	{
		if(n == 0)
			return ddouble(1);

		ddouble a = this;
		if(n < 0)
		{
			a = a.inverse;
			n = -n;
		}

		ddouble r = ddouble(1);
		for(; n != 0; n >>= 1, a = a.sqr)
			if(n & 1)
				r *= a;
		return r;
	}

	/** opOpAssign for convenience */
	ddouble opOpAssign(string op, S)(S b) pure @safe nothrow
	{
		this = this.opBinary!op(b);
		return this;
	}

	/** 1/this */
	ddouble inverse() const pure @safe nothrow
	{
		double s1, s2;
		double q1, q2;

		q1 = 1 / hi;
		ddouble r = this * q1;
		twoSum(1, -r.hi, s1, s2);
		s2 -= r.lo;
		s2 += 0;

		/* get next approximation */
		q2 = (s1 + s2) / hi;

		/* renormalize */
		twoSumQuick(q1, q2, q1, q2);
		return ddouble(q1, q2);
	}

	/** absolute value */
	ddouble abs() const pure @safe nothrow
	{
		if(this < 0)
			return -this;
		else
			return this;
	}

	/** efficiently (and accurately) calculate this*2^n */
	ddouble ldexp(int n) const pure @safe nothrow
	{
		return ddouble(std.math.ldexp(hi, n), std.math.ldexp(lo, n));
	}

	/** square of this (more efficient than this*this) */
	ddouble sqr() const pure @safe nothrow
	{
		double p1, p2;
		twoSqr(hi, p1, p2);
		twoSumQuick(p1, p2 + 2*hi*lo, p1, p2);
		return ddouble(p1, p2);
	}

	/** square root of this */
	ddouble sqrt() const pure @safe nothrow
	{
		double x = 1/std.math.sqrt(hi);
		return this*x  + (this - (this*x).sqr) * x * 0.5; // TODO: optimize
	}


	//////////////////////////////////////////////////////////////////////
	/// comparison
	//////////////////////////////////////////////////////////////////////

	bool opEquals(double b) const pure @safe nothrow
	{
		return hi == b && lo == 0;
	}

	int opCmp(double b) const pure @safe nothrow
	{
		if(hi > b) return 1;
		if(hi < b) return -1;
		if(lo > 0) return 1;
		if(lo < 0) return -1;
		return 0;
	}

	int opCmp(ddouble b) const pure @safe nothrow
	{
		if(hi > b.hi) return 1;
		if(hi < b.hi) return -1;
		if(lo > b.lo) return 1;
		if(lo < b.lo) return -1;
		return 0;
	}

	bool opEquals(ddouble b) const pure @safe nothrow
	{
		return hi == b.hi && lo == b.lo;
	}
}

unittest
{
	assert(format("%.31e", ddouble(2).sqrt) == "1.414213562373095048801688724209");
	assert(format("%.29e", (ddouble(6)/ddouble(7))^^100) == "2.0198589217018753306533514440e-7");
	assert(format("%.31e", ddouble.pi.inverse) == "3.183098861837906715377675267450e-1");
}


//////////////////////////////////////////////////////////////////////
/// Basic operations with exact results. No rounding whatsoever.
//////////////////////////////////////////////////////////////////////

/**
 * NOTE: CTFE is explicitly disabled because during compile-time all doubles
 * are promoted to reals, which breaks pretty much all algorithms in this module.
 */

/** compute (hi+lo) = a + b */
private void twoSum(double a, double b, ref double hi, ref double lo) pure @safe nothrow
{
	if(__ctfe)
		assert(false, "compile-time double-double computations are not supported");

	hi = a+b;
	double v = hi-a;
	lo = (a-(hi-v))+(b-v);
}

/** compute (hi+lo) = a + b assuming abs(a) >= abs(b) */
private void twoSumQuick(double a, double b, ref double hi, ref double lo) pure @safe nothrow
{
	if(__ctfe)
		assert(false, "compile-time double-double computations are not supported");

	hi = a+b;
	lo = b-(hi-a);
}

/** compute (hi+lo) = a such that hi and lo only have ~27 bits of mantissa each */
private void split(double a, ref double hi, ref double lo) pure @safe nothrow
{
	if(__ctfe)
		assert(false, "compile-time double-double computations are not supported");

	// NOTE: this fails due to overflow for abs(a) > 2^996
	enum double s = 134217729.0; // 2^27+1
	double t = s*a;
	hi = t-(t-a);
	lo = a-hi;
	debug(ddouble) assert(lo + hi == a);
}

/** compute (hi+lo) = a * b */
private void twoProd(double a, double b, ref double hi, ref double lo) pure @safe nothrow
{
	if(__ctfe)
		assert(false, "compile-time double-double computations are not supported");

	// TODO: use hardware FMA when available
	double ahi, alo, bhi, blo;
	split(a, ahi, alo);
	split(b, bhi, blo);
	hi = a*b;
	lo = ((ahi*bhi - hi) + ahi*blo + alo*bhi) + alo*blo;
}

/** compute (hi+lo) = a * a */
private void twoSqr(double a, ref double hi, ref double lo) pure @safe nothrow
{
	if(__ctfe)
		assert(false, "compile-time double-double computations are not supported");

	// TODO: use hardware FMA when available
	double ahi, alo;
	split(a, ahi, alo);
	hi = a*a;
	lo = ((ahi*ahi - hi) + 2*ahi*alo) + alo*alo;
}
