module math.numerics;

/**
 * Some basic numerics helper functions.
 */

private import std.complex;
private import std.math;
private import std.algorithm : min, max;
private import std.traits;


/** Thrown when some method does not converge */
class NumericsException : Exception
{
    this(string s = "some numerical method did not converge")
    {
        super(s);
    }
}

/** returns true if a = b with relative error <= eps */
bool approxEqual(T)(T a, T b, T eps = 4*RealTypeOf!T.epsilon) pure nothrow
{
    return abs(a-b) <= eps * max(abs(a), abs(b));
}

template RealTypeOf(T)
{
	static if(is(T : Complex!R, R))
		alias RealTypeOf = R;
	else
		alias RealTypeOf = T;
}

static assert(is(RealTypeOf!(Complex!float) == float));
static assert(is(RealTypeOf!float == float));

private ref int asInt(ref float x) pure nothrow
{
    return *cast(int*)&x;
}

private ref long asInt(ref double x) pure nothrow
{
    return *cast(long*)&x;
}

private int reverseNegative(int x) pure nothrow
{
    return x ^ ((x >> 31) & int.max);
}

private long reverseNegative(long x) pure nothrow
{
    return x ^ ((x >> 63) & long.max);
}

/**
 * Return a number that is in the middle of x and y with respect to their
 * floating point representation. Main use case for binary search.
 * TODO: probably replace by std.math.ieeeMean when available.
 */
T ieeeMean(T)(T x, T y) pure nothrow
{
    if(isNaN(x) || isNaN(y))
        return T.nan;

    // in the remaining cases (including de-normalized and infinite numbers)
    // the correct result will be obtained by re-interpreting floats as ints,
    // though the order of negative numbers is reversed.

    x.asInt = reverseNegative(x.asInt);
    y.asInt = reverseNegative(y.asInt);
    T z;
    z.asInt = (x.asInt >> 1) + (y.asInt >> 1) + (x.asInt & 1);
    z.asInt = reverseNegative(z.asInt);
    return z;
}

/** phase factor of x. +-1 for real x, e^(i arg(x)) for complex x */
auto phase(T)(auto ref const T x)
{
	// NOTE: result for x=0 can be arbitrary but should not be nan
	static if(is(T : Complex!R, R))
	{
		auto a = abs(x);
		if(a == 0)
			return T(1);
		else
			return x/a;
	}
	else
	{
		if(x >= 0)
			return 1;
		else
			return -1;
	}
}

T conj(T)(T x) pure nothrow
	if(!is(T : Complex!R, R))
{
	return x;
}

/**
 * Compute roots of x^2 + px + q. Returns the root with larger magnitude.
 * The other one should be obtained as q/root, as this is numerically more
 * stable than getting both by the direct formula.
 */
Complex!T polyRoot(T)(T p, T q)
	if(!is(T : Complex!R, R))
{
	auto d = p*p/4 - q;
	if(d < 0)
		return Complex!T(-p/2, sqrt(-d));

	if(p < 0)
		return Complex!T(-p/2+sqrt(d));
	else
		return Complex!T(-p/2-sqrt(d));
}

/** ditto */
T polyRoot(T)(T p, T q)
	if(is(T : Complex!R, R))
{
	auto d = sqrt(p*p/4 - q);
	if((p*std.complex.conj(d)).re < 0)
        return -p/2 + d;
    else
        return -p/2 - d;
}
