module math.numerics;

/**
 * Some basic numerics helper functions.
 */

private import std.complex;
private import std.math;

/** Thrown when some method does not converge */
class NumericsException : Exception
{
    this(string s)
    {
        super(s);
    }
}

/** returns true if a = b with relative error <= eps */
bool approxEqual(T)(T a, T b, T eps = 4*RealTypeOf!T.epsilon) pure nothrow
{
    return abs(a-b) <= eps * max(abs(a), abs(b));
}

template isComplex(T)
{
	enum isComplex = is(T : Complex!A, A);
}

template RealTypeOf(T)
{
	static if(is(T : Complex!R, R))
		alias RealTypeOf = R;
	else
		alias RealTypeOf = T;
}

template ComplexTypeOf(T)
{
	static if(is(T : Complex!R, R))
		alias ComplexTypeOf = T;
	else
		alias ComplexTypeOf = Complex!T;
}

static assert(isComplex!(Complex!real));
static assert(!isComplex!float);
static assert(is(RealTypeOf!(Complex!float) == float));
static assert(is(RealTypeOf!float == float));
static assert(is(ComplexTypeOf!(Complex!float) == Complex!float));
static assert(is(ComplexTypeOf!float == Complex!float));

/** phase factor of x. +-1 for real x, e^(i arg(x)) for complex x */
auto phase(T)(auto ref const T x)
{
	// NOTE: result for x=0 can be arbitrary but should not be nan
	static if(isComplex!T)
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
    if(!isComplex!T)
{
    return x;
}

/**
 * Compute roots of x^2 + px + q. Returns the root with larger magnitude.
 * The other one should be obtained as q/root, as this is numerically more
 * stable than getting both by the direct formula.
 */
Complex!T polyRoot(T)(T p, T q)
	if(!isComplex!T)
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
	if(isComplex!T)
{
	auto d = sqrt(p*p/4 - q);
	if((p*std.complex.conj(d)).re < 0)
        return -p/2 + d;
    else
        return -p/2 - d;
}
