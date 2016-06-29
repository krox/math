module math.complex;

import std.complex;
import std.math;

private import std.algorithm;

template isComplex(T)
{
	enum isComplex = is(T : Complex!A, A);
}

static assert(isComplex!(Complex!float));
static assert(isComplex!(Complex!double));
static assert(isComplex!(Complex!real));
static assert(!isComplex!float);
static assert(!isComplex!double);
static assert(!isComplex!real);

template RealTypeOf(T)
{
	static if(is(T : Complex!R, R))
		alias RealTypeOf = R;
	else
		alias RealTypeOf = T;
}

static assert(is(RealTypeOf!(Complex!float) == float));
static assert(is(RealTypeOf!(Complex!double) == double));
static assert(is(RealTypeOf!(Complex!real) == real));
static assert(is(RealTypeOf!float == float));
static assert(is(RealTypeOf!double == double));
static assert(is(RealTypeOf!real == real));

auto phase(T)(auto ref const T x)
{
	static if(isComplex!T)
		if(x == 0)
			return T(1);
		else
			return x/abs(x);
	else
		if(x >= 0)
			return 1;
		else
		 return -1;
}

/** roots of quadratic polynomial */
T[2] polyRoot(T)(T c, T b, T a)
{
	b /= a;
	c /= a;

	T[2] x;
	x[0] = (-b+sqrt(b*b-4*c))/a;
	x[1] = (-b-sqrt(b*b-4*c))/a;
	if(abs(x[0]) > abs(x[1]))
		swap(x[0], x[1]);

	return x;
}
