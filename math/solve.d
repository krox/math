module math.solve;

/**
 * Numerically solve an equations f(x) = 0.
 * The template type of all functions is assumed to be float/double/real/Floating/...
 */

private import std.math : abs, signbit, isNaN, nextUp, nextDown;
private import std.functional : unaryFun;
private import std.algorithm;
private import math.numerics;

/**
 * General purpose method for solving f(x) = 0. Uses secant method and falls
 * back to bisecion if necessary. Result is accurate to full double precision.
 */
double solve(F)(F f, double a, double b)
{
	assert(!isNaN(a) && !isNaN(b));

	double fa = f(a);
	double fb = f(b);
	assert(!isNaN(fa) && !isNaN(fb));
	if(fa == 0)
		return a;
	if(fb == 0)
		return b;
	assert(signbit(fa) != signbit(fb));
	double c = b;
	double fc = fb;

	// a should be the best guess
	if(abs(fb) < abs(fa))
	{
		swap(a, b);
		swap(fa, fb);
	}

	for(int iter = 0; iter < 100; ++iter)
	{
		// choose new point x
		double x = (b*fa - a*fb) / (fa - fb);   // secant method
		if(!(a < x && x < c) && !(c < x && x < a)) // outside bracket (or nan) -> fall back to bisection
		{
			x = ieeeMean(a, c);
			if(x == a || x == c) // there is no further floating point number between a and b -> we are done
			{
				if(abs(fc) < abs(fa))
					return c;
				else
					return a;
			}
		}

		// evaluate f at new point
		b = a;
		fb = fa;
		a = x;
		fa = f(x);
		assert(!isNaN(fa));
		if(fa == 0)
			return a;

		// update brackets
		if(signbit(fa) != signbit(fb))
		{
			c = b;
			fc = fb;
		}
	}

	// TODO: this can be avoided by falling back to bisection not only when
	// when secant method leads outside of bracket, but also if previous step(s)
	// have been bad.
	throw new NumericsException;
}

unittest
{
	import std.math;
	// NOTE: these are accurate to full double precision
	assert(solve((double x)=>x*x-2, 0, 10) == cast(double)sqrt(2.0));
	assert(solve((double x)=>sin(x), 3, 4) == cast(double)PI);
}

/**
 * Newton-Raphson method. Unsafe without bracketing. Only use this if you know
 * that the initial guess is good and the function is nice enough.
 */
T solveNewton(T, F, FD)(F f, FD fd, T x)
{
	T fx = f(x);

	// NOTE: If everything is good the method converges quadratically so
	// 10 iterations should be more than sufficient for any practical type T
	for(int i = 0; i < 10; ++i)
	{
		// NOTE: either the sequence becomes consant, or it oscillates between
		// two values. In the latter case, choose the smaller abs(f(x)).
		T y = x - fx/fd(x);
		if(y == x)
			return x;
		T fy = f(y);
		if(abs(fy) >= abs(fx))
			return x;
		x = y;
		fx = fy;
	}
	throw new NumericsException;
}
