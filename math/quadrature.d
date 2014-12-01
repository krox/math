module math.quadrature;

private import std.functional;
private import std.exception;
private import math.orthpoly;

/** quadrature rules of the form ∫f ≈ ∑ w_i * f(x_i) */
class Quadrature
{
	// quadrature points and weights normalized for the interval [0,1]
	immutable(double)[] x;
	immutable(double)[] w;

	this(immutable(double)[] x, immutable(double)[] w)
	{
		assert(x.length == w.length);
		this.x = x;
		this.w = w;
	}

	/** integrate function fun from a to b */
	final double integrate(alias fun)(double a, double b)
	{
		double sum = 0;
		for(int i = 0; i < x.length; ++i)
			sum += w[i]*unaryFun!fun(a + (b-a) * x[i]);
		return sum*(b-a);
	}
}

/** only for debugging/demonstration */
class NaiveQuadrature : Quadrature
{
	/** create a quadrature rule using n points */
	this(int n)
	{
		assert(n >= 0);

		auto x = new double[n];
		auto w = new double[n];

		for(int i = 0; i < n; ++i)
		{
			x[i] = (i+0.5)/n;
			w[i] = 1.0/n;
		}

		super(assumeUnique(x), assumeUnique(w));
	}
}

/** quadrature using expansion in legendre polynomials */
class LegendreQuadrature : Quadrature
{
	/** create a quadrature rule with n points (exact up to order 2*n+1) */
	this(int n)
	{
		assert(n >= 0);

		auto x = cast(double[])legendreRoots!double(n);

		auto w = new double[n];
		for(int i = 0; i < n; ++i)
		{
			auto l = legendre!double(n, 1, x[i]);
			w[i] = 1/((1-x[i]*x[i])*l*l);
			x[i] = (x[i]+1)/2;
		}

		super(assumeUnique(x), assumeUnique(w));
	}
}
