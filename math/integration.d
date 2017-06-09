module math.integration;

/**
 * Various algorithms for numerical integration.
 * Currently only for 1-dimensional double -> double functions.
 */

private import std.functional;
private import std.exception;
private import std.math;
private import math.statistics;


//////////////////////////////////////////////////////////////////////
/// direct quadrature formulas
//////////////////////////////////////////////////////////////////////

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

/+
/** quadrature using expansion in legendre polynomials */
class LegendreQuadrature : Quadrature
{
	/** create a quadrature rule with n points (exact up to order 2*n+1) */
	this(int n)
	{
		assert(n >= 0);

		auto x = legendreRoots!double(n).release;

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
+/


//////////////////////////////////////////////////////////////////////
/// adaptive integration
//////////////////////////////////////////////////////////////////////

/**
 * simple adaptive quadrature
 */
double integrate(alias fun)(double a, double b, double eps, int maxDepth = 15)
{
	alias f = unaryFun!fun;

	return integrateImpl!fun(a, b, f(a), f(b), eps, maxDepth);
}

double integrateImpl(alias fun)(double a, double b, double fa, double fb, double eps, int maxDepth)
{
	alias f = unaryFun!fun;

	double m = (a+b)/2;
	double fm = f(m);

	double i1 = (b-a) * fm;
	double i2 = (b-a)/2 * (fa + fb);
	double i3 = (b-a)/6 * (fa + 4*fm + fb);

	if(maxDepth == 0 || abs(i1-i2) <= eps)
		return i3;
	else
		return integrateImpl!f(a, m, fa, fm, eps/2, maxDepth-1) + integrateImpl!f(m, b, fm, fb, eps/2, maxDepth-1);
}


//////////////////////////////////////////////////////////////////////
/// monte-carlo integration
//////////////////////////////////////////////////////////////////////

/**
 * Monte Carlo integration with fixed distribution.
 */
Var integrateMC(alias _f, Dist)(Dist dist, int n)
{
    alias f = unaryFun!(_f, "x");

    double sum = 0;
    double sum2 = 0;

    for(int i = 0; i < n; ++i)
    {
        auto x = dist.sample();
        double fx = f(x)/dist.weight(x);

        sum += fx;
        sum2 += fx*fx;
    }

    sum /= n;
    sum2 /= n;

    return Var(sum, (sum2 - sum*sum)/(n-1));
}

Var integrateUniformMC(alias _f)(double a, double b, int n)
{
	return integrateMC!_f(UniformDistribution(a, b), n);
}

Var integrateNormalMC(alias _f)(double mu, double sigma, int n)
{
	return integrateMC!_f(NormalDistribution(mu, sigma), n);
}

Var integrateExponentialMC(alias _f)(double lambda, int n)
{
	return integrateMC!_f(ExponentialDistribution(lambda), n);
}

Var integrateSimplexMC(alias _f, size_t d)(int n)
{
	return integrateMC!(_f,SimplexDistribution!d)(SimplexDistribution!d.init, n);
}
