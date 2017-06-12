module math.integration;

/**
 * Various algorithms for numerical integration.
 * Currently only for 1-dimensional double -> double functions.
 */

private import std.stdio;
private import std.functional;
private import std.exception;
private import std.math;
private import math.mat;
private import math.statistics;
private import math.cuba;


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
Var integrateMC(alias _f, Dist)(Dist dist, long n)
{
    alias f = unaryFun!(_f, "x");

    double sum = 0;
    double sum2 = 0;

    for(long i = 0; i < n; ++i)
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

Var integrateUniformMC(alias _f)(double a, double b, long n)
{
	return integrateMC!_f(UniformDistribution(a, b), n);
}

Var integrateNormalMC(alias _f)(double mu, double sigma, long n)
{
	return integrateMC!_f(NormalDistribution(mu, sigma), n);
}

Var integrateExponentialMC(alias _f)(double lambda, long n)
{
	return integrateMC!_f(ExponentialDistribution(lambda), n);
}

Var integrateBoxMC(alias _f, size_t d)(long n)
{
	return integrateMC!(_f,BoxDistribution!d)(BoxDistribution!d.init, n);
}

Var integrateSimplexMC(alias _f, size_t d)(long n)
{
	return integrateMC!(_f,SimplexDistribution!d)(SimplexDistribution!d.init, n);
}

//////////////////////////////////////////////////////////////////////
/// monte-carlo integration with CUBA library
//////////////////////////////////////////////////////////////////////


private extern(C) int cubaIntegrand(size_t d)(const(int)* ndim, const(double)* xs, const(int)* ncomp, double* ys, void* userdata)
{
	assert(*ndim == d);
	assert(*ncomp == 1);

	Vec!(double, d) x = *cast(Vec!(double,d)*)xs;
	auto f = *cast(double delegate(Vec!(double,d))*)userdata;
	ys[0] = f(x);
	return 0;
}

private long dummy;

enum Cuba
{
	Vegas, Suave, Divonne, Cuhre
}

/**
 * Integrate f over [0,1]^d using the CUBA library.
 */
Var integrateCuba(size_t d)(double delegate(Vec!(double,d)) f, Cuba method, double epsabs = 1e-6, long maxeval = 10_000, ref long neval = dummy)
{
	// TODO: parameters need to be exposed for tweaking
	
	// general parameters
	int ndim = cast(int)d; // dimension of x
	int ncomp = 1; // dimension of y
	long nvec = 1; // vectorization
	int flags = 0; // 0-3 for verbosity
	int seed = 0; // 0 = sobol numbers
	long mineval = 100;
	double epsrel = 0;

	// vegas parameters
	long nstart = 1000; // evaluations per iteration
	long ninc = 500;
	long nbatch = 1000;
	int gridno = 0; // grid re-use

	// suave parameters
	long nnew = 1000;
	long nmin = 2;
	double flatness = 25.0;

	// divonne parameters
	int key1 = 47;
	int key2 = 1;
	int key3 = 1;
	int maxpass = 5;
	double border = 0;
	double maxchisq = 10;
	double mindeviation = 0.25;
	long ngiven = 0;
	int ldxgiven = ndim;
	double* xgiven = null;
	long nextra = 0;
	peakfinder_t peakfinder = null;

	// cuhre paramters
	int key = 0;

	// results
	int fail; // 0 = accurate result, >0 = inaccurate result, <0 = error
	double integral, error, prob;
	int nregions;

	switch(method)
	{
		case Cuba.Vegas:
			llVegas(ndim, ncomp, &cubaIntegrand!d, &f, nvec, epsrel, epsabs,
			flags, seed, mineval, maxeval,
			nstart, ninc, nbatch, gridno,
			null, null, &neval, &fail, &integral, &error, &prob);
			break;

		case Cuba.Suave:
			llSuave(ndim, ncomp, &cubaIntegrand!d, &f, nvec, epsrel, epsabs,
			flags, seed, mineval, maxeval,
			nnew, nmin, flatness,
			null, null, &nregions, &neval, &fail, &integral, &error, &prob);
			break;

		case Cuba.Divonne:
			llDivonne(ndim, ncomp, &cubaIntegrand!d, &f, nvec, epsrel, epsabs,
			flags, seed, mineval, maxeval,
			key1, key2, key3, maxpass, border, maxchisq, mindeviation,
			ngiven, ldxgiven, xgiven, nextra, peakfinder,
			null, null, &nregions, &neval, &fail, &integral, &error, &prob);
			break;

		case Cuba.Cuhre:
			llCuhre(ndim, ncomp, &cubaIntegrand!d, &f, nvec, epsrel, epsabs,
			flags, mineval, maxeval,
			key,
			null, null, &nregions, &neval, &fail, &integral, &error, &prob);
			break;

		default: assert(0);
	}

	if(fail < 0)
		throw new Exception("CUBA integration failed");
	if(prob > 0.95)
		stderr.writefln("WARNING: bad chi^2 probability in CUBA integration: %s", prob);

	return Var(integral, error^^2);
}
