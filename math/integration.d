module math.integration;

/**
 * Various algorithms for numerical integration.
 * Currently only for 1-dimensional double -> double functions.
 */

private import std.stdio;
private import std.functional;
private import std.exception;
private import std.math;
private import std.random;
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
/// Monte-Carlo integration
//////////////////////////////////////////////////////////////////////

/**
 * NOTE: these methods integrate over the hypercube [0,1]^d. If you need a
 * different integration region, put the transformation inside the integrand.
 * TODO: possibly add some adapters to help with this.
 */

/** propotype of function to be integrated */
alias Integrand = double delegate(const(double)[] x);

/** some global configuration */
private enum long minEval = 100; // never stop with fewer evaluations

/**
 * Integrate f using Monte-Carlo with fixed distribution.
 */
Var integrateMC(Integrand f, size_t dim, double eps, long maxEval)
{
	assert(dim >= 1);
	assert(eps >= 0);
	assert(maxEval >= minEval);

	auto x = new double[dim];

	long nEval = 0;
	double mean = 0;
	double var = 0;

	while(nEval < minEval || (nEval < maxEval && var > eps*eps*nEval*(nEval-1)))
	{
		// evaluate the function at a random point
		foreach(ref xi; x)
			xi = uniform01();
		double fx = f(x);

		// update mean and variance
		++nEval;
		double delta = fx - mean;
		mean += delta/nEval;
		double delta2 = fx - mean;
		var += delta*delta2;
	}

	var /= nEval-1;

	return Var(mean, var/nEval);
}

private extern(C) int cubaIntegrand(const(int)* ndim, const(double)* xs, const(int)* ncomp, double* ys, void* userdata)
{
	assert(*ncomp == 1);

	auto x = xs[0..*ndim];
	auto f = *cast(Integrand*)userdata;
	ys[0] = f(x);
	return 0;
}

enum Cuba
{
	Vegas, Suave, Divonne, Cuhre
}

/**
 * Integrate f over [0,1]^d using the CUBA library.
 * NOTE: needs to be templated, so that CUBA don't need to be linked if not used
 */
Var integrateCuba(Cuba method)(Integrand f, int dim, double eps, long maxEval)
{
	// TODO: parameters need to be exposed for tweaking

	// general parameters
	int ncomp = 1; // dimension of y
	long nvec = 1; // vectorization
	int flags = 0; // 0-3 for verbosity
	int seed = 0; // 0 = sobol numbers
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
	int ldxgiven = cast(int)dim;
	double* xgiven = null;
	long nextra = 0;
	peakfinder_t peakfinder = null;

	// cuhre paramters
	int key = 0;

	// results
	int fail; // 0 = accurate result, >0 = inaccurate result, <0 = error
	double integral, error, prob;
	int nregions;
	long nEval;

	switch(method)
	{
		case Cuba.Vegas:
			llVegas(dim, ncomp, &cubaIntegrand, &f, nvec, epsrel, eps,
			flags, seed, minEval, maxEval,
			nstart, ninc, nbatch, gridno,
			null, null, &nEval, &fail, &integral, &error, &prob);
			break;

		case Cuba.Suave:
			llSuave(dim, ncomp, &cubaIntegrand, &f, nvec, epsrel, eps,
			flags, seed, minEval, maxEval,
			nnew, nmin, flatness,
			null, null, &nregions, &nEval, &fail, &integral, &error, &prob);
			break;

		case Cuba.Divonne:
			llDivonne(dim, ncomp, &cubaIntegrand, &f, nvec, epsrel, eps,
			flags, seed, minEval, maxEval,
			key1, key2, key3, maxpass, border, maxchisq, mindeviation,
			ngiven, ldxgiven, xgiven, nextra, peakfinder,
			null, null, &nregions, &nEval, &fail, &integral, &error, &prob);
			break;

		case Cuba.Cuhre:
			llCuhre(dim, ncomp, &cubaIntegrand, &f, nvec, epsrel, eps,
			flags, minEval, maxEval,
			key,
			null, null, &nregions, &nEval, &fail, &integral, &error, &prob);
			break;

		default: assert(0);
	}

	if(fail < 0)
		throw new Exception("CUBA integration failed");
	if(prob > 0.95)
		stderr.writefln("WARNING: bad chi^2 probability in CUBA integration: %s", prob);

	return Var(integral, error^^2);
}
