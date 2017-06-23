module math.integration;

/**
 * Various algorithms for numerical integration.
 * Currently only for 1-dimensional double -> double functions.
 */

private import std.stdio;
private import std.functional;
private import std.algorithm;
private import std.exception;
private import std.math;
private import std.random;
private import jive.array;
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
	Average avg;

	while(nEval < minEval || (nEval < maxEval && avg.mean.var > eps*eps))
	{
		// evaluate the function at a random point
		foreach(ref xi; x)
			xi = uniform01();
		double fx = f(x);
		++nEval;

		// update mean and variance
		avg.add(fx);
	}

	return avg.mean;
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

struct MetropolisIntegration
{
	// the integral to do
	const Integrand f;
	const size_t dim;

	// parameters of method
	private double alpha = 1; // will be dynamically adjusted
	long burn = 0;
	long step = 1;
	double probArt = 0.5; // proposition-probability for real -> art
	double trimMean = 0.0;

	// result with error estimate
	Var estimate = Var(0, double.infinity);

	// current point (and a temporary proposal point)
	bool state;
	double[] x, x2;
	double fx;

	// statistics
	long nEval = 0; // # of evaluations of f
	long nReal = 0; // # of samples of real integral
	long nArt = 0; // # of samples of artificial integral
	long nAccept = 0; // # of accepted transitions
	long nReject = 0; // # of rejected transitions
	double fMin = double.infinity;
	double fMax = -double.infinity;

	/** constructor */
	this(Integrand f, size_t dim)
	{
		assert(f !is null);
		assert(dim >= 1);
		this.f = f;
		this.dim = dim;
		x = new double[dim];
		x2 = new double[dim];
	}

	private bool accept(double p)
	{
		if(p >= 1 || uniform01() < p)
		{
			++nAccept;
			return true;
		}
		else
		{
			++nReject;
			return false;
		}
	}

	private double eval(const(double)[] x)
	{
		double fx = f(x);
		++nEval;
		fMin = min(fMin, fx);
		fMax = max(fMax, fx);
		return fx;
	}

	/** Do one step of the Markov process */
	private void doStep()
	{
		if(state)
		{
			if(uniform01() < probArt) // real -> artificial
			{
				if(accept(abs(alpha/fx) / (probArt/1.0)))
					state = false;
			}
			else // real -> real (gibbs sampling)
			{
				x2[] = x[];
				x2[uniform(0, dim)] = uniform01();
				double fx2 = eval(x2);

				if(accept(abs(fx2/fx) / (0.5/0.5)))
				{
					swap(x, x2);
					swap(fx, fx2);
				}
			}
		}
		else // artificial -> real (new uniform point)
		{
			foreach(ref xi; x)
				xi = uniform01();
			fx = eval(x);
			if(accept(abs(fx/alpha) / (1.0/probArt)))
				state = true;
		}
	}

	/**
	 * Run n steps of the Markov process.
	 * NOTE: this can return +-infinity
	 */
	private double runBatch(long n)
	{
		state = false; // start in artificial integral

		long z0 = 0; // weights of the artificial...
		long z1 = 0; // ...and real integral

		for(long i = 0; i < burn; ++i)
			doStep();

		for(long k = 0; k < n; ++k)
		{
			for(long i = 0; i < step; ++i)
				doStep();

			if(state)
			{
				z1 += sgn(fx);
				++nReal;
			}
			else
			{
				z0 += 1;
				++nArt;
			}
		}

		if(z0 == 0)
		{
			writefln("WARNING: weight(0) = 0");;
			return z1*double.infinity;
		}

		return alpha*z1/z0;
	}

	/** run batches until accuary or max evaluations are reached */
	Var run(double eps, long maxEval)
	{
		assert(eps >= 0);
		assert(maxEval > 0);

		Array!double rs;
		Var estimate;

		long batchSize = 1000;

		nEval = 0;
		while(nEval < maxEval)
		{
			double r = runBatch(batchSize);
			if(isFinite(r))
				rs.pushBack(r);

			// estimate integral and error
			estimate = Average(rs[0..$], trimMean).mean;
			if(estimate.stddev <= eps)
				break;

			// adjust alpha to estimated integral
			if(isFinite(estimate.var))
				alpha = min(2.0*alpha, max(0.5*alpha, abs(estimate.mean)));
		}

		return estimate;
	}

	void printStats() const
	{
		writefln("fMin = %.2s (%.2s * alpha)", fMin, fMin/alpha);
		writefln("fMax = %.2s (%.2s * alpha)", fMax, fMax/alpha);
		writefln("art-rate = %.3f", cast(double)nArt/(nReal+nArt));
		writefln("acc-rate = %.3f", cast(double)nAccept/(nAccept+nReject));
	}
}

/**
 * Integrate f using the Metropolis-Hastings algorithm.
 */
Var integrateMH(Integrand f, size_t dim, double eps, long maxEval)
{
	auto state = MetropolisIntegration(f, dim);
	auto r = state.run(eps, maxEval);
	//state.printStats();
	return r;
}
