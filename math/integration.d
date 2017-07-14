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
private import jive.priorityqueue;
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

/** flags for Monte-Carlo methods */
enum MC
{
	none = 0, // default
	stats = 1, // print some stats
}


/**
 * Integrate f over using uniform Monte-Carlo sampling.
 */
Var integrateMC(Integrand f, size_t dim, long n, MC flags = MC.none)
{
	auto a = new double[dim];
	auto b = new double[dim];
	a[] = 0;
	b[] = 1;
	auto x = new double[dim];
	return estimate(f, a, b, n, x);
}

struct MiserIntegration
{
	/** integral to be performed (regions defaults to [0,1]^d) */
	Integrand f;
	double[] left;
	double[] right;

	/** parameters of method (defaults set in constructor) */
	double estimateFrac; // fraction of evals used on variance-estimate
	long minEstimate; // minimum evals used for estimation
	long minBisect; // minimum evals to perform bisection
	long minEval; // minimum evals for any region total

	/** statistics */
	long nRegions = 0; // number of sub-regions the integration was split into
	long[] nSplit; // number of splits per dimension
	long nEval = 0; // total number of evals used
	long nRealEval = 0; // number of evals used for actual integration (without variance-estimation)

	/** temporary workspace */
	private double[] x;
	private double[] mid;
	private Average[] estLeft;
	private Average[] estRight;

	/** constructor */
	this(Integrand f, size_t dim)
	{
		assert(dim >= 1);
		this.f = f;

		// default parameters (partially taken from GSL)
		estimateFrac = 0.1;
		minEstimate = 16*dim;
		minBisect = 32*dim;
		minEval = 4*dim;

		// init region to [0,1]^d
		left = new double[dim];
		left[] = 0;
		right = new double[dim];
		right[] = 1;

		// allocate workspace
		x = new double[dim];
		mid = new double[dim];
		estLeft = new Average[dim];
		estRight = new Average[dim];
		nSplit = new long[dim];
	}

	/** recursively run MISER using a total of n function evaluations */
	Var run(long n)
	{
		assert(n >= minEval);

		// not enough points for bisection -> run a primitive MC integration
		if(n < minBisect)
		{
			nRegions += 1;
			nEval += n;
			nRealEval += n;

			double volume = 1;
			for(size_t i = 0; i < x.length; ++i)
				volume *= right[i] - left[i];

			Average avg;
			for(long k = 0; k < n; ++k)
			{
				for(size_t i = 0; i < x.length; ++i)
					x[i] = left[i] + uniform01() * (right[i] - left[i]);
				avg.add(f(x));
			}

			return avg.mean * volume;
		}

		// estimate variance on each half-space
		foreach(ref a; estLeft)
			a.clear;
		foreach(ref a; estRight)
			a.clear;
		for(size_t i = 0; i < x.length; ++i)
			mid[i] = left[i] + 0.5 * (right[i] - left[i]); // TODO: some "dither" for symmetry-breaking?

		long nEst = max(minEstimate, cast(long)(estimateFrac*n));
		n -= nEst;
		nEval += nEst;
		estimate(f, left, right, nEst, mid, estLeft, estRight, x);

		// choose "best" dimension
		double bestVar = double.infinity;
		size_t bestDim = -1;
		for(size_t i = 0; i < x.length; ++i)
		{
			if(estLeft[i].n < 2 || estRight[i].n < 2)
				throw new Exception("Need at least 2 points for variance-estimation.");

			double var = estLeft[i].var + estRight[i].var;
			if(var < bestVar)
			{
				bestVar = var;
				bestDim = i;
			}
		}
		assert(0 <= bestVar && bestVar < double.infinity);
		++nSplit[bestDim];

		// determine splitting of points to the half-spaces
		double sLeft = sqrt(estLeft[bestDim].var);
		double sRight = sqrt(estRight[bestDim].var);
		double leftFrac = sLeft / (sLeft + sRight);
		long nLeft = cast(long)(leftFrac * n);
		nLeft = max(nLeft, minEval);
		nLeft = min(nLeft, n-minEval);

		// run recursively on both half-spaces
		// NOTE: need to restore left/right afterwards
		double l = left[bestDim];
		double r = right[bestDim];
		double m = mid[bestDim];
		left[bestDim] = l;
		right[bestDim] = m;
		Var intLeft = run(nLeft);
		left[bestDim] = m;
		right[bestDim] = r;
		Var intRight = run(n-nLeft);
		left[bestDim] = l;
		right[bestDim] = r;

		return intLeft + intRight;
	}

	void printStats()
	{
		writefln("# of regions = %s", nRegions);
		writefln("# of splits = %s", nSplit[]);
		writefln("# of evals = %s (%.2f real)", nEval, cast(double)nRealEval/nEval);
	}
}

/**
 * Integrate f over [0,1]^d using the MISER algorithm.
 * NOTE: eps != 0 is not supported
 */
Var integrateMiser(Integrand f, size_t dim, long n, MC flags = MC.none)
{
	auto state = MiserIntegration(f, dim);
	auto r = state.run(n);
	if(flags & MC.stats)
		state.printStats();
	return r;
}

/**
 * Adaptive Monte-Carlo integration.
 * TODO: find the name of this method in literature
 * TODO: improve memory consumption (at least do fewer/bigger allocs)
 * TODO: fix the error estimation
 */
Var integrateOrange(Integrand f, size_t dim, long n, MC flags = MC.none)
{
	assert(dim >= 1);

	long batchSize = 32*dim; // TODO: tweak this

	PriorityQueue!Region regs;
	auto zero = new double[dim];
	auto one = new double[dim];
	zero[] = 0;
	one[] = 1;
	regs.pushBack(Region(zero, one, Var(0, double.infinity)));

	auto estLeft = new Average[dim];
	auto estRight = new Average[dim];
	auto mid = new double[dim];
	auto x = new double[dim];

	for(long nEval = 0; nEval < n; )
	{
		auto r = regs.pop();
		double vol = 1;
		for(int i = 0; i < dim; ++i)
		{
			estLeft[i].clear();
			estRight[i].clear();
			mid[i] = 0.5*(r.a[i] + r.b[i]);
			vol *= r.b[i] - r.a[i];
		}

		estimate(f, r.a, r.b, batchSize, mid, estLeft, estRight, x);
		nEval += batchSize;

		// choose "best" dimension
		double bestVar = double.infinity;
		size_t bestDim = -1;
		double worstVar = 0;
		size_t worstDim = -1;
		for(size_t i = 0; i < x.length; ++i)
		{
			if(estLeft[i].n < 2 || estRight[i].n < 2)
				throw new Exception("Need at least 2 points for variance-estimation.");

			double var = estLeft[i].var + estRight[i].var;
			if(var < bestVar)
			{
				bestVar = var;
				bestDim = i;
			}
			if(var > worstVar)
			{
				worstVar = var;
				worstDim = i;
			}
		}
		assert(0 <= bestVar && bestVar < double.infinity);

		auto m1 = new double[dim];
		auto m2 = new double[dim];
		m1[] = r.b[];
		m1[bestDim] = mid[bestDim];

		m2[] = r.a[];
		m2[bestDim] = mid[bestDim];

		regs.pushBack(Region(r.a, m1, 0.5*vol*estLeft[bestDim].mean));
		regs.pushBack(Region(m2, r.b, 0.5*vol*estRight[bestDim].mean));
	}

	auto est = Var(0,0);
	auto rs = regs.release;
	foreach(r; rs[])
		est += r.estimate;
	return est;
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
Var integrateCuba(Cuba method)(Integrand f, int dim, long maxEval, MC fl = MC.none)
{
	// TODO: parameters need to be exposed for tweaking

	// general parameters
	int ncomp = 1; // dimension of y
	long nvec = 1; // vectorization
	int flags = (fl & MC.stats) ? 1 : 0; // 0-3 for verbosity
	int seed = 0; // 0 = sobol numbers
	double eps = 0;
	double epsrel = 0;
	int minEval = 100;

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

struct MetropolisIntegration(bool simplex = false)
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

	private static newPoint(double[] x)
	{
		static if(simplex)
		{
			double s = 0;
			for(size_t i = 0; i < x.length; ++i)
			{
				x[i] = log(uniform01());
				s += x[i];
			}
			x[] /= s;
		}
		else
		{
			for(size_t i = 0; i < x.length; ++i)
				x[i] = uniform01();
		}
	}

	private static changePoint(double[] x)
	{
		static if(simplex)
		{
			size_t i = uniform(0, x.length);
			size_t j = (i+uniform(1, x.length)) % x.length;
			double z = uniform01();
			double s = x[i]+x[j];
			x[i] = z*s;
			x[j] = (1-z)*s;
		}
		else
		{
			size_t i = uniform(0,x.length);
			x[i] = uniform01();
		}
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
				changePoint(x2);
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
			newPoint(x);
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
				z1 += cast(int)sgn(fx);
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
			writefln("WARNING: weight(0) = 0");
			return z1*double.infinity;
		}

		return alpha*z1/z0;
	}

	/** run batches using ~n evaluations */
	Var run(long n)
	{
		assert(n > 0);

		Array!double rs;
		Var estimate;

		long batchSize = 1000;

		nEval = 0;
		while(nEval < n)
		{
			double r = runBatch(batchSize);
			if(isFinite(r))
				rs.pushBack(r);

			// estimate integral and error
			estimate = Average(rs[0..$], trimMean).mean;

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
Var integrateMH(bool simplex = false)(Integrand f, size_t dim, long n, MC flags = MC.none)
{
	auto state = MetropolisIntegration!simplex(f, dim);
	auto r = state.run(n);
	if(flags & MC.stats)
		state.printStats();
	return r;
}


//////////////////////////////////////////////////////////////////////
/// internal helper functions
//////////////////////////////////////////////////////////////////////

private:

/**
 * Most basic Monte-Carlo integration: Uniformly sample n point in [a,b] and
 * compute average with error. x is used as temporary workspace.
 */
Var estimate(Integrand f, const(double)[] a, const(double)[] b, long n, double[] x)
{
	double vol = 1;
	for(size_t i = 0; i < x.length; ++i)
		vol *= b[i] - a[i];

	Average avg;
	for(long k = 0; k < n; ++k)
	{
		for(size_t i = 0; i < x.length; ++i)
			x[i] = uniform(a[i], b[i]);
		avg.add(f(x));
	}

	return vol * avg.mean;
}

/** Sample f at n points in [a,b] and track average for all half-spaces. */
void estimate(Integrand f, const(double)[] a, const(double)[] b, long n, const(double)[] mid, Average[] estLeft, Average[] estRight, double[] x)
{
	for(long k = 0; k < n; ++k)
	{
		// generate a point
		for(size_t i = 0; i < x.length; ++i)
		{
			if(i == (k/2) % x.length) // stratify one dimension
			{
				if(k % 2 == 0)
					x[i] = a[i] + uniform01()*(mid[i] - a[i]); // force left
				else
					x[i] = mid[i] + uniform01()*(b[i] - mid[i]); // force right
			}
			else // completely uniform in other dimensions
				x[i] = a[i] + uniform01()*(b[i] - a[i]);
		}

		// sample the point
		double fx = f(x);

		// take statistics
		for(size_t i = 0; i < x.length; ++i)
			if(x[i] < mid[i])
				estLeft[i].add(fx);
			else
				estRight[i].add(fx);
	}
}

struct Region
{
	const(double)[] a, b;
	Var estimate = Var(0, double.infinity);

	int opCmp(Region other) const pure
	{
		if(estimate.var < other.estimate.var)
			return 1;
		if(estimate.var > other.estimate.var)
			return -1;
		return 0;
	}
}
