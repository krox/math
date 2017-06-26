module math.statistics;

private import std.math;
private import std.mathspecial;
private import std.random;
private import std.algorithm;
private import std.stdio;
private import std.format;
private import jive.array;
private import jive.internal;
private import math.mat;

struct Histogram
{
    double low = 0, high = 0;
    int nBins = 0;

    Array!long hist;
    long count = 0;
    long countLow = 0;
    long countHigh = 0;
    double sum = 0;
    double sum2 = 0;
    double min = double.infinity;
    double max = -double.infinity;

    /** create empty histogram */
    this(double low, double high, int nBins)
    {
        this.low = low;
        this.high = high;
        this.nBins = nBins;
        hist.resize(nBins, 0);
    }

    /** create histogram from data with automatic binning */
    this(const(double)[] xs)
    {
        double a = double.infinity;
        double b = -double.infinity;

        foreach(x; xs)
        {
            a = std.algorithm.min(a, x);
            b = std.algorithm.max(b, x);
        }

        a -= 0.0001*(b-a);
        b += 0.0001*(b-a);

        this(a, b, std.algorithm.min(std.algorithm.max(xs.length/10, 5), 50));

        foreach(x; xs)
            add(x);
    }

    void add(double x, long n = 1)
    {
        count += n;
        sum += n*x;
        sum2 += n*x*x;
        min = std.algorithm.min(min, x);
        max = std.algorithm.max(max, x);

        if(x < low)
            countLow += n;
        else if(x >= high)
            countHigh += n;
        else
            hist[cast(int)((x-low)/(high-low)*nBins)] += n;
    }

    double avg() const nothrow @property @safe
    {
        return sum / count;
    }

    double var() const nothrow @property @safe
    {
        return sum2/count - sum/count*sum/count;
    }

    /** iterate over center-of-bin / count-of-bin */
    int opApply(int delegate(double, double) dg) const
    {
        int r = 0;
        foreach(i, x; hist)
            if((r = dg(low + (i+0.5)*(high-low)/nBins, x)) != 0)
                break;
        return r;
    }

    void write() const
    {
        if(countLow)
            writefln("<:\t%s", countLow);
        foreach(x, y; this)
            writefln("%s:\t%s", x, y);
        if(countHigh)
            writefln(">:\t%s", countHigh);
        writefln("all:\t%s", count);
        writefln("avg = %s +- %s", avg, sqrt(var));
    }
}

/** put n-dimensional data points in, get average/variance/covariance out */
struct Statistics(size_t n = 1)
{
    double count = 0;
    double[n] sum;
    double[n][n] sum2;

    /** add a new data point */
    void add(Times!(n, double) xs, double w = 1)
    {
        // TODO: figure out a way to do this statically
        if(count == 0)
        {
            foreach(ref x; sum)
                x = 0;
            foreach(ref l; sum2)
                foreach(ref x; l)
                    x = 0;
        }

        count += w;
        foreach(i, x; xs)
            sum[i] += w*x;
        foreach(i, x; xs)
            foreach(j, y; xs)
                sum2[i][j] += w*x*y;
    }

    /** average in dimension i */
    double avg(size_t i = 0)() const nothrow @property @safe
    {
        return sum[i]/count;
    }

    /** variance in dimension i */
    double var(size_t i = 0)() const nothrow @property @safe
    {
        return cov!(i,i);
    }

    /** standard deviation in dimension i */
    double stddev(size_t i = 0)() const nothrow @property @safe
    {
        return sqrt(var!i);
    }

    /** covariance between diemnsions i and j */
    double cov(size_t i = 0, size_t j = 1)() const nothrow @property @safe
    {
        return sum2[i][j]/count - sum[i]/count * sum[j]/count;
    }

    /** correlation between diemnsions i and j */
    double corr(size_t i = 0, size_t j = 1)() const nothrow @property @safe
    {
        return cov!(i,j) / sqrt(var!i*var!j);
    }
}

/**
 * A (mean, variance) tuple. Can be interpreted as a random variable,
 * or a measurement with error. Standard arithmetic is overloaded to propagate
 * the error term, but there are severe limitations:
 *    -  but it assumes there is no correlation between variables,
 * so use carefully. For example "x+x" is not the same as "2*x" (and the latter
 * one is generally the correct one).
 */
struct Var
{
    double mean = double.nan;
    double var = 0;

	this(double mean, double var = 0) pure nothrow @safe
	{
		this.mean = mean;
		this.var = var;
	}

    /** standard deviation sqrt(variance) */
    double stddev() const pure nothrow @property @safe
    {
        return sqrt(var);
    }

    Var opBinary(string op)(double b) const pure nothrow @safe
    {
        switch(op)
        {
            case "+": return Var(mean + b, var);
            case "-": return Var(mean - b, var);
            case "*": return Var(mean*b, var*b*b);
            case "/": return Var(mean/b, var/(b*b));
            default: assert(false);
        }
    }

    Var opBinaryRight(string op)(double a) const pure nothrow @safe
    {
        switch(op)
        {
            case "+": return Var(a + mean, var);
            case "-": return Var(a - mean, var);
            case "*": return Var(a * mean, a*a*var);
            default: assert(false);
        }
    }

    Var opBinary(string op)(Var b) const pure nothrow @safe
    {
        switch(op)
        {
            case "+": return Var(mean + b.mean, var + b.var);
            case "-": return Var(mean - b.mean, var + b.var);
            case "*": return Var(mean * b.mean, mean*mean*b.var + b.mean*b.mean*var + var*b.var);
            default: assert(false);
        }
    }

    void opOpAssign(string op, S)(S b) pure nothrow @safe
    {
      this = this.opBinary!op(b);
    }

    /** returns human readable string "mean +- stddev" */
    string toString() const @property @safe
    {
        return format("%s +- %s", mean, stddev);
    }
}

/**
 * Estimate mean and variance of population as samples are coming in. This is
 * the same as the standard formula "Var(x) = n/(n-1) (E(x^2) - E(x)^2)" but
 * numerically more stable.
 */
struct Average
{
	long n = 0;
	double m = 0; // = 1/n ∑ x_i
	double m2 = 0; // = ∑ (x_i - m)^2

	/**
	 * Constructor taking data points. Trimming some fraction of the data
	 * makes it more robust, but potentially less accurate.
	 */
	this(const(double)[] xs, double trim = 0) pure
	{
		assert(0 <= trim && trim <= 1);
		Array!double arr;
		if(trim != 0)
		{
			arr = Array!double(xs[]);
			sort(arr[]);
			size_t cut = min(xs.length-1, cast(size_t)(trim*xs.length));
			xs = arr[cut/2 .. $-cut/2];
		}

		foreach(x; xs)
			add(x);
	}

	/** add a sample. NaN is silently ignored. */
	void add(double x) pure
	{
		if(isNaN(x))
			return;
		n += 1;
		double delta = x - m;
		m += delta/n;
		double delta2 = x - m;
		m2 += delta*delta2;
	}

	/** estimate of population mean */
	Var mean() const pure
	{
		if(n < 2)
			return Var(m, double.infinity);
		return Var(m, var/n);
	}

	/** estimate of population variance */
	double var() const pure
	{
		if(n < 2)
			return double.nan;
		return m2/(n-1);
	}

	/** reset everything */
	void clear() pure
	{
		n = 0;
		m = 0;
		m2 = 0;
	}
}


//////////////////////////////////////////////////////////////////////
/// Model fitting
//////////////////////////////////////////////////////////////////////

/** fit the constant function f() = a */
struct ConstantFit
{
	Var a;
	double prob;

	this(const(Var)[] ys, double exact = double.nan) pure
	{
		long ndf;

		if(isNaN(exact))
		{
			// compute weighted mean
			a = Var(0,0);
			foreach(y; ys)
			{
				a.mean += y.mean/y.var;
				a.var += 1/y.var;
			}
			a.mean /= a.var;
			a.var = 1/a.var;
			ndf = ys.length-1;
		}
		else
		{
			a = Var(exact, 0);
			ndf = ys.length;
		}

		// compute chi^2 test
		double chi2 = 0;
		foreach(y; ys)
			chi2 += (y.mean-a.mean)^^2/y.var;
		prob = chi2cdf(ndf, chi2);
	}

	string toString() const
	{
		return format("%s (χ²-prob = %.2f)", a, prob);
	}
}

/** fit a linear function f(x) = ax + b using the Theil-Sen estimator */
struct RobustLinearFit
{
	double a, b;

	this(const(double)[] xs, const(double)[] ys)
	{
		assert(xs.length == ys.length);

		Array!double arr;

		if(xs.length <= 20)
		{
			// gather all slopes in O(n^2)
			for(size_t i = 0; i < xs.length; ++i)
				for(size_t j = i+1; j < xs.length; ++j)
				{
					double slope = (ys[i]-ys[j])/(xs[i]-xs[j]);
					if(!isNaN(slope)) // ignore case of duplicate x-coords
						arr.pushBack(slope);
				}
		}
		else
		{
			// gather a random sample of slopes in "O(600)"
			for(int k = 0; k < 600; ++k)
			{
				size_t i = uniform(0, xs.length);
				size_t j = uniform(0, xs.length);
				double slope = (ys[i]-ys[j])/(xs[i]-xs[j]);
				if(!isNaN(slope)) // ignore case of duplicate x-coords
					arr.pushBack(slope);
			}
		}

		if(arr.empty)
		{
			a = b = double.nan;
			return;
		}

		// determine a = median of slopes
		sort(arr[]);
		a = arr[$/2];

		// determine b = median of (y - ax)
		arr.clear();
		for(size_t i = 0; i < xs.length; ++i)
			arr.pushBack(ys[i] - a*xs[i]);
		sort(arr[]);
		b = arr[$/2];
	}
}

/** cumulative distribution function of chi-square distribution */
double chi2cdf(long ndf, double chi2) pure
{
	if(ndf <= 0 || isNaN(chi2))
		return double.nan;
	if(chi2 <= 0)
		return 0;
	return gammaIncomplete(0.5*ndf, 0.5*chi2);
}

//////////////////////////////////////////////////////////////////////
/// probability distributions to be used in monte-carlo algorithms
//////////////////////////////////////////////////////////////////////

// TODO: restructure this into transformations, possibly without parameters
// TODO: decide on weight-function outside support (error/undefined/zero)

/** uniform distribution in the interval [a,b) */
struct UniformDistribution
{
	double a = 0;
	double b = 1;

	this(double a, double b)
	{
		assert(a < b);
		this.a = a;
		this.b = b;
	}

	double sample() const
	{
		return a + uniform01()*(b-a);
	}

	double weight(double x) const
	{
		if(x < a || x >= b)
			return 0;
		return 1/(b-a);
	}
}

/** normal distribution with mean mu and deviation sigma */
struct NormalDistribution
{
	double mu = 0;
	double sigma = 1;
	double c = 1/sqrt(2*PI);

	this(double mu, double sigma)
	{
		assert(sigma > 0);
		this.mu = mu;
		this.sigma = sigma;
		this.c = 1/sqrt(2*PI)/sigma;
	}

	double sample() const
	{
        double u = uniform01();
        double v = uniform01();

        //double x = std.mathspecial.normalDistributionInverse(u);  // correct, but slow (?)
        double x = sqrt(-2*log1p(-u)) * sin(2*PI*v);
        //double y = sqrt(-2*log1p(-u)) * cos(2*PI*v);

        return mu + sigma * x;
	}

	double weight(double x) const
	{
		return c*exp(-(x*x)/(sigma*sigma)/2);
	}
}

/** exponential distribution with parameter lambda */
struct ExponentialDistribution
{
	double lambda = 1;

	this(double lambda)
	{
		assert(lambda > 0);
		this.lambda = lambda;
	}

	double sample() const
	{
		return -log1p(-uniform01())/lambda;
	}

	double weight(double x) const
	{
		if(x < 0)
			return 0;
		return lambda*exp(-lambda*x);
	}
}

/** uniform distribution from the hypercube [0,1)^N */
struct BoxDistribution(size_t N)
{
	Vec!(double, N) sample() const
	{
		Vec!(double, N) r;
		for(int i = 0; i < N; ++i)
			r[i] = uniform01();
		return r;
	}

	double weight(Vec!(double, N) x) const
	{
		return 1;
	}
}

/** uniform sampling from the standard N-simplex */
struct SimplexDistribution(size_t N)
{
  private ExponentialDistribution exp;

  Vec!(double, N) sample() const
  {
    Vec!(double, N) r;
    double s = exp.sample;
    for(int i = 0; i < N; ++i)
    {
      r[i] = exp.sample();
      s += r[i];
    }
    for(int i = 0; i < N; ++i)
      r[i] /= s;

    return r;
  }

  double weight(Vec!(double, N) x) const
  {
    double r = 1;
    for(int i = 2; i <= N; ++i)
      r *= i;
    return r;
  }
}

/**
 * Generator for pseudo-random vectors in [0,1]^d distributed proportional to
 * some function g which does not need to be normalized. The algorithm is quite
 * efficient, but (1) initial values can be bad and (2) successive values are
 * correlated. You can use the "burn" and "step" parameters to mitigate these
 * problems to some degree.
 */
struct Metropolis
{
	private double[] x; // current point
	private double[] x2; // temporary
	private double delegate(const(double)[]) g;
	private double gx; // = g(x)
	private size_t step;

	this() @disable;

	this(double delegate(const(double)[]) g, size_t dim, size_t burn = -1, size_t step = 1)
	{
		if(burn == -1)
			burn = 2*dim;
		assert(dim >= 1);
		assert(burn >= 0);
		assert(step >= 1);

		this.g = g;
		this.step = step;
		x = new double[dim];
		x2 = new double[dim];

		foreach(ref xi; x)
			xi = uniform01();
		gx = g(x);

		for(int i = 0; i < burn; ++i)
			next();
	}

	private void next()
	{
		// new proposal point
		x2[] = x[];
		size_t k = uniform(0, x.length);
		x2[k] = uniform01();

		// transition with some amplitude
		double gx2 = g(x2);
		double p = abs(gx2/gx);
		if(p >= 1 || uniform01() < p)
		{
			swap(x, x2);
			swap(gx, gx2);
		}
	}

	double weight() const pure @property
	{
		return gx;
	}

	const(double)[] front() const pure @property
	{
		return x;
	}

	void popFront()
	{
		for(int i = 0; i < step; ++i)
			next();
	}
}
