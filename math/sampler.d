module math.sampler;

import std.algorithm;
import std.range;
import std.math;
import jive.array;
import mir.random;
import math.solve;
import math.integration;

/**
 * Adaptive Rejection Sampler for (nearly) arbitrary (log-) distributions.
 */
class LogSampler(F, FD, FDD)
{
	F f; // log-probability density
	FD fd; // first derivative of f
	FDD fdd; // second derivative of f

	int maxRegs; // maximum number of regions to create

	double a, b;
	Array!Region regs; // generally un-ordered

	Array!double sums;

	static struct Region
	{
		double a, b; // left/right bounds of region

		// upper/lower bound of f in the interval [a,b]
		// alpha + beta*x >= f(a+x(b-a)) >= gamma + delta*x
		double alpha, beta;
		double gamma, delta;

		/** area of upper bound */
		double area()
		{
			return (b-a)/beta*exp(alpha)*expm1(beta);
		}
	}

	private Region makeRegion(double a, double b)
	{
		assert(a < b);
		Region r;
		r.a = a;
		r.b = b;

		// secant through endpoints
		r.alpha = f(a);
		r.beta = f(b) - f(a);
		if(r.alpha == -double.infinity)
			r.beta = 0;

		// tangent to midpoint
		r.delta = fd(0.5*(a+b))*(b-a);
		r.gamma = f(0.5*(a+b)) - 0.5*r.delta;

		if(r.alpha < r.gamma)
		{
			swap(r.alpha, r.gamma);
			swap(r.beta, r.delta);
		}

		assert(r.alpha >= r.gamma);
		assert(r.alpha + r.beta >= r.gamma + r.delta);
		return r;
	}

	private void makePartialSums()
	{
		sums.resize(regs.length);
		double sum = 0;
		for(int i = 0; i < regs.length; ++i)
		{
			sum += regs[i].area;
			sums[i] = sum;
		}
	}

	this(F f, FD fd, FDD fdd, double a, double b, int maxRegs = 50)
	{
		assert(a < b);
		this.a = a;
		this.b = b;
		this.f = f;
		this.fd = fd;
		this.fdd = fdd;
		this.maxRegs = maxRegs;

		Array!double xs;
		Array!double fddxs;

		// sample f at 100 points
		xs.pushBack(a);
		for(int i = 1; i < 99; ++i)
			xs.pushBack(a+(b-a)*i/99);
		xs.pushBack(b);
		for(int i = 0; i < xs.length; ++i)
			fddxs.pushBack(fdd(xs[i]));

		// find inflection points to use as initial region bounds
		Array!double ys;
		ys.pushBack(xs[0]);
		for(int i = 1; i < xs.length; ++i)
			if(fddxs[i-1] * fddxs[i] < 0)
				ys.pushBack(solve!(double, x=>fdd(x))(xs[i-1], xs[i]));
		ys.pushBack(xs[$-1]);

		for(int i = 1; i < ys.length; ++i)
			regs.pushBack(makeRegion(ys[i-1], ys[i]));

		makePartialSums;
	}

	long nTries = 0;
	long nAccepts = 0;
	long nEvals = 0;

	double opCall(Rng)(ref Rng rng)
		if(isSaturatedRandomEngine!Rng)
	{
		while(true)
		{
			++nTries;

			// choose a region
			size_t i = sums[].assumeSorted.lowerBound(sums[$-1]*rand!double(rng).abs).length;

			double u = rand!double(rng).abs*(1 - exp(-regs[i].beta)) + exp(-regs[i].beta); // uniform in interval [exp(-beta), 1]
			double x = log(u)/regs[i].beta+1; // exponential in interval [0,1]
			if(regs[i].beta == 0)
				x = rand!double(rng).abs;
			assert(0 <= x && x <= 1);
			double y = rand!double(rng).abs*exp(regs[i].alpha + regs[i].beta*x);

			if(y > exp(regs[i].gamma + regs[i].delta*x))
			{
				++nEvals;
				if(y > exp(f(regs[i].a + (regs[i].b-regs[i].a)*x)))
				{
					// subdivide the interval
					if(regs.length < maxRegs)
					{
						double a = regs[i].a;
						double b = regs[i].b;
						double m = 0.5*(a+b);
						regs[i] = makeRegion(a,m);
						regs.pushBack(makeRegion(m,b));
						makePartialSums();
					}

					// try again
					continue;
				}
			}

			++nAccepts;
			return regs[i].a + (regs[i].b-regs[i].a)*x;
		}
	}

	double accProb() const @property
	{
		return cast(double)nAccepts/nTries;
	}

	double evalProb() const @property
	{
		return cast(double)nEvals/nTries;
	}

	/** plot f and its current approximation (for debugging) */
	void plot()() @property
	{
		import math.gnuplot;
		Array!double plotX, plotY, plotY2;
		sort!"a.a < b.a"(regs[]);

		foreach(const ref r; regs)
		{
			plotX.pushBack(r.a);
			plotX.pushBack(r.b);
			plotY.pushBack(r.gamma);
			plotY.pushBack(r.gamma + r.delta);
			plotY2.pushBack(r.alpha);
			plotY2.pushBack(r.alpha + r.beta);
		}

		auto plot = new Gnuplot;
		plot.plotFunction(f, a, b, 1000, "log-prob");
		plot.plotData(plotX[], plotY[], "min", "lines");
		plot.plotData(plotX[], plotY2[], "max", "lines");
	}
}
