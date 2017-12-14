module math.sampler;

import std.algorithm;
import std.math;
import jive.array;
import mir.random;
import math.solve;
import math.integration;


/**
 * Random variable with arbitrary prabability density.
 */
class Sampler(F, FD)
{
	F f; // probability density
	FD fd;
	int n;
	double a, b; // boundaries
	Array!double xs;
	Array!double mins, maxs;

	this(F f, FD fd, double a, double b, int n, double eps = 1.0e-12, int maxIter = 100)
	{
		this.f = f;
		this.fd = fd;
		this.a = a;
		this.b = b;
		this.n = n;

		// initial equidistant subdivision
		xs.resize(n+1);
		xs[0] = a;
		xs[$-1] = b;
		for(int i = 1; i < n; ++i)
			xs[i] = a+(b-a)/n*i;

		// adapt the grid multiple times
		auto area = Array!double(n);
		auto ys = Array!double(n+1);
		ys[0] = a;
		ys[$-1] = b;
		for(int iter = 0; ; ++iter)
		{
			double norm = 0;
			assert(xs[0] == a && xs[$-1] == b);

			for(int i = 0; i < n; ++i)
			{
				assert(xs[i] < xs[i+1]);
				area[i] = integrateGK(f, xs[i], xs[i+1]);
				assert(area[i] > 0);
				norm += area[i];
			}

			double error = 0;
			for(int i = 0; i < n; ++i)
				error = max(error, abs(area[i]*n/norm-1));
			if(error < eps)
				break;
			if(iter > maxIter)
				throw new Exception("Sampler grid did not converge.");

			if(error > 0.01)
			{
				int j = 1;
				double sum = 0;
				for(int i = 0; i < n && j < n; ++i)
				{
					while(j*norm/n < sum+area[i] && j < n)
					{
						ys[j] = xs[i] + (xs[i+1]-xs[i])* (j*norm/n - sum) / area[i];
						ys[j] = max(ys[j], a);
						ys[j] = min(ys[j], b);
						++j;
					}
					sum += area[i];
				}
				assert(j==n);
			}
			else
			{
				double sum = 0;
				for(int i = 1; i < n; ++i)
				{
					sum += area[i-1];
					ys[i] = xs[i] + (norm/n*i-sum)/f(xs[i]);
				}
			}

			swap(xs,ys);
		}

		// determine minima and maxima in each interval
		mins.resize(n);
		maxs.resize(n);
		for(int i = 0; i < n; ++i)
		{
			// NOTE: a little leeway here does not impact correctness,
			// but only decreases performance slightly.
			mins[i] = minimize(f, fd, xs[i], xs[i+1])*0.99999;
			maxs[i] = maximize(f, fd, xs[i], xs[i+1])*1.00001;
		}
	}

	long nTries = 0;
	long nAccepts = 0;
	long nEvals = 0;

	double opCall(Rng)(ref Rng rng)
		if(isSaturatedRandomEngine!Rng)
	{
		int i = randIndex(rng, cast(uint)n);
		for(int hit = 0; hit < 50; ++hit)
		{
			++nTries;
			double x = xs[i] + (xs[i+1]-xs[i])*rand!double(rng).abs;
			double z = rand!double(rng).abs*maxs[i];
			if(z < mins[i])
			{
				++nAccepts;
				return x;
			}

			++nEvals;
			double fx = f(x);
			assert(mins[i] <= fx && fx <= maxs[i]);
			if(z < fx)
			{
				++nAccepts;
				return x;
			}
		}

		// good samplers will have acc-probs of 0.9 to 0.99, so 50 rejected
		// hits in a row are a very strong sign that something is wrong
		throw new Exception("Sampler could not generate a new random number.");
	}

	double accProb() const @property
	{
		return cast(double)nAccepts/nTries;
	}

	double evalProb() const @property
	{
		return cast(double)nEvals/nTries;
	}

	void plot()() const @property
	{
		import math.gnuplot;
		Array!double plotX, plotY, plotY2;
		for(int i = 0; i < n; ++i)
		{
			plotX.pushBack(xs[i]);
			plotX.pushBack(xs[i+1]);
			plotY.pushBack(mins[i]);
			plotY.pushBack(mins[i]);
			plotY2.pushBack(maxs[i]);
			plotY2.pushBack(maxs[i]);
		}

		auto plot = new Gnuplot;
		plot.plotFunction(f, a, b);
		plot.plotData(plotX[], plotY[], "min", "lines");
		plot.plotData(plotX[], plotY2[], "max", "lines");
	}
}

// minimizer for monotone/convex/concave functions
private double minimize(F, FD)(F f, FD fd, double a, double b, double sign = 1)
{
	assert(a < b);
	double r = min(sign*f(a), sign*f(b));
	if(fd(a)*fd(b) <= 0)
	{
		double m = solve!(double, x=>fd(x))(a,b);
		r = min(r, sign*f(m));
	}
	return r/sign;
}

private double maximize(F, FD)(F f, FD fd, double a, double b, double sign = 1)
{
	return minimize!F(f,fd,a,b,-sign);
}
