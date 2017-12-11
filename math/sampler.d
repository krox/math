module math.sampler;

import std.algorithm;
import std.math;
import jive.array;
import mir.random;
import math.integration;


/**
 * Random variable with arbitrary prabability density.
 */
class Sampler(F)
{
	F f; // probability density
	int n;
	double a, b; // boundaries
	Array!double xs;
	Array!double mins, maxs;

	this(F f, double a, double b, int n, double eps = 1.0e-12, int maxIter = 100)
	{
		this.f = f;
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
			mins[i] = minimize(f, xs[i], xs[i+1]);
			maxs[i] = maximize(f, xs[i], xs[i+1]);
		}
	}

	long nTries = 0;
	long nAccepts = 0;
	long nEvals = 0;

	double opCall(Rng)(ref Rng rng)
		if(isSaturatedRandomEngine!Rng)
	{
		int i = randIndex(rng, cast(uint)n);
		while(true)
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

// bisection-style minimizer that works for monotone and concave functions
private double minimize(F)(F f, double a, double b, double sign = 1)
{
	assert(a < b);
	const double eps = (b-a)*1.0e-10;

	double m = 0.5*(a+b);
	double fa = sign*f(a);
	double fb = sign*f(b);
	double fm = sign*f(m);

	while(b-a > eps)
	{
		double m2;
		if(m-a > b-m) // split left
			m2 = 0.5*(a+m);
		else // split right
			m2 = 0.5*(m+b);
		double fm2 = sign*f(m2);
		if(m2 < m)
		{
			swap(m, m2);
			swap(fm, fm2);
		}

		assert(a < m && m < m2 && m2 < b);
		if(fm < fm2)
		{
			b = m2;
			fb = fm2;
		}
		else
		{
			a = m;
			fa = fm;
			m = m2;
			fm = fm2;
		}
	}
	return fm/sign;
}

private double maximize(F)(F f, double a, double b, double sign = 1)
{
	return minimize!F(f,a,b,-sign);
}
