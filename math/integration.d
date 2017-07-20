module math.integration;

/**
 * Numerical integration based on quadrature schemes such as Gauss-Legendre.
 */

private import std.exception;
private import std.math;
private import std.format;
private import math.solve;

//////////////////////////////////////////////////////////////////////
/// One-Dimensional integration
//////////////////////////////////////////////////////////////////////

/**
 * Gauss-Legendre quadrature using n function evaluations.
 */
double integrateGL(F)(F f, double a, double b, int n)
{
	double mid = (a+b)/2;
	double halfSize = (b-a)/2;

	auto q = legendreRoots!double(n);
	double sum = 0;
	if(n%2)
		sum += q.w[0]*f(mid);
	for(size_t i = n%2; i < q.x.length; ++i)
		sum += q.w[i] * (f(mid - halfSize * q.x[i]) + f(mid + halfSize * q.x[i]));
	return sum * halfSize;
}

//////////////////////////////////////////////////////////////////////
/// Multidimensionl integration
//////////////////////////////////////////////////////////////////////

/**
 * NOTE: these methods scale very badly with number of dimensions. In fact
 * for d larger than 2 or 3, chances are good that the randomized methods
 * in math.integration_mc are more suitable.
 */

/**
 * Gauss-Legendre quadrature using n^d function evaluations.
 */
double integrateGL(F)(F f, const(double)[] a, const(double)[] b, int n)
{
	size_t dim = a.length;
	assert(dim == b.length && 1 <= dim && dim <= 63);

	auto mid = new double[dim];
	mid[] = 0.5*(a[]+b[]);
	auto halfSize = new double[dim];
	halfSize[] = 0.5*(b[]-a[]);

	auto x = new double[dim];
	auto i = new size_t[dim];
	i[] = 0;

	auto q = legendreRoots!double(n);
	double sum = 0;
	while(i[$-1] < (n+1)/2)
	{
		// sample the function
		double fx = 0;
		outer: for(ulong signs = 0; signs < (1UL<<dim); ++signs)
		{
			for(size_t k = 0; k < dim; ++k)
			{
				if(signs & (1UL<<k))
					if(i[k] == 0 && n%2==1)
						continue outer;
				if(signs & (1UL<<k))
					x[k] = mid[k] + q.x[i[k]] * halfSize[k];
				else
					x[k] = mid[k] - q.x[i[k]] * halfSize[k];
			}
			fx += f(x[]);
		}

		// multiply with Legendre-weights
		for(size_t k = 0; k < dim; ++k)
			fx *= q.w[i[k]];
		sum += fx;

		// next point
		i[0] += 1;
		for(size_t k = 0; k < dim-1 && i[k] == (n+1)/2; ++k)
		{
			i[k] = 0;
			i[k+1] += 1;
		}
	}

	// normalize to integration volume
	for(size_t k = 0; k < dim; ++k)
		sum *= halfSize[k];
	return sum;
}

//////////////////////////////////////////////////////////////////////
/// Backend (orthogonal polynomials)
//////////////////////////////////////////////////////////////////////

/** Legendre polynomial */
struct Legendre
{
	int n;

	this(int n)
	{
		assert(n >= 0);
		this.n = n;
	}

	T opCall(T)(T x) pure
	{
		T a = T(0);
		T b = T(1);

		for(int k = 1; k <= n; ++k)
		{
			T c = ((2*k-1)*x*b - (k-1)*a) / k;
			a = b;
			b = c;
		}

		return b;
	}
}

/** derivative of Legendre polynomial */
struct LegendreD
{
	int n;

	this(int n)
	{
		assert(n >= 0);
		this.n = n;
	}

	T opCall(T)(T x) pure
	{
		if(n == 0)
			return T(0);

		T a = 0;
		T b = 1;

		for(int k = 2; k <= n; ++k)
		{
			T c = ((2*k-1)*x*b - k*a) / (k-1);
			a = b;
			b = c;
		}

		return b;
	}
}

struct LegendreQuadrature(T)
{
	immutable(T)[] x;
	immutable(T)[] w;

	void toString(scope void delegate(const(char)[]) sink, FormatSpec!char fmt) const
	{
		for(int i = 0; i < x.length; ++i)
		{
			formatValue(sink, x[i], fmt);
			sink(" ");
			formatValue(sink, w[i], fmt);
			sink("\n");
		}
	}

	string toString() const
	{
		return format("%.15e", this);
	}
}

/**
 * Compute Legendre roots and weights of order n. Because of symmetry, only
 * n/2 points need to be computed. Results are cached.
 */
LegendreQuadrature!T legendreRoots(T)(int n)
{
	assert(n >= 0);

	static LegendreQuadrature!T[int] cache;
	if(n in cache)
		return cache[n];

	auto x = new T[(n+1)/2];
	auto w = new T[(n+1)/2];

	for(int k = 0; k < (n+1)/2; ++k)
	{
		int i = (n-1)/2 - k;
		double a = PI*(4*i+3)/(4*n+2);
		if(n%2==1 && k == 0)
			x[k] = 0; // Newton would take rather long in this special case
		else
		{
			double guess = (1 - 1.0/(8.0*n*n))*cos(a);
			x[k] = solveNewton!T(Legendre(n), LegendreD(n), T(guess));
		}
		T d = LegendreD(n)(x[k]);
		w[k] = 2/((1-x[k]*x[k])*d*d);

		// check known bound
		assert(cos(PI*(4*i+2)/(4*n+2)) > x[k] && x[k] > cos(PI*(4*i+4)/(4*n+2)));
	}

	auto q = LegendreQuadrature!T(assumeUnique(x), assumeUnique(w));
	cache[n] = q;
	return q;
}

unittest
{
	immutable(double)[] x20 = [0.0765265211334973337546404,0.2277858511416450780804962,0.3737060887154195606725482,0.5108670019508270980043641,0.6360536807265150254528367,0.7463319064601507926143051,0.8391169718222188233945291,0.9122344282513259058677524,0.9639719272779137912676661,0.9931285991850949247861224];
	immutable(double)[] w20 = [0.1527533871307258506980843,0.1491729864726037467878287,0.1420961093183820513292983,0.1316886384491766268984945,0.1181945319615184173123774,0.1019301198172404350367501,0.0832767415767047487247581,0.0626720483341090635695065,0.0406014298003869413310400,0.0176140071391521183118620];

	for(size_t i = 0; i < x20.length; ++i)
	{
		assert(abs(x20[i] - cast(double)legendreRoots!double(20).x[i]) < 1e-15);
		assert(abs(w20[i] - cast(double)legendreRoots!double(20).w[i]) < 1e-15);
	}
}
