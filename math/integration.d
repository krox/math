module math.integration;

/**
 * Numerical integration based on quadrature schemes such as Gauss-Legendre.
 */

import std.exception;
import std.math;
import std.format;
import jive.priorityqueue;
import math.solve;

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

/**
 * Adaptive Gauss-Kronrod qudrature using 15/31 points.
 */
double integrateGK(F)(F f, double a, double b)
{
	static struct Region
	{
		double a, b;
		double val, err;
		this(F f, double a, double b)
		{
			auto est = integrateKronrod(f, a, b, GK15);
			this.a = a;
			this.b = b;
			this.val = est[1];
			this.err = abs(est[1]-est[0]);
		}
	}

	PriorityQueue!(Region, "a.err > b.err") q;

	auto reg = Region(f, a, b);
	q.pushBack(reg);
	double val = reg.val;
	double err = reg.err;
	while(abs(err/val) > 1.0e-10)
	{
		if(q.length >= 100)
			throw new Exception("Gauss-Kronrod adaptive integral did not converge.");

		reg = q.pop;
		auto regLeft = Region(f, reg.a, 0.5*(reg.a+reg.b));
		auto regRight = Region(f, 0.5*(reg.a+reg.b), reg.b);
		val += regLeft.val + regRight.val - reg.val;
		err += regLeft.err + regRight.err - reg.err;
		q.push(regLeft);
		q.push(regRight);
	}
	return val;
}

private double[2] integrateKronrod(F)(F f, double a, double b, GaussKronrodQuadrature!double q)
{
	double mid = (a+b)/2;
	double half = (b-a)/2;

	assert(q.x[0] == 0); // only odd GK rules are supported
	double f0 = f(mid);
	double sumG = q.wG[0]*f0;
	double sumK = q.wK[0]*f0;
	for(size_t i = 1; i < q.wG.length; ++i)
	{
		f0 = f(mid - half * q.x[i]) + f(mid + half * q.x[i]);
		sumG += q.wG[i] * f0;
		sumK += q.wK[i] * f0;
	}
	for(size_t i = q.wG.length; i < q.wK.length; ++i)
	{
		f0 = f(mid - half * q.x[i]) + f(mid + half * q.x[i]);
		sumK += q.wK[i] * f0;
	}

	sumG *= half;
	sumK *= half;
	return [sumG, sumK];
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

struct GaussKronrodQuadrature(T)
{
	// wG.length < wK.length = x.length
	immutable(T)[] x;
	immutable(T)[] wG;
	immutable(T)[] wK;
}

// values taken from https://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/
private static GK15 = GaussKronrodQuadrature!double(
	[ // Gauss nodes
	0.000000000000000000000000000000000e+00,
	2.011940939974345223006283033945962e-01,
	3.941513470775633698972073709810455e-01,
	5.709721726085388475372267372539106e-01,
	7.244177313601700474161860546139380e-01,
	8.482065834104272162006483207742169e-01,
	9.372733924007059043077589477102095e-01,
	9.879925180204854284895657185866126e-01,
	// additional Kronrod nodes
	1.011420669187174990270742314473923e-01,
	2.991800071531688121667800242663890e-01,
	4.850818636402396806936557402323506e-01,
	6.509967412974169705337358953132747e-01,
	7.904185014424659329676492948179473e-01,
	8.972645323440819008825096564544959e-01,
	9.677390756791391342573479787843372e-01,
	9.980022986933970602851728401522712e-01,
	],[ // Gauss weights
	2.025782419255612728806201999675193e-01,
	1.984314853271115764561183264438393e-01,
	1.861610000155622110268005618664228e-01,
	1.662692058169939335532008604812088e-01,
	1.395706779261543144478047945110283e-01,
	1.071592204671719350118695466858693e-01,
	7.036604748810812470926741645066734e-02,
	3.075324199611726835462839357720442e-02,
	],[ // Kronrod weights
	1.013300070147915490173747927674925e-01,
	9.917359872179195933239317348460313e-02,
	9.312659817082532122548687274734572e-02,
	8.308050282313302103828924728610379e-02,
	6.985412131872825870952007709914748e-02,
	5.348152469092808726534314723943030e-02,
	3.534636079137584622203794847836005e-02,
	1.500794732931612253837476307580727e-02,

	1.007698455238755950449466626175697e-01,
	9.664272698362367850517990762758934e-02,
	8.856444305621177064727544369377430e-02,
	7.684968075772037889443277748265901e-02,
	6.200956780067064028513923096080293e-02,
	4.458975132476487660822729937327969e-02,
	2.546084732671532018687400101965336e-02,
	5.377479872923348987792051430127650e-03,
	]);
