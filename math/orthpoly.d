module math.orthpoly;

/**
 * systems of orthogonal polynomials
 */

private import std.traits;
private import math.matrix;

import jive.array;

/** evaluate the m'th derivative of the n'th legendre polynomial at position x */
T legendre(T)(int n, int m, T x)
{
	assert(n >= 0 && m >= 0);

	if(n < m)
		return 0;

	T a = 0;
	T b = 1;
	T c;
	for(int k = m+1; k <= 2*m; ++k)
		b *= k/2.0;

	for(int k = m+1; k <= n; ++k)
	{
		c = ((2*k-1)*x*b - (k-1+m)*a) / (k-m);
		a = b;
		b = c;
	}

	return b;
}

/** returns the n roots of the n'th legendre polynomial */
Array!T legendreRoots(T)(int n)
	if(isFloatingPoint!T)
{
	assert(n >= 0);

	// special cases for very low order
	if(n == 0) return Array!T(0);
	if(n == 1) return Array!T([T(0)]);

	// compute roots as eigenvalues of a symmetric tridiagonal matrix
	auto m = Matrix!T.buildSymmetricBand!("i==j?0:j/sqrt(4*cast("~T.stringof~")j*j-1)")(n,1);
	return m.hermitianEigenvalues();
}
