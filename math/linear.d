module math.linear;

/**
 * linear algebra backend
 */

import jive.array;
private import std.math;
private import std.traits;
private import std.complex;
private import std.algorithm : move, swap;
private import jive.bitarray;
private import math.complex;


//////////////////////////////////////////////////////////////////////
/// LU decomposition
//////////////////////////////////////////////////////////////////////

/** compute LU decomposition of m inplace, put row permutation into p */
void denseComputeLU(T)(Slice2!T m, int[] p)
{
	int n = cast(int)p.length;
	if(m.size[0] != n || m.size[1] != n)
		throw new Exception("matrix dimension mismatch");

	for(int i = 0; i < n; ++i)
		p[i] = i;

	for(int k = 0; k < n; ++k)
	{
		// find a good pivot in column k
		static if(isFloatingPoint!T)
		{
			int pivot = k;
			for(int i = k+1; i < n; ++i)
				if(abs(m[i, k]) > abs(m[pivot, k]))
					pivot = i;
		}
		else
		{
			int pivot = k;
			for(int i = k; i < n; ++i)
				if(m[i, k] != 0)
				{
					pivot = i;
					break;
				}
		}
		if(m[pivot, k] == 0)
			throw new Exception("matrix not invertible");

		// swap the pivot row with row k
		swap(p[k], p[pivot]);
		for(int i = 0; i < n; ++i)
			swap(m[k, i], m[pivot,i]);

		// eliminate all entries below the pivot (which is now in m[k,k])
		for(int l = k+1; l < n; ++l)
		{
			m[l,k] /= m[k,k]; // this is now part of the L matrix
			for(int i = k+1; i < n; ++i)
				m[l,i] -= m[l,k] * m[k,i];
		}
	}
}

/** solve linear equation of lower triangular matrix with implicit 1's on diagonal */
void denseSolveL(T)(const Slice2!T m, Slice2!T b)
{
	int n = cast(int)m.size[0];
	if(m.size[1] != n || b.size[0] != n)
		throw new Exception("matrix dimension mismatch");

	for(int k = 0; k < n; ++k)
		for(int l = k+1; l < n; ++l)
			for(int i = 0; i < b.size[1]; ++i)
				b[l,i] -= m[l,k] * b[k,i];
}

/** solve linear equation of upper triangular matrix */
void denseSolveU(T)(Slice2!(const(T)) m, Slice2!T b)
{
	int n = cast(int)m.size[0];
	if(m.size[1] != n || b.size[0] != n)
		throw new Exception("matrix dimension mismatch");

	for(int k = n-1; k >= 0; --k)
	{
		for(int i = 0; i < b.size[1]; ++i)
			b[k,i] /= m[k,k];

		for(int l = k-1; l >= 0; --l)
			for(int i = 0; i < b.size[1]; ++i)
				b[l,i] -= m[l,k] * b[k,i];
	}
}

/** solve linear equation after LU decomposition was computed */
void denseSolveLU(T)(const Slice2!T m, Slice2!T b)
{
	int n = cast(int)m.size[0];
	if(m.size[1] != n || b.size[0] != n)
		throw new Exception("matrix dimension mismatch");

	denseSolveL!T(m, b);
	denseSolveU!T(m, b);
}


//////////////////////////////////////////////////////////////////////
/// QR decomposition (real numbers only implemented)
//////////////////////////////////////////////////////////////////////

/** same as scalarProduct(v, v) */
RealTypeOf!T norm2(T)(Slice!(const(T)) a)
{
	RealTypeOf!T r = 0;
	for(int i = 0; i < a.length; ++i)
		r += sqAbs(a[i]);
	return r;
}

/** scalar product of two vectors */
T scalarProduct(T, bool conjugate = true)(Slice!(const(T)) a, Slice!(const(T)) b)
{
	if(a.length != b.length)
		throw new Exception("vector dimension mismatch");

	T r = 0;
	for(int i = 0; i < a.length; ++i)
		static if(isComplex!T && conjugate)
			r += conj(a[i])*b[i];
		else
			r += a[i]*b[i];
	return r;
}

/** compute householder reflection inplace, returns beta factor */
RealTypeOf!T makeHouseholderReflection(T)(Slice!T v)
{
	T c = -phase(v[0])*sqrt(norm2!T(v));
	for(int i = 1; i < v.length; ++i)
		v[i] /= v[0] - c;
	v[0] = c;
	return 2/(1 + norm2(v[1..$]));
}

/** a = (1 - beta*|v><v|) * a */
void householderReflection(T)(Slice2!T a, Slice!(const(T)) v, RealTypeOf!T beta)
{
	int n = cast(int)v.size[0]+1;
	if(a.size[0] != n)
		throw new Exception("vector dimension mismatch");

	for(int k = 0; k < a.size[1]; ++k)
	{
		T s = beta * (a[0,k] + scalarProduct(v, a[1..$,k]));
		a[0, k] -= s;
		for(int i = 1; i < n; ++i)
			a[i, k] -= s * v[i-1];
	}
}

/** a = a * (1 - beta*|v><v|) */
void householderReflectionRight(T)(Slice2!T a, Slice!(const(T)) v, RealTypeOf!T beta)
{
	int n = cast(int)v.size[0]+1;
	if(a.size[1] != n)
		throw new Exception("vector dimension mismatch");

	for(int k = 0; k < a.size[0]; ++k)
	{
		T s = beta * (a[k,0] + scalarProduct!(T,false)(v, a[k,1..$]));
		a[k, 0] -= s;
		for(int i = 1; i < n; ++i)
		{
			static if(isComplex!T)
				a[k, i] -= s * conj(v[i-1]);
			else
				a[k, i] -= s * v[i-1];
		}
	}
}

/** compute QR decomposition of m inplace */
void denseComputeQR(T)(Slice2!T m, RealTypeOf!T[] beta)
{
	int n = cast(int)beta.length;
	if(m.size[0] != n || m.size[1] != n)
		throw new Exception("matrix dimension mismatch");

	for(int k = 0; k < n; ++k)
	{
		// make householder reflection for current column
		beta[k] = makeHouseholderReflection(m[k..$, k]);

		// update remaining columns
		householderReflection!T(m[k..$, k+1..$], m[k+1..$, k], beta[k]);
	}
}

/** solve linear equation after QR decomposition was computed */
void denseSolveQR(T)(Slice2!(const(T)) m, RealTypeOf!T[] beta, Slice2!T b)
{
	int n = cast(int)beta.length;
	if(m.size[0] != n || m.size[1] != n || b.size[0] != n)
		throw new Exception("matrix dimension mismatch");

	for(int k = 0; k < n; ++k)
		householderReflection!T(b[k..$,0..$], m[k+1..$,k], beta[k]);

	denseSolveU!T(m, b);
}

/** compute Hessenberg decomposition of m inplace */
void denseComputeHessenberg(T)(Slice2!T m, RealTypeOf!T[] beta)
{
	int n = cast(int)beta.length;
	if(m.size[0] != n || m.size[1] != n)
		throw new Exception("matrix dimension mismatch");

	for(int k = 0; k < n-1; ++k)
	{
		// make householder reflection for current column
		beta[k] = makeHouseholderReflection(m[k+1..$, k]);

		// update remaining columns
		householderReflection!T(m[k+1..$, k+1..$], m[k+2..$, k], beta[k]);

		// update all(!) rows (TODO: do a little less work on hermitian matrices)
		householderReflectionRight!T(m[0..$, k+1..$], m[k+2..$, k], beta[k]);
	}
}
