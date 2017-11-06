module math.linear;

/**
 * linear algebra backend
 */

import std.math;
import std.traits;
import std.complex;
import std.algorithm : swap, min, max;
import mir.ndslice;
import jive.bitarray;
import math.numerics;

/** Matrix decompositions

name				decomposition 		notes

LU 			PA = LU				L has 1's on the diagonal
QR			A = QR
Hessenberg	A = QHQ^-1
Schur		A = QUQ^-1			does not always exist for real matrices
Cholesky	A = LDL^T			only for Hermitian A

*/

//////////////////////////////////////////////////////////////////////
/// LU and LDL decomposition
//////////////////////////////////////////////////////////////////////

/** compute LU decomposition of m inplace, put row permutation into p */
void denseComputeLU(T)(ContiguousMatrix!T m, int[] p)
{
	int n = cast(int)p.length;
	if(m.length!0 != n || m.length!1 != n)
		throw new Exception("matrix dimension mismatch");

	for(int i = 0; i < n; ++i)
		p[i] = i;

	for(int k = 0; k < n; ++k)
	{
		// find a good pivot in column k
		static if(isFloatingPoint!T || is(T : Complex!R, R))
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
		each!swap(m[k], m[pivot]);

		// eliminate all entries below the pivot (which is now in m[k,k])
		for(int l = k+1; l < n; ++l)
		{
			m[l, k] /= m[k, k]; // this is now part of the L matrix
			m[l, k+1 .. n] -= m[l,k] * m[k, k+1 .. n];
		}
	}
}

/** solve linear equation of lower triangular */
void denseSolveL(T, bool implicitOne, M, B)(M m, B b)
	if(isMatrix!M && isMatrix!B)
{
	int n = cast(int)m.length!0;
	if(m.length!1 != n || b.length!0 != n)
		throw new Exception("matrix dimension mismatch");

	for(int k = 0; k < n; ++k)
	{
		static if(!implicitOne)
			b[k, 0..$] /= m[k,k];
		for(int l = k+1; l < n; ++l)
			b[l, 0..$] -= m[l,k] * b[k, 0..$];
	}
}

/** solve linear equation of upper triangular matrix */
void denseSolveU(T, bool implicitOne, M, B)(M m, B b)
	if(isMatrix!M && isMatrix!B)
{
	denseSolveL!(T, implicitOne)(m.reversed!(0,1), b.reversed!0);
}

/** solve linear equation after LU decomposition was computed */
void denseSolveLU(T)(ContiguousMatrix!(const(T)) m, ContiguousMatrix!T b)
{
	int n = cast(int)m.length!0;
	if(m.length!1 != n || b.length!0 != n)
		throw new Exception("matrix dimension mismatch");

	denseSolveL!(T,true)(m, b);
	denseSolveU!(T,false)(m, b);
}

/** compute LDL^T decomposition inplace. Only uses lower triangle including diagonal */
void denseComputeLDL(T)(ContiguousMatrix!T m)
{
	int n = cast(int)m.length!0;
	if(m.length!1 != n)
		throw new Exception("matrix dimension mismatch");

	for(int j = 0; j < n; ++j)
	{
		for(int k = 0; k < j; ++k)
			m[j,j] -= m[j,k].sqAbs * m[k,k];
		for(int i = j+1; i < n; ++i)
		{
			for(int k = 0; k < j; ++k)
				m[i,j] -= m[i,k]*m[j,k].conj*m[k,k];
			m[i,j] /= m[j,j];
		}
	}
}

/** solve linear equation after LDL decomposition was computed */
void denseSolveLDL(T)(ContiguousMatrix!(const(T)) m, ContiguousMatrix!T b)
{
	int n = cast(int)m.length!0;
	if(m.length!1 != n || b.length!0 != n)
		throw new Exception("matrix dimension mismatch");

	denseSolveL!(T,true)(m, b);
	for(int i = 0; i < n; ++i)
		b[i, 0..$] /= m[i,i];
	denseSolveU!(T,true)(m.adjoint, b);
}


//////////////////////////////////////////////////////////////////////
/// QR and Hessenberg decomposition
//////////////////////////////////////////////////////////////////////

/** compute householder reflection inplace, returns beta factor */
RealTypeOf!T makeHouseholder(T, V)(V v)
	if(isVector!V)
{
	T c = -phase(v[0])*sqrt(sqNorm2(v));

	if(c == 0) // typically this indicate a singular matrix
		return RealTypeOf!T(0);

	v[1..$] /= v[0] - c;
	v[0] = c;
	return 2/(1 + sqNorm2(v[1..$]));
}

/** a = (1 - beta*|v><v|) * a */
void applyHouseholder(T, A, V)(A a, V v, RealTypeOf!T beta)
	if(isMatrix!A && isVector!V)
{
	int n = cast(int)v.length!0+1;
	if(a.length!0 != n)
		throw new Exception("vector dimension mismatch");

	for(int k = 0; k < a.length!1; ++k)
	{
		T s = beta * (a[0,k] + dot(v.conj, a[1..$,k]));
		a[0, k] -= s;
		a[1..n, k] -= s * v[0..n-1];
	}
}

/** a = a * (1 - beta*|v><v|) */
void applyHouseholderRight(T, A, V)(A a, V v, RealTypeOf!T beta)
	if(isMatrix!A && isVector!V)
{
	applyHouseholder!T(a.transposed, v.conj, beta);
}

/** compute givens rotation */
void makeGivens(T)(ref T a, ref T b, ref RealTypeOf!T c, ref T s)
{
	if(b == 0)
	{
		c = 1;
		s = 0;
		//a = a;
		//b = 0;
	}
	else if(a == 0)
	{
		c = 0;
		s = phase(conj(b));
		a = abs(b);
		b = 0;
	}
	else
	{
		auto l = sqrt(sqAbs(a)+sqAbs(b));
		c = abs(a)/l;
		s = phase(a)*conj(b)/l;
		a = phase(a)*l;
		b = 0;
	}
}

void applyGivens(T)(UniversalMatrix!T m, int i, int j, RealTypeOf!T c, ref T s)
{
	for(int k = 0; k < m.length!1; ++k)
	{
		auto tmp = c*m[i,k] + s*m[j,k];
		m[j,k] = -conj(s)*m[i,k] + c*m[j,k];
		m[i,k] = tmp;
	}
}

void applyGivensRight(T)(UniversalMatrix!T m, int i, int j, RealTypeOf!T c, ref T s)
{
	for(int k = 0; k < m.length!1; ++k)
	{
		auto tmp = c*m[k,i] + conj(s)*m[k,j];
		m[k,j] = -s*m[k,i] + c*m[k,j];
		m[k,i] = tmp;
	}
}

/** compute QR decomposition of m inplace */
void denseComputeQR(T)(ContiguousMatrix!T m, RealTypeOf!T[] beta)
{
	int n = cast(int)beta.length;
	if(m.length!0 != n || m.length!1 != n)
		throw new Exception("matrix dimension mismatch");

	for(int k = 0; k < n; ++k)
	{
		// make householder reflection for current column
		beta[k] = makeHouseholder!T(m[k..$, k]);

		// update remaining columns
		applyHouseholder!T(m[k..$, k+1..$], m[k+1..$, k], beta[k]);
	}
}

/** solve linear equation after QR decomposition was computed */
void denseSolveQR(T)(ContiguousMatrix!(const(T)) m, const(RealTypeOf!T)[] beta, ContiguousMatrix!T b)
{
	int n = cast(int)beta.length;
	if(m.length!0 != n || m.length!1 != n || b.length!0 != n)
		throw new Exception("matrix dimension mismatch");

	for(int k = 0; k < n; ++k)
		applyHouseholder!T(b[k..$,0..$], m[k+1..$,k], beta[k]);

	denseSolveU!(T,false)(m, b);
}

/** compute Hessenberg decomposition of m inplace */
void denseComputeHessenberg(T)(ContiguousMatrix!T m, RealTypeOf!T[] beta)
{
	int n = cast(int)beta.length;
	if(m.length!0 != n || m.length!1 != n)
		throw new Exception("matrix dimension mismatch");

	for(int k = 0; k < n-1; ++k)
	{
		// make householder reflection for current column
		beta[k] = makeHouseholder!T(m[k+1..$, k]);

		// update remaining columns
		applyHouseholder!T(m[k+1..$, k+1..$], m[k+2..$, k], beta[k]);

		// update all(!) rows
		applyHouseholderRight!T(m[0..$, k+1..$], m[k+2..$, k], beta[k]);
	}
}


//////////////////////////////////////////////////////////////////////
/// Eigenvalue decomposition
//////////////////////////////////////////////////////////////////////

/**
 * Schur decomposition: A = Q * U * Q^-1
 *   - m is replaced by the triangular matrix U (diagonal in case of A being hermitian)
 *   - if A is real with complex eigenvalues, U is only block-triangular
 */
void denseComputeSchur(T)(ContiguousMatrix!T m, ContiguousMatrix!T v = ContiguousMatrix!T.init)
{
	int n = cast(int)m.length!0;
	bool trafo = v.length!0 != 0 || v.length!1 != 0;
	if(m.length!1 != n || trafo && (v.length!0 != n || v.length!1 != n))
		throw new Exception("matrix dimension mismatch");

	// reduce m to Hessenberg matrix
	auto alpha = new RealTypeOf!T[n];
	denseComputeHessenberg!T(m, alpha);
	if(trafo)
		for(int i = 0; i < n-1; ++i)
			applyHouseholderRight!T(v[0..$,i+1..$], m[i+2..$,i], alpha[i]);
	delete alpha;
	for(int j = 0; j < n; ++j)
		m[min(j+2,n)..n, j] = T(0);

	// Francis algorithm
	auto eps = 4*RealTypeOf!T.epsilon;
	int steps = 0;
	for(size_t b = n; b > 0; )
	{
		// find trailing block of non-zero off-diagonal
		long a;
		for(a = b-1; a > 0; --a)
			if(abs(m[a, a-1]) <= eps*(abs(m[a, a]) + abs(m[a-1, a-1])))
			{
				m[a, a-1] = T(0);
				break;
			}

		// single value -> nothing to be done
		if(a == b-1)
		{
			steps = 0;
			b -= 1;
			continue;
		}

		// 2x2 block -> solve it directly
		if(a == b-2)
		{
			// directly compute one eigenvalue
			auto p = -m[a,a] - m[a+1,a+1];
			auto q = m[a,a]*m[a+1,a+1] - m[a,a+1]*m[a+1,a];
			auto root = polyRoot(p, q);

			T x,y;
			RealTypeOf!T beta;

			// complex roots of real matrix -> leave the block as it is
			// TODO: I think there is a normal-form for this case
			if(!is(T : Complex!R, R) && root.im != 0)
			{
				steps = 0;
				b -= 2;
				continue;
			}

			x = m[a,a] - root.re;
			y = m[a+1,a];
			x += phase(x)*sqrt(sqAbs(x) + sqAbs(y));
			beta = 2/(sqAbs(x) + sqAbs(y));

			for(size_t i = a; i < n; ++i)
			{
				auto s = conj(x)*m[a,i] + conj(y)*m[a+1,i];
				m[a,i] -= beta*s*x;
				m[a+1,i] -= beta*s*y;
			}
			for(size_t i = 0; i <= a+1; ++i)
			{
				auto s = x*m[i,a] + y*m[i,a+1];
				m[i,a] -= beta*s*conj(x);
				m[i,a+1] -= beta*s*conj(y);
			}
			if(trafo)
			for(size_t i = 0; i < n; ++i)
			{
				auto s = x*v[i,a] + y*v[i,a+1];
				v[i,a] -= beta*s*conj(x);
				v[i,a+1] -= beta*s*conj(y);
			}

			// In precise arithmetic, the current 2x2 block would be solved now.
			// But that seems to be unreliable in practice. So don't decrease b,
			// which means doing another step when neccessary. TODO: figure this out.

			//steps = 0;
			//b -= 2;
			continue;
		}

		// larget block -> do a francis step
		if(++steps > 30)
		{
			import std.stdio;
			writefln("!!!!!!!!!!!!!!!!!!!!! TIMEOUT n = %s !!!!!!!!!!!!!!!!!!!",n);
			return;
		}

		// chararacteristic polynomial of trailing 2x2 matrix
		auto p = -m[b-1,b-1]-m[b-2,b-2];
		auto q = m[b-1,b-1]*m[b-2,b-2] - m[b-1,b-2]*m[b-2,b-1];

		// first column of p(A)
		T[3] _x;
		auto x = _x[].sliced(3);
		x[0] = m[0,0]*m[0,0] + m[0,1]*m[1,0] + p*m[0,0] + q;
		x[1] = m[1,0]*m[0,0] + m[1,1]*m[1,0] + p*m[1,0];
		x[2] = m[2,1]*m[1,0];

		// make a householder reflection of x (i.e. create the bulge)
		auto beta = makeHouseholder!T(x);
		applyHouseholder!T(m[0..3,0..$], x[1..$], beta);
		applyHouseholderRight!T(m[0..min(4,$),0..3], x[1..$], beta);
		if(trafo)
			applyHouseholderRight!T(v[0..$,0..3], x[1..$], beta);

		// make the matrix Hessenberg again (i.e. chase the bulge down)
		for(int i = 0; i < n-2; ++i)
		{
			beta = makeHouseholder!T(m[i+1..min(i+4,$), i]);
			applyHouseholder!T(m[i+1..min(i+4,$),i+1..$], m[i+2..min(i+4,$),i], beta);
			applyHouseholderRight!T(m[0..$,i+1..min(i+4,$)], m[i+2..min(i+4,$),i], beta);
			if(trafo)
				applyHouseholderRight!T(v[0..$,i+1..min(i+4,$)], m[i+2..min(i+4,$),i], beta);
			m[i+2 .. min(i+4,n), i] = T(0);
		}
	}
}


//////////////////////////////////////////////////////////////////////
/// private helpers
//////////////////////////////////////////////////////////////////////

/** complex conjugate of a scalar */
private T conj(T)(T x)
	if(!isVector!T)
{
	static if(is(T : Complex!R, R))
		return std.complex.conj(x);
	else
		return x;
}

/** complex conjugate of a vector */
private auto conj(V)(V a)
	if(isVector!V)
{
	alias T = DeepElementType!V;
	static if(is(T : Complex!R, R))
		return a.map!(std.complex.conj);
	else
		return a;
}

/** adjoint of a matrix */
private auto adjoint(M)(M a)
{
	alias T = DeepElementType!M;
	static if(is(T : Complex!R, R))
		return a.transposed.map!(std.complex.conj);
	else
		return a.transposed;
}

/** dot-product of two vectors without any complex conjugation */
private auto dot(V,W)(V a, W b)
	if(isVector!V && isVector!W)
{
	alias T = typeof(DeepElementType!V.init * DeepElementType!W.init);
	return reduce!"a+b*c"(T(0), a, b);
}

/** same as dot(v.conj, v) */
private RealTypeOf!(DeepElementType!V) sqNorm2(V)(V a)
	if(isVector!V)
{
	alias R = RealTypeOf!(DeepElementType!V);
	return R(0).reduce!"a+b"(a.map!(std.complex.sqAbs));
}
