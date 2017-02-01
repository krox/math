module math.linear;

/**
 * linear algebra backend
 */

private import jive.array;
private import std.math;
private import std.traits;
private import std.complex;
private import std.algorithm;
private import jive.bitarray;
private import math.numerics;

/** Matrix decompositions

name				decomposition 		notes

LU 			PA = LU				L has 1's on the diagonal
QR			A = QR
Hessenberg	A = QHQ^-1
Schur		A = QUQ^-1			does not always exist for real matrices
Cholesky	A = LDL^T			only for Hermitian A

*/

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
void denseSolveL(T)(Slice2!(const(T)) m, Slice2!T b)
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
void denseSolveLU(T)(Slice2!(const(T)) m, Slice2!T b)
{
	int n = cast(int)m.size[0];
	if(m.size[1] != n || b.size[0] != n)
		throw new Exception("matrix dimension mismatch");

	denseSolveL!T(m, b);
	denseSolveU!T(m, b);
}


//////////////////////////////////////////////////////////////////////
/// QR decomposition
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
		static if(conjugate)
			r += conj(a[i])*b[i];
		else
			r += a[i]*b[i];
	return r;
}

/** compute householder reflection inplace, returns beta factor */
RealTypeOf!T makeHouseholder(T)(Slice!T v)
{
	T c = -phase(v[0])*sqrt(norm2!T(v));

	if(c == 0)	// typically this indicate a singular matrix
		return 0;
		//throw new Exception("zero in householder trafo");

	for(int i = 1; i < v.length; ++i)
		v[i] /= v[0] - c;
	v[0] = c;
	return 2/(1 + norm2(v[1..$]));
}

/** a = (1 - beta*|v><v|) * a */
void applyHouseholder(T)(Slice2!T a, Slice!(const(T)) v, RealTypeOf!T beta)
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
void applyHouseholderRight(T)(Slice2!T a, Slice!(const(T)) v, RealTypeOf!T beta)
{
	int n = cast(int)v.size[0]+1;
	if(a.size[1] != n)
		throw new Exception("vector dimension mismatch");

	for(int k = 0; k < a.size[0]; ++k)
	{
		T s = beta * (a[k,0] + scalarProduct!(T,false)(v, a[k,1..$]));
		a[k, 0] -= s;
		for(int i = 1; i < n; ++i)
			a[k, i] -= s * conj(v[i-1]);
	}
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

void applyGivens(T)(Slice2!T m, int i, int j, RealTypeOf!T c, ref T s)
{
	for(int k = 0; k < m.size[1]; ++k)
	{
		auto tmp = c*m[i,k] + s*m[j,k];
		m[j,k] = -conj(s)*m[i,k] + c*m[j,k];
		m[i,k] = tmp;
	}
}

void applyGivensRight(T)(Slice2!T m, int i, int j, RealTypeOf!T c, ref T s)
{
	for(int k = 0; k < m.size[1]; ++k)
	{
		auto tmp = c*m[k,i] + conj(s)*m[k,j];
		m[k,j] = -s*m[k,i] + c*m[k,j];
		m[k,i] = tmp;
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
		beta[k] = makeHouseholder!T(m[k..$, k]);

		// update remaining columns
		applyHouseholder!T(m[k..$, k+1..$], m[k+1..$, k], beta[k]);
	}
}

/** solve linear equation after QR decomposition was computed */
void denseSolveQR(T)(Slice2!(const(T)) m, const(RealTypeOf!T)[] beta, Slice2!T b)
{
	int n = cast(int)beta.length;
	if(m.size[0] != n || m.size[1] != n || b.size[0] != n)
		throw new Exception("matrix dimension mismatch");

	for(int k = 0; k < n; ++k)
		applyHouseholder!T(b[k..$,0..$], m[k+1..$,k], beta[k]);

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
		beta[k] = makeHouseholder(m[k+1..$, k]);

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
void denseComputeSchur(T)(Slice2!T m, Slice2!T v = Slice2!T.init)
{
	size_t n = m.size[0];
	if(m.size[1] != n || v.ptr !is null && (v.size[0] != n || v.size[1] != n))
		throw new Exception("matrix dimension mismatch");

	// initialize q to identity matrix
	if(v.ptr !is null)
		for(int j = 0; j < n; ++j)
			for(int i = 0; i < n; ++i)
				v[i,j] = i==j?1:0;

	// reduce m to Hessenberg matrix
	auto alpha = new RealTypeOf!T[n];
	denseComputeHessenberg!T(m, alpha);
	if(v.ptr !is null)
		for(int i = 0; i < n-1; ++i)
			applyHouseholderRight!T(v[0..$,i+1..$], m[i+2..$,i], alpha[i]);
	for(int j = 0; j < n; ++j)
		for(int i = j+2; i < n; ++i)
			m[i,j] = 0;

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
				m[a, a-1] = 0;
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
			if(!isComplex!T && root.im != 0)
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
			if(v.ptr !is null)
			for(size_t i = 0; i < n; ++i)
			{
				auto s = x*v[i,a] + y*v[i,a+1];
				v[i,a] -= beta*s*conj(x);
				v[i,a+1] -= beta*s*conj(y);
			}

			// In precise arithmetic, the current 2x2 block would be solved now.
			// But that seems to be kinda. So don't decrease b, which means
			// doing another step when neccessary. TODO: figure this out.

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
		auto x = Slice!T(3);
		x[0] = m[0,0]*m[0,0] + m[0,1]*m[1,0] + p*m[0,0] + q;
		x[1] = m[1,0]*m[0,0] + m[1,1]*m[1,0] + p*m[1,0];
		x[2] = m[2,1]*m[1,0];

		// make a householder reflection of x (i.e. create the bulge)
		auto beta = makeHouseholder(x);
		applyHouseholder!T(m[0..3,0..$], x[1..$], beta);
		applyHouseholderRight!T(m[0..min(4,$),0..3], x[1..$], beta);
		if(v.ptr !is null)
			applyHouseholderRight!T(v[0..$,0..3], x[1..$], beta);

		// make the matrix Hessenberg again (i.e. chase the bulge down)
		for(int i = 0; i < n-2; ++i)
		{
			beta = makeHouseholder(m[i+1..min(i+4,$), i]);
			applyHouseholder!T(m[i+1..min(i+4,$),i+1..$], m[i+2..min(i+4,$),i], beta);
			applyHouseholderRight!T(m[0..$,i+1..min(i+4,$)], m[i+2..min(i+4,$),i], beta);
			if(v.ptr !is null)
				applyHouseholderRight!T(v[0..$,i+1..min(i+4,$)], m[i+2..min(i+4,$),i], beta);
			for(int j = i+2; j < min(i+4,n); ++j)
				m[j,i] = 0;
		}
	}
}
