module math.linear;

/**
 * linear algebra backend
 */

import jive.array;
private import std.math;
private import std.traits;
private import std.algorithm : move, swap;
private import jive.bitarray;

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
			T s = m[l,k] = m[l,k] / m[k,k]; // this is now part of the L matrix
			for(int i = k+1; i < n; ++i)
				m[l,i] = m[l,i] - s * m[k,i];
		}
	}
}

/** solve linear equation after LU decomposition was computed */
void denseSolveLU(T)(const Slice2!T m, Slice2!T b)
{
	int n = cast(int)m.size[0];
	if(m.size[1] != n || b.size[0] != n)
		throw new Exception("matrix dimension mismatch");

	for(int k = 0; k < n; ++k)
	{
		// no division here as L has implicit 1's on diagonal

		for(int l = k+1; l < n; ++l)
			for(int i = 0; i < b.size[1]; ++i)
				b[l,i] = b[l,i] - m[l,k] * b[k,i];
	}

	for(int k = n-1; k >= 0; --k)
	{
		for(int i = 0; i < b.size[1]; ++i)
			b[k,i] = b[k,i] / m[k,k];

		for(int l = k-1; l >= 0; --l)
			for(int i = 0; i < b.size[1]; ++i)
				b[l,i] = b[l,i] - m[l,k] * b[k,i];
	}
}
