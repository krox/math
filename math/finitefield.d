module math.finitefield;

/**
 * Finite fields of order p^n with p prime.
 *
 * Should be more efficient, and possibly simpler to use than an explicit
 * construction using math.integer, math.factorring and math.polynomial,
 * which would be along the lines of Coset!(Polynomial!(Coset!Integer)).
 *
 * Note that while an individual Element takes very little space, the field
 * itself precomputes tables of size O(n*p^n), which make calculations very 
 * fast, but make this module unsuitable for very large fields.
 */

private import std.conv : to;
private import std.exception : assumeUnique;
private import std.typecons : Rebindable;
private import std.random : uniform;

private import math.conway;
private import math.prime;

class FiniteField
{
	immutable int p, q, n; // q = p^^n
	immutable FFE z;

	immutable int[] poly;
	immutable int[] expTable; // more like int[][]
	immutable int[] logTable; // entry 0 not used

	private this(int p, int n) immutable pure
	{
		this.p = p;
		this.n = n;
		this.q = p^^n;

		poly = conwayPolynomial(p,n);
		assert(poly.length == n+1);
		assert(poly[$-1] == 1);
		{
			auto tmp = new int[n*(q-1)];
			tmp[0] = 1;

			for(int i = 1; i < q-1; ++i)
			{
				tmp[i*n+1 .. (i+1)*n] = tmp[(i-1)*n .. i*n-1];
				tmp[i*n   .. (i+1)*n] -= tmp[i*n-1] * poly[0..$-1];
				tmp[i*n   .. (i+1)*n] %= p;
				tmp[i*n   .. (i+1)*n] += p;
				tmp[i*n   .. (i+1)*n] %= p;
			}

			expTable = assumeUnique(tmp);
		}

		{
			auto tmp = new int[q];

			tmp[0] = -1; // log(0) = -1 is a little bit correct

			for(int i = 1; i < q-1; ++i)
			{
				int j = 0;
				foreach_reverse(x; expTable[i*n .. (i+1)*n])
				{
					j *= p;
					j += x;
				}
				tmp[j] = i;
			}

			logTable = assumeUnique(tmp);
		}

		z = FFE(this, 1);
	}

	private static immutable(FiniteField)[int] cache;

	static immutable(FiniteField) opCall(int p, int n)
	{
		if(!isPrime(p) || n < 1)
			throw new Exception("creating finite field with invalid parameters");
		int q = p^^n;

		if(q in cache)
			return cache[q];

		return cache[q] = new immutable(FiniteField)(p, n);
	}

	FFE opCall(int x) immutable pure
	{
		x %= p;
		x += p;
		x %= p;
		if(x == 0)
			return FFE(this, -1);
		else
			return FFE(this, log(x));
	}

	FFE random() immutable @property
	{
		return FFE(this, uniform(-1,q-1));
	}

	private immutable(int)[] exp(int r) immutable pure
	{
		assert(r >= 0);
		assert(r < q-1);
		return expTable[r*n..(r+1)*n];
	}

	// log(0) = -1, its kinda right
	private int log(int p) immutable pure
	{
		assert(p >= 0);
		return logTable[p];
	}
}

struct FFE
{
	Rebindable!(immutable(FiniteField)) field;
	int r; // r = -1 means zero

	private this(immutable(FiniteField) field, int r) pure
	{
		this.field = field;
		this.r = r;
	}

	pure invariant()
	{
		if(field is null)
			return;
		if(r < -1 || r >= field.q-1)
			throw new Exception("finite field operations are broken");
	}

	FFE opUnary(string op)() const pure
		if(op == "-")
	{
		if(r == -1) // -0 = 0
			return this;

		auto x = field.exp(r);
		int z = 0;
		for(int i = cast(int)x.length-1; i >= 0; --i)
		{
			z *= field.p;
			z += (field.p - x[i]) % field.p;
		}

		return FFE(field, field.log(z));
	}

	FFE inverse() const @property pure
	{
		if(r == -1)
			throw new Exception("division by zero (in a finite field)");

		return FFE(field, (field.q-1-r)%(field.q-1));
	}

	FFE opBinary(string op)(int b) const pure
		if(op == "+" || op == "-" || op == "*" || op == "/")
	{
		// TODO: + and - could be faster without this lowering
		return opBinary!op(field(b));
	}

	FFE opBinary(string op)(FFE b) const pure
		if(op == "+")
	{
		if(field !is b.field)
			throw new Exception("tried to do "~op~" on elements of different finite fields");

		if(r == -1)
			return b;
		if(b.r == -1)
			return this;

		auto x = field.exp(r);
		auto y = field.exp(b.r);
		assert(x.length == y.length);
		int z = 0;
		for(int i = cast(int)x.length-1; i >= 0; --i)
		{
			z *= field.p;
			z += (x[i] + y[i]) % field.p;
		}

		return FFE(field, field.log(z));
	}

	FFE opBinary(string op)(FFE b) const pure
		if(op == "-")
	{
		if(field !is b.field)
			throw new Exception("tried to do "~op~" on elements of different finite fields");

		if(r == -1)
			return -b;
		if(b.r == -1)
			return this;

		auto x = field.exp(r);
		auto y = field.exp(b.r);
		assert(x.length == y.length);
		int z = 0;
		for(int i = cast(int)x.length-1; i >= 0; --i)
		{
			z *= field.p;
			z += (x[i] - y[i] + field.p) % field.p;
		}

		return FFE(field, field.log(z));
	}

	FFE opBinary(string op)(FFE b) const pure
		if(op == "*")
	{
		if(field !is b.field)
			throw new Exception("tried to do "~op~" on elements of different finite fields");

		if(r == -1 || b.r == -1)
			return FFE(field, -1);
		return FFE(field, (r + b.r) % (field.q - 1));
	}

	FFE opBinary(string op)(FFE b) const pure
		if(op == "/")
	{
		if(field !is b.field)
			throw new Exception("tried to do "~op~" on elements of different finite fields");

		if(b.r == -1)
			throw new Exception("division by zero (in a finite field)");
		if(r == -1)
			return FFE(field, -1);
		return FFE(field, (r - b.r + field.q - 1) % (field.q - 1));
	}

	FFE opBinary(string op)(int e) const pure
		if(op == "^^")
	{
		if(r == -1)
			if(e == 0) // TODO: e% == 0, maybe make it undefined
				return FFE(field, 0);
			else
				return FFE(field, -1);
		e = e % (field.q-1) + (field.q-1);
		return FFE(field, (r*e) % (field.q-1));
	}

	bool opEquals(int b) const pure
	{
		if(r == -1)
			return b == 0;

		auto x = field.exp(r);
		foreach(c; x[1..$])
			if(c != 0)
				return false;
		b %= field.p;
		b += field.p;
		b %= field.p;
		return x[0] == b;
	}

	bool opEquals(FFE b) const pure
	{
		if(field !is b.field)
			throw new Exception("tried to compare elements of different finite fields");

		return r == b.r;
	}

	/**
	 * Finite fields can not be ordered in a way compatible with the field
	 * structure, but its still useful to have a canonical order.
	 */
	int opCmp(FFE b) const pure
	{
		return r - b.r;
	}

	string toString() const @property pure
	{
		if(r == -1)
			return "0";

		if(field.n == 1)
			return "["~to!string(field.exp(r)[0])~"]";
		else
			return "z^"~to!string(r);
	}
}
