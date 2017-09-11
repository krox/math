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
private import math.numtheory;

/**
 * General implementation ideas:
 *
 * The finite field F_(p^n) can be constructed as F_p[x]/f(x) where f is an
 * irreducible polynomial of degree n, so the straight forward way to represent
 * an element of the field would be to store an array of coefficients of the
 * representative polynomial. A more sophisticated approach is to note that the
 * unit group of a finite group is always cyclic. Therefore we can choose one
 * fixed generator z and represent all elements (except zero) as z^r for some z
 * in [0,p^n-2]. By choosing f as a Conway polynomial we can assure that z=x
 * is always a valid choice, i.e. x is a generating element of F_p[x]/f(x).
 *
 * In this module we choose the second representation.
 *   + individual elements are just are trivial (just an integer)
 *   + converting to polynomial representation is easy ("binary exponentiation")
 *   + multiplication is trivial (just add the exponents)
 *   - converting from polynomial representation is hard ("discrete logarithm"),
 *     so to make it efficient, we need big tables
 *
 * We use the magic value r=-1 for the zero element, and thus the convention
 * log(0) = -1, which makes at least a little sense, as in real arithmetic
 * log(0) = -infinity.
 */

/** actually only immutable instances of this should ever be created */
class FiniteField
{
	int p, q, n; // order of the field q = p^^n
	FFE z;	// generator of the unit group

	immutable(int)[] poly;	// conway polynomial used for construction
	private int[] expTable; // more like int[][]
	private int[] logTable; // entry 0 not used

	/**
	 * constrcutor which computes the static exp/log tables. Should not be
	 * called directly. Use FiniteField.opCall for a cached version.
	 */
	private this(int p, int n) pure
	{
		if(!isPrime(p) || n < 1)
			throw new Exception("creating finite field with invalid parameters");

		this.p = p;
		this.n = n;
		this.q = p^^n;

		poly = conwayPolynomial(p,n);
		assert(poly.length == n+1);
		assert(poly[$-1] == 1);

		expTable = new int[n*(q-1)];
		logTable = new int[q];
		expTable[0] = 1;
		logTable[0] = -1;

		for(int i = 1; i < q-1; ++i)
		{
			expTable[i*n+1 .. (i+1)*n] = expTable[(i-1)*n .. i*n-1];
			expTable[i*n   .. (i+1)*n] -= expTable[i*n-1] * poly[0..$-1];
			expTable[i*n   .. (i+1)*n] %= p;
			expTable[i*n   .. (i+1)*n] += p;
			expTable[i*n   .. (i+1)*n] %= p;

			int j = 0;
			foreach_reverse(x; expTable[i*n .. (i+1)*n])
				j = p*j+x;
			logTable[j] = i;
		}

		if(q == 2)
			z = FFE.fromExp(cast(immutable)this, 0);
		else
			z = FFE.fromExp(cast(immutable)this, 1);
	}

	private static Rebindable!(immutable(FiniteField))[int] cache;

	/** create a finite field of order p^n or retrieve one from cache */
	static immutable(FiniteField) opCall(int p, int n)
	{
		int q = p^^n;

		if(q in cache)
			return cache[q];

		auto f = new immutable(FiniteField)(p, n);
		cache[q] = f;
		return f;
	}

	FFE random() immutable @property
	{
		return FFE.fromExp(this, uniform(-1,q-1));
	}

	private immutable(int)[] exp(int r) immutable pure nothrow
	{
		assert(r >= 0);
		assert(r < q-1);
		return expTable[r*n..(r+1)*n];
	}

	// log(0) = -1, its kinda right
	private int log(int j) immutable pure nothrow
	{
		assert(0 <= j && j < q);
		return logTable[j];
	}
}

struct FFE
{
    // power of generating element.
    // r = -1 is the "zero" element
    // r = 0 is the "one" element
    // in these two cases, field may be null. This simplifies some usecases
    // as it enables the constructor calls FFE(0) and FFE(1).
    int r = -1;
    Rebindable!(immutable(FiniteField)) field = null;

	//////////////////////////////////////////////////////////////////////
	/// constructors
	//////////////////////////////////////////////////////////////////////

	/**
	 * constructor that embeds integers into the finite field. Note that not
	 * all elements are constructable this way (unless n = 1).
	 */
	this(int v, immutable(FiniteField) field = null) pure nothrow
	{
        this.field = field;

		if(v == 0)
			r = -1;
		else if(v == 1)
			r = 0;
		else
		{
			assert(field !is null, "generic finite field elements may only be zero or one");
			v %= field.p;
			v += field.p;
			v %= field.p;
			r = field.log(v);
		}
	}

    /** z^r, where z is the designated generator of the unit group */
    static fromExp(immutable(FiniteField) field, int r) pure nothrow
    {
        FFE x;
        if(field is null)
            assert(r == -1 || r == 0);
        else
            assert(-1 <= r && r <= field.q-2);
		x.field = field;
		x.r = r;
        return x;
	}

	//////////////////////////////////////////////////////////////////////
	/// unary operations
	//////////////////////////////////////////////////////////////////////

	/** additive inverse */
	FFE opUnary(string op)() const pure
		if(op == "-")
	{
		if(r == -1)
			return this; // -0 = 0

		assert(field !is null);

		auto p = field.exp(r);
		int j = 0;
		foreach_reverse(x; p)
			j = field.p*j + (field.p-x)%field.p;

		return fromExp(field, field.log(j));
	}

    /** multiplicative inverse */
    FFE inverse() const @property pure
	{
		if(r == -1)
			throw new Exception("division by zero (in a finite field)");
		if(r == 0)
			return this; // 1/1 = 1

		assert(field !is null);
		return fromExp(field, (field.q-1-r)%(field.q-1));
	}

    /** power */
    FFE opBinary(string op)(int e) const pure
        if(op == "^^")
	{
		if(r == -1)
		{
			if(e == 0)
				return fromExp(field, 0);	// 0^0 = 1 by convention
			else
				return fromExp(field, -1);	// 0^e = 0 for e != 0
		}

        assert(field !is null);
		e %= field.q-1;
		return fromExp(field, (r*e) % (field.q-1));
	}

	//////////////////////////////////////////////////////////////////////
	/// binary operations
	//////////////////////////////////////////////////////////////////////

	FFE opBinary(string op)(FFE b) const pure
		if(op == "+" || op == "-")
	{
        auto f = common(this.field, b.field);

		if(b.r == -1) // x +- 0 = x
			return fromExp(f, r);

		if(r == -1)
		{
			static if(op == "+")
				return fromExp(f, b.r);   // 0 + x = x
			else
				return -fromExp(f, b.r);  // 0 - x = -x
		}


        assert(f !is null);

		auto x = f.exp(r);
		auto y = f.exp(b.r);
		assert(x.length == y.length);
		int z = 0;
		for(int i = cast(int)x.length-1; i >= 0; --i)
		{
			static if(op == "+")
				z = f.p*z + (x[i] + y[i]) % f.p;
			else
				z = f.p*z + (x[i] - y[i] + f.p) % f.p;
		}

		return fromExp(f, f.log(z));
	}

	FFE opBinary(string op)(FFE b) const pure
		if(op == "*")
	{
		auto f = common(this.field, b.field);

		if(r == -1 || b.r == -1)
			return fromExp(f, -1);

        assert(f !is null);
        return fromExp(f, (r + b.r) % (f.q - 1));
	}

	FFE opBinary(string op)(FFE b) const pure
		if(op == "/")
	{
		return this*b.inverse;
	}


	//////////////////////////////////////////////////////////////////////
	/// convenience wrappers for opAssign ans FFE<->int operations
	//////////////////////////////////////////////////////////////////////

	FFE opOpAssign(string op)(FFE b) pure
	{
		return this = this.opBinary!op(b);
	}

	FFE opOpAssign(string op)(int b) pure
	{
		return this = this.opBinary!op(b);
	}

    FFE opBinary(string op)(int b) const pure
        if(op == "+" || op == "-" || op == "*" || op == "/")
	{
		// TODO: + and - could be faster without this lowering
		return this.opBinary!op(FFE(b, field));
	}

    FFE opBinaryRight(string op)(int a) const pure
    {
        return FFE(a, field).opBinary!op(this);
    }


    //////////////////////////////////////////////////////////////////////
    /// misc stuff
    //////////////////////////////////////////////////////////////////////

	bool opEquals(int b) const pure nothrow
	{
        return this.opEquals(FFE(b, field));
	}

	bool opEquals(FFE b) const pure nothrow
	{
		cast(void)common(this.field, b.field);
		return r == b.r;
	}

	/**
	 * Finite fields can not be ordered in a way compatible with the field
	 * structure, but its still useful to have a canonical order.
	 */
	int opCmp(FFE b) const pure nothrow
	{
		cast(void)common(this.field, b.field);
		return r - b.r;
	}

	/** (somewhat) human readable form */
	string toString() const pure @property
	{
		if(r == -1)
			return "0";
		if(r == 0)
			return "1";

        assert(field !is null);
		if(field.n == 1)
			return "["~to!string(field.exp(r)[0])~"]";
		else
			return "z^"~to!string(r);
	}

    private static T common(T)(T a, T b) pure nothrow
    {
        if(a is null)
            return b;
        if(b is null)
            return a;
        if(a !is b)
        	throw new Error("incompatible finite fields");
        return a;
    }
}

unittest
{
	auto F = FiniteField(3,2);
	auto x = FFE(1,F);
	assert((x+x)^^3 == FFE.fromExp(F, 4));
}
