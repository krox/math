module math.integer;

private import std.string : toStringz;
private import std.conv : to;
private import std.typecons;
private import std.exception : assumeUnique;
private import std.random : unpredictableSeed;

private import math.gmp;

/**
 * BigInteger type using the GMP library.
 * It is implemented with copy-on-write semantics.
 * Should be mostly compatible with int. Differences include:
 *   * There is a NaN value, to which all variables are initialized.
 *     Note that it is an error to use a NaN in computations.
 *   * Integer division rounds the quotient towards -infinity.
 *     In particular, (a % b) always has the same sign as b.
 */
struct Integer
{
	//////////////////////////////////////////////////////////////////////
	/// constructors
	//////////////////////////////////////////////////////////////////////

	static enum nan = Integer.init;

	this(immutable GmpInteger z)
	{
		this.z = z;
	}

	/** constructor for given value */
	this(int v)
	{
		if(0 <= v && v < cacheSize)
		{
			this(cache[v]);
			return;
		}

		auto r = new GmpInteger;
		__gmpz_set_si(r.ptr, v);
		this(cast(immutable)r);
	}

	/** ditto */
	this(double v)
	{
		auto r = new GmpInteger;
		__gmpz_set_d(r.ptr, v);
		this(cast(immutable)r);
	}

	/** ditto */
	this(string v)
	{
		// TODO: throw exception on bad strings
		auto r = new GmpInteger;
		__gmpz_init_set_str(r.ptr, toStringz(v), 0);
		this(cast(immutable)r);
	}

	static Integer random(Integer n)
	{
		auto r = new GmpInteger;
		__gmpz_urandomm(r.ptr, &rand, n.ptr);
		return Integer(cast(immutable)r);
	}

	static Integer randomBits(size_t b)
	{
		auto r = new GmpInteger;
		__gmpz_urandomb(r.ptr, &rand, b);
		return Integer(cast(immutable)r);
	}


	//////////////////////////////////////////////////////////////////////
	/// conversion
	//////////////////////////////////////////////////////////////////////

	/** returns value as (decimal) string */
	string toString() const @property
	{
		auto buflen = __gmpz_sizeinbase(ptr, 10)+2;	// one for sign, one for \0
		auto buf = new char[buflen];
		return to!string(__gmpz_get_str(buf.ptr, 10, ptr));
	}

	/** returns value as int, assuming it is small enough */
	int opCast(T)() const
		if(is(T == int))
	{
		return __gmpz_get_si(z);
	}

	/** test if bit i is set */
	bool opIndex(size_t i) const
	{
		return __gmpz_tstbit(ptr, i) != 0;
	}


	//////////////////////////////////////////////////////////////////////
	/// size metric and comparisons
	//////////////////////////////////////////////////////////////////////

	bool isNan() const @property
	{
		return z is null;
	}

	bool opEquals(int b) const
	{
		return __gmpz_cmp_si(ptr, b) == 0;
	}

	bool opEquals(Integer b) const
	{
		return __gmpz_cmp(ptr, b.ptr) == 0;
	}

	int opCmp(int b) const
	{
		return __gmpz_cmp_si(ptr, b);
	}

	int opCmp(Integer b) const
	{
		return __gmpz_cmp(ptr, b.ptr);
	}

	/** returns -1 / 0 / +1, possibly faster than actual compare */
	int sign() const @property
	{
		return z.z._mp_size < 0 ? -1 : z.z._mp_size > 0;
	}

	/** number of bits */
	size_t length() const @property
	{
		if(sign < 0)
			throw new Exception("negative integer dont have a length (actually more like infinity");
		if(sign == 0)
			return 0;
		return __gmpz_sizeinbase(ptr, 2);
	}


	//////////////////////////////////////////////////////////////////////
	/// arithmetic operations
	//////////////////////////////////////////////////////////////////////

	Integer opUnary(string op)() const
		if(op == "-")
	{
		auto r = new GmpInteger;
		__gmpz_neg(r.ptr, ptr);
		return Integer(cast(immutable)r);
	}

	Integer opBinary(string op)(int b) const
		if(op != "%")
	{
		auto r = new GmpInteger;

		static if(op == "+")
		     if(b >= 0) __gmpz_add_ui(r.ptr, ptr, b);
		     else       __gmpz_sub_ui(r.ptr, ptr, -b);
		else static if(op == "-")
			if(b >= 0) __gmpz_sub_ui(r.ptr, ptr, b);
			else       __gmpz_add_ui(r.ptr, ptr, -b);
		else static if(op == "*")
			__gmpz_mul_si(r.ptr, ptr, b);
		else static if(op == "/")
		{
			if(b>0)
				__gmpz_fdiv_q_ui(r.ptr, ptr, b);
			else
				throw new Exception("TODO");
		}
		else static if(op == "^^")
		{
			if(b < 0)
				throw new Exception("negative powers of integers dont exist");
			__gmpz_pow_ui(r.ptr, ptr, b);
		}
		else static if(op == ">>")
		{
			if(b < 0)
				throw new Exception("negative shift amounts not supported");
			__gmpz_fdiv_q_2exp(r.ptr, ptr, b);
		}
		else static if(op == "<<")
		{
			if(b < 0)
				throw new Exception("negative shift amounts not supported");
			__gmpz_mul_2exp(r.ptr, ptr, b);
		}

		else static assert(false, "binary '"~op~"' is not defined");

		return Integer(cast(immutable)r);
	}

	int opBinary(string op)(int b) const
		if(op == "%")
	{
		if(b <= 0)
			throw new Exception("only positive modulus is supported");
		scope r = new GmpInteger;
		return cast(int)__gmpz_fdiv_r_ui(r.ptr, ptr, b);
	}

	Integer opBinaryRight(string op)(int a) const
		if(op == "+" || op == "*")
	{
		return opBinary!op(a); // commutative operators can simply be forwarded
	}

	Integer opBinaryRight(string op)(int a) const
		if(op == "-")
	{
		auto r = new GmpInteger;
		if(a >= 0)
			__gmpz_ui_sub(r.ptr, a, ptr);
		else
		{
			__gmpz_add_ui(r.ptr, ptr, -a);
			__gmpz_neg(r.ptr, r.ptr);
		}

		return Integer(cast(immutable)r);
	}

	Integer opBinary(string op)(Integer b) const
	{
		auto r = new GmpInteger;

		     static if(op == "+") __gmpz_add(r.ptr, ptr, b.ptr);
		else static if(op == "-") __gmpz_sub(r.ptr, ptr, b.ptr);
		else static if(op == "*") __gmpz_mul(r.ptr, ptr, b.ptr);
		else static if(op == "/") __gmpz_fdiv_q(r.ptr, ptr, b.ptr);
		else static if(op == "%") __gmpz_fdiv_r(r.ptr, ptr, b.ptr);
		else static assert(false, "binary '"~op~"' is not defined");

		return Integer(cast(immutable)r);
	}

	/** returns this/b. Faster, but only works if the division is exact (i.e. no rounding) */
	Integer divExact(Integer b) const
	{
		auto r = new GmpInteger;
		__gmpz_divexact(r.ptr, ptr, b.ptr);
		return Integer(cast(immutable)r);
	}


	//////////////////////////////////////////////////////////////////////
	/// internals
	//////////////////////////////////////////////////////////////////////

	Rebindable!(immutable(GmpInteger)) z;

	immutable(mpz_t)* ptr() const @property
	{
		return &z.z;
	}

	private enum cacheSize = 10;

	private static const immutable(GmpInteger)[cacheSize] cache;
	private static mp_randstate_t rand;

	static this()
	{
		for(int i = 0; i < cacheSize; ++i)
		{
			auto z = new GmpInteger;
			__gmpz_init_set_si(z.ptr, i);
			cache[i] = cast(immutable)z;
		}

		__gmp_randinit_default(&rand);
		__gmp_randseed_ui(&rand, unpredictableSeed);
	}
}

/**
 * calculate modular inverse a^-1 % m
 */
Integer inverseMod(Integer a, Integer m)
{
	auto r = new GmpInteger;
	__gmpz_invert(r.ptr, a.ptr, m.ptr);
	return Integer(cast(immutable)r);
}

/**
 * returns true if a is a square number (this includes 0 and 1)
 */
bool isPerfectSquare(Integer a)
{
	return __gmpz_perfect_square_p(a.ptr) != 0;
}

/**
 * checks wether a is a prime number returns
 *   - 0 on composite
 *   - 1 on probable prime
 *   - 2 definite prime
 * chance of error should be negligible
 */
int isPrime(Integer a)
{
	// 25 is the number of round suggested by the GMP manual. Error probability < 2^-50.
	return __gmpz_probab_prime_p(a.ptr, 25);
}

/** returns floor(sqrt(a)) */
Integer isqrt(Integer a)
{
	auto r = new GmpInteger;
	__gmpz_sqrt(r.ptr, a.ptr);
	return Integer(cast(immutable)r);
}

/** greatest common divisor */
Integer gcd(Integer a, int b)
{
	auto r = new GmpInteger;
	if(b<0)
		b = -b;
	__gmpz_gcd_ui(r.ptr, a.ptr, b);
	return Integer(cast(immutable)r);
}

/** ditto */
Integer gcd(Integer a, Integer b)
{
	auto r = new GmpInteger;
	__gmpz_gcd(r.ptr, a.ptr, b.ptr);
	return Integer(cast(immutable)r);
}

/** ditto */
Integer gcd(Integer a, Integer b, Integer c)
{
	auto r = new GmpInteger;
	__gmpz_gcd(r.ptr, a.ptr, b.ptr);
	__gmpz_gcd(r.ptr, r.ptr, c.ptr);
	return Integer(cast(immutable)r);
}

/** least common multiple */
Integer lcm(Integer a, int b)
{
	auto r = new GmpInteger;
	if(b<0)
		b = -b;
	__gmpz_lcm_ui(r.ptr, a.ptr, b);
	return Integer(cast(immutable)r);
}

/** ditto */
Integer lcm(Integer a, Integer b)
{
	auto r = new GmpInteger;
	__gmpz_lcm(r.ptr, a.ptr, b.ptr);
	return Integer(cast(immutable)r);
}

/** ditto */
Integer lcm(Integer a, Integer b, Integer c)
{
	auto r = new GmpInteger;
	__gmpz_lcm(r.ptr, a.ptr, b.ptr);
	__gmpz_lcm(r.ptr, r.ptr, c.ptr);
	return Integer(cast(immutable)r);
}

/** calculate a^e % m */
Integer powMod(Integer a, Integer e, Integer m)
{
	auto r = new GmpInteger;
	__gmpz_powm(r.ptr, a.ptr, e.ptr, m.ptr);
	return Integer(cast(immutable)r);
}

/** ditto */
Integer powMod(int a, Integer e, Integer m)
{
	auto r = new GmpInteger;
	__gmpz_set_si(r.ptr, a);
	__gmpz_powm(r.ptr, r.ptr, e.ptr, m.ptr);
	return Integer(cast(immutable)r);
}

/** ditto */
Integer powMod(Integer a, int e, Integer m)
{
	if(e < 0)
		throw new Exception("TODO");
	auto r = new GmpInteger;
	__gmpz_powm_ui(r.ptr, a.ptr, e, m.ptr);
	return Integer(cast(immutable)r);
}

// TODO: figure out the reason for the three different routines jacobi/legendre/kronecker

int legendreSymbol(Integer a, Integer b)
{
	return __gmpz_legendre(a.ptr, b.ptr);
}

int legendreSymbol(int a, Integer b)
{
	return __gmpz_si_kronecker(a, b.ptr);
}

int legendreSymbol(Integer a, int b)
{
	return __gmpz_kronecker_si(a.ptr, b);
}

/**
 * Returns an integer x between 0 and p-1 such that x^2 ≡ n (mod p).
 * Only works for p prime. For details of the algorithm, see:
 * http://en.wikipedia.org/wiki/Tonelli-Shanks_algorithm
 */
Integer sqrtMod(Integer n, Integer p)
{
	n = n % p;

	//assert(p.isPrime);
	//assert(p == 2 || a == 0 || kroneckerSymbol(a, p) == 1);

	if(n == 0)
		return Integer(0);
	if(n == 1)
		return Integer(1);
	if(p%4 == 3)
		return n.powMod((p+1)/4, p);

	// NOTE: at this point, we already have p >= 5

	// rewrite p-1 = s*2^e
	int s = cast(int)__gmpz_scan1(p.ptr, 1);
	auto q = p >> s; // this will round down

	// find a non-square mod p
	int z = 2;
	while(legendreSymbol(z, p) != -1)
		z += 1;
	auto c = z.powMod(q, p);

	// invariant: x^2 = ab (mod p)

	auto r = n.powMod((q+1)/2, p);
	auto t = n.powMod(q, p);
	auto m = s;

	while(true)
	{
		if(t == 1)
			return r;

		int i = 0;
		auto t2 = t;
		while(t2 != 1)
		{
			++i;
			t2 = t2.powMod(2, p);
		}

		auto b = c;
		for(int j = 0; j < m-i-1; ++j)
			b = (b*b)%p;
		r = b*r;
		c = b*b;
		t = c*t;
		m = i;
	}
}
