module math.integer;

private import std.string : toStringz;
private import std.format;
private import std.conv : to;
private import std.typecons;
private import std.exception : assumeUnique;
private import std.random : unpredictableSeed;
private import std.algorithm;
private import jive.array;
private import math.gmp;
private import math.numtheory;

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

	this(immutable GmpInteger z) pure nothrow
	{
		this.z = z;
	}

	/** constructor for given value */
	this(int v) pure nothrow
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
	this(double v) pure nothrow
	{
		auto r = new GmpInteger;
		__gmpz_set_d(r.ptr, v);
		this(cast(immutable)r);
	}

	/** ditto */
	this(string v) pure
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
	string toString() const pure nothrow @property
	{
		auto buflen = __gmpz_sizeinbase(ptr, 10)+2;	// one for sign, one for \0
		auto buf = new char[buflen];
		return to!string(__gmpz_get_str(buf.ptr, 10, ptr));
	}

	/** returns value as int, asserting it is small enough */
	int opCast(T)() const pure
		if(is(T == int))
	{
		assert(int.min <= this && this <= int.max);
		return cast(int)__gmpz_get_si(ptr);
	}

	double opCast(T)() const pure
		if(is(T == double))
	{
		return __gmpz_get_d(ptr);
	}

	/** test if bit i is set */
	bool opIndex(size_t i) const pure nothrow
	{
		return __gmpz_tstbit(ptr, i) != 0;
	}


	//////////////////////////////////////////////////////////////////////
	/// size metric and comparisons
	//////////////////////////////////////////////////////////////////////

	bool isNan() const pure nothrow @property
	{
		return z is null;
	}

	bool opEquals(int b) const pure nothrow
	{
		return __gmpz_cmp_si(ptr, b) == 0;
	}

	bool opEquals(Integer b) const pure nothrow
	{
		return __gmpz_cmp(ptr, b.ptr) == 0;
	}

	int opCmp(int b) const pure nothrow
	{
		return __gmpz_cmp_si(ptr, b);
	}

	int opCmp(Integer b) const pure nothrow
	{
		return __gmpz_cmp(ptr, b.ptr);
	}

	/** returns -1 / 0 / +1, possibly faster than actual compare */
	int sign() const @property pure nothrow
	{
		return z.z._mp_size < 0 ? -1 : z.z._mp_size > 0;
	}

	/** number of bits */
	size_t length() const pure @property
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

	Integer opUnary(string op)() const pure nothrow
		if(op == "-")
	{
		auto r = new GmpInteger;
		__gmpz_neg(r.ptr, ptr);
		return Integer(cast(immutable)r);
	}

	Integer opBinary(string op)(int b) const pure
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

	int opBinary(string op)(int b) const pure
		if(op == "%")
	{
		if(b <= 0)
			throw new Exception("only positive modulus is supported");
		scope r = new GmpInteger;
		return cast(int)__gmpz_fdiv_r_ui(r.ptr, ptr, b);
	}

	Integer opBinaryRight(string op)(int a) const pure
		if(op == "+" || op == "*")
	{
		return opBinary!op(a); // commutative operators can simply be forwarded
	}

	Integer opBinaryRight(string op)(int a) const pure
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

	Integer opBinary(string op)(Integer b) const pure
		if(op != "^^")
	{
		auto r = new GmpInteger;

		     static if(op == "+") __gmpz_add(r.ptr, ptr, b.ptr);
		else static if(op == "-") __gmpz_sub(r.ptr, ptr, b.ptr);
		else static if(op == "*") __gmpz_mul(r.ptr, ptr, b.ptr);
		else static if(op == "/") __gmpz_fdiv_q(r.ptr, ptr, b.ptr);
		else static if(op == "%") __gmpz_fdiv_r(r.ptr, ptr, b.ptr);
		else static if(op == "^^") return this^^cast(int)b;
		else static assert(false, "binary '"~op~"' is not defined");

		return Integer(cast(immutable)r);
	}

	Integer opBinary(string op)(Integer b) const pure
		if(op == "^^")
	{
		return this^^cast(int)b;
	}

	/** returns this/b. Faster, but only works if the division is exact (i.e. no rounding) */
	Integer divExact(Integer b) const pure
	{
		auto r = new GmpInteger;
		__gmpz_divexact(r.ptr, ptr, b.ptr);
		return Integer(cast(immutable)r);
	}


	//////////////////////////////////////////////////////////////////////
	/// internals
	//////////////////////////////////////////////////////////////////////

	Rebindable!(immutable(GmpInteger)) z;

	immutable(mpz_t)* ptr() const pure nothrow @property
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
Integer inverseMod(Integer a, Integer m) pure
{
	auto r = new GmpInteger;
	__gmpz_invert(r.ptr, a.ptr, m.ptr);
	return Integer(cast(immutable)r);
}

/**
 * returns true if a is a square number (this includes 0 and 1)
 */
bool isPerfectSquare(Integer a) pure
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

/**
 * checks wether 2^n-1 is a (Mersenne) prime.
 * Implemented using Lucas-Lehmer test. A simple table would be faster (and
 * not too big) for realistic inputs, but that's not the point here.
 */
bool isMersennePrime(int n) pure
{
	assert(n >= 0);
	if(n == 2)
		return true;
	if(!math.numtheory.isPrime(n))
		return false;

	// Lucas-Lehmer only works for odd primes n
	Integer p = Integer(2)^^n-1;
	Integer s = 4;
	for(int i = 0; i < n-2; ++i)
		s = (s*s-2)%p;
	return s == 0;
}

/**
 * checks wether 2^2^n+1 is a (Fermat) prime.
 * Implemented using Pepin's test. In practice, it returns
 * true for 0 <= n <= 4, and false (or timeout) for larger inputs.
 */
bool isFermatPrime(int n) pure
{
	if(n == 0)
		return true;
	Integer p = Integer(2)^^(Integer(2)^^n)+1;
	return powMod(3, (p-1)/2, p) == p-1;
}

/** returns floor(sqrt(a)) */
Integer isqrt(Integer a) pure
{
	auto r = new GmpInteger;
	__gmpz_sqrt(r.ptr, a.ptr);
	return Integer(cast(immutable)r);
}

/** greatest common divisor */
Integer gcd(Integer a, int b) pure
{
	auto r = new GmpInteger;
	if(b<0)
		b = -b;
	__gmpz_gcd_ui(r.ptr, a.ptr, b);
	return Integer(cast(immutable)r);
}

/** ditto */
Integer gcd(Integer a, Integer b) pure
{
	auto r = new GmpInteger;
	__gmpz_gcd(r.ptr, a.ptr, b.ptr);
	return Integer(cast(immutable)r);
}

/** ditto */
Integer gcd(Integer a, Integer b, Integer c) pure
{
	auto r = new GmpInteger;
	__gmpz_gcd(r.ptr, a.ptr, b.ptr);
	__gmpz_gcd(r.ptr, r.ptr, c.ptr);
	return Integer(cast(immutable)r);
}

/** least common multiple */
Integer lcm(Integer a, int b) pure
{
	auto r = new GmpInteger;
	if(b<0)
		b = -b;
	__gmpz_lcm_ui(r.ptr, a.ptr, b);
	return Integer(cast(immutable)r);
}

/** ditto */
Integer lcm(Integer a, Integer b) pure
{
	auto r = new GmpInteger;
	__gmpz_lcm(r.ptr, a.ptr, b.ptr);
	return Integer(cast(immutable)r);
}

/** ditto */
Integer lcm(Integer a, Integer b, Integer c) pure
{
	auto r = new GmpInteger;
	__gmpz_lcm(r.ptr, a.ptr, b.ptr);
	__gmpz_lcm(r.ptr, r.ptr, c.ptr);
	return Integer(cast(immutable)r);
}

/** calculate a^e % m */
Integer powMod(Integer a, Integer e, Integer m) pure
{
	auto r = new GmpInteger;
	__gmpz_powm(r.ptr, a.ptr, e.ptr, m.ptr);
	return Integer(cast(immutable)r);
}

/** ditto */
Integer powMod(int a, Integer e, Integer m) pure
{
	auto r = new GmpInteger;
	__gmpz_set_si(r.ptr, a);
	__gmpz_powm(r.ptr, r.ptr, e.ptr, m.ptr);
	return Integer(cast(immutable)r);
}

/** ditto */
Integer powMod(Integer a, int e, Integer m) pure
{
	if(e < 0)
		throw new Exception("TODO");
	auto r = new GmpInteger;
	__gmpz_powm_ui(r.ptr, a.ptr, e, m.ptr);
	return Integer(cast(immutable)r);
}

// TODO: figure out the reason for the three different routines jacobi/legendre/kronecker

int legendreSymbol(Integer a, Integer b) pure
{
	return __gmpz_legendre(a.ptr, b.ptr);
}

int legendreSymbol(int a, Integer b) pure
{
	return __gmpz_si_kronecker(a, b.ptr);
}

int legendreSymbol(Integer a, int b) pure
{
	return __gmpz_kronecker_si(a.ptr, b);
}

/**
 * Returns an integer x between 0 and p-1 such that x^2 ≡ n (mod p).
 * Only works for p prime. For details of the algorithm, see:
 * http://en.wikipedia.org/wiki/Tonelli-Shanks_algorithm
 */
Integer sqrtMod(Integer n, Integer p) pure
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

/**
 * convenience wrapper around Array!(Tuple!(Integer,int)) for integer factorization
 */
struct Factorization
{
	Array!(Tuple!(Integer,int)) factors;
	alias factors this;

	/**
	 * puts the product in human readable form
	 */
	string toString() const @property
	{
		if(factors.empty)
			return "1";
		string r;
		foreach(f; factors)
		{
			if(f[1] == 1)
				r ~= format("(%s)", f[0]);
			else
				r ~= format("(%s)^%s", f[0], f[1]);
		}

		return r;
	}

	/**
	 * sort factors and collect duplicates
	 */
	void normalize()
	{
		sort!"a[0] < b[0]"(factors[]);
		foreach(i, f, ref bool rem; &factors.prune)
			if(i+1 < factors.length && f[0] == factors[i+1][0])
			{
				factors[i+1][1] += f[1];
				rem = true;
			}
	}

	/**
	 * multiply the factorization again. Mostly for checking.
	 */
	Integer multiply()const @property
	{
		Integer r = 1;
		foreach(f; factors)
			for(int i = 0; i < f[1]; ++i)
				r = r*f[0];
		return r;
	}

	/**
	 * checks if all factors are prime
	 */
	bool allPrime()const @property
	{
		foreach(f; factors)
			if(!isPrime(f[0]))
				return false;
		return true;
	}
}

/**
 * factorize a number using pollard rho
 */
Factorization factor(Integer n)
{
	assert(n > 0);

	Factorization f;
	if(n == 1)
		return f;
	f.pushBack(tuple(n,1));

	for(int i = 0; i < f.length; ++i)
	{
		while(!isPrime(f[i][0]))
		{
			Integer d = f[i][0];
			for(int c = 1; d == f[i][0]; ++c)
				d = findFactor(f[i][0], Integer(0), c);

			f[i][0] = f[i][0].divExact(d);
			f.pushBack(tuple(d, 1));
		}
	}

	f.normalize();
	assert(f.multiply == n);
	//assert(f.allPrime);
	return f;
}

/**
 * find a factor of n using basic pollard rho
 * returns either a proper factor of n (which is not necessarily prime),
 * or n if none was found. in the latter case try using a different value for c
 */
Integer findFactor(Integer n, Integer x0, int c) pure
{
	assert(n > 0);
	assert(0 < c && c < n);

	auto x = new GmpInteger;
	auto y = new GmpInteger;
	auto d = new GmpInteger;
	__gmpz_set(x.ptr, x0.ptr);
	long runLength = 1;

	while(true)
	{
		__gmpz_set(y.ptr, x.ptr);

		for(long i = 0; i < runLength; ++i)
		{
			// x = (x * x + c) % n
			__gmpz_mul(x.ptr, x.ptr, x.ptr);
			__gmpz_add_ui(x.ptr, x.ptr, c);
			__gmpz_fdiv_r(x.ptr, x.ptr, n.ptr);

			// d = gcd(x-y, n)
			__gmpz_sub(d.ptr, x.ptr, y.ptr);
			__gmpz_gcd(d.ptr, d.ptr, n.ptr); // ignores sign of d

			if(__gmpz_cmp_ui(d.ptr, 1) != 0)
			{
				delete x;
				delete y;
				return Integer(cast(immutable)d);
			}
		}

		runLength *= 2;
	}
}
