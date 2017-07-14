module math.numtheory;

/**
 * Number Theory and combinatorical functions for 64-bit integers and smaller
 *
 * For similar functionality using arbitrary large numbers, see math.integer.
 * Also, there is no special consideration for architectures without a native
 * 64-bit integer type, so performance might be bad on such CPUs, even in
 * cases where 32-bit computions would suffice.
 */

 /+
TODO (or decide not to, if not worth it)
- dont use the slow overflow-safe modular functions when modulus is small
- prime number generation of segments not starting at 2
- "BitArray isPrime" is smaller than "long[] primes" for realisitc sizes. Should that be used more?
- skip even numbers in the tables of MultuplicativeFunction's. Saves half the memory and should be trivial to recompute on the fly.
+/

private import core.bitop : bsf;
private import std.math : log, sqrt, cbrt;
private import std.typecons;
private import std.algorithm;
private import std.range;
private import std.exception;
private import std.format;
private import std.functional : unaryFun, binaryFun;
private import jive.array;
private import jive.bitarray;
private import math.numberfield : Quadratic;


//////////////////////////////////////////////////////////////////////
/// modular arithmetic
//////////////////////////////////////////////////////////////////////

/**
 * calculate (a + b) % m
 * (without the overflow problems of the naive expression)
 * conditions: m > 0 and 0 <= a,b,result < m
 */
long addmod(long a, long b, long m) pure nothrow
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if(b < m-a)
		return a+b;
	else
		return a+b-m;
}

/**
 * calculate (a - b) % m
 * (without the overflow problems of the naive expression)
 * conditions: m > 0 and 0 <= a,b,result < m
 */
long submod(long a, long b, long m) pure nothrow
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if(a >= b)
		return a-b;
	else
		return a-b+m;
}

/**
 * calculate (-a) % m
 * conditions: m > 0 and 0 <= a,result < m
 */
long negmod(long a, long m) pure nothrow
{
	assert(m > 0);
	assert(0 <= a && a < m);

	if(a == 0)
		return 0;
	else
		return m - a;
}

/**
 * calculate (a * b) % m
 * (without the overflow problems of the naive expression)
 * conditions: m > 0 and 0 <= a,b,result < m
 */
long mulmod(long a, long b, long m) pure nothrow
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(0 <= b && b < m);

	if(a > b)
		swap(a, b);

	long r = 0;
	for (; a != 0; a >>= 1, b = addmod(b, b, m))
		if(a & 1)
			r = addmod(r, b, m);
	return r;
}

/**
 * calculate (a ^ b) % m using binary exponentiation
 * uses modular inverse for negative powers
 * conditions: m > 0 and 0 <= a,result < mod
 * triggers division by zero if b < 0 and gcd(x, m) != 1
 */
long powmod(long a, long b, long m) pure nothrow
{
	assert(m > 0);
	assert(0 <= a && a < m);

	if(b < 0)
	{
		a = invmod(a, m);
		b = -b;
	}

	long r = 1;
	if(m <= int.max) // no risk for overflow -> use faster direct expressions
	{
		for(; b != 0; b >>= 1, a = a*a%m)
			if(b & 1)
				r = r*a%m;
	}
	else
	{
		for(; b != 0; b >>= 1, a = mulmod(a, a, m))
			if(b & 1)
				r = mulmod(r, a, m);
	}
	return r;
}

/**
 * calculate the modular inverse a^-1 % m
 * conditions: m > 0 and 0 <= a,result < m and gcd(a, m) = 1
 */
long invmod(long a, long m) pure nothrow
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(gcd(a,m) == 1);

	long a0 = m;
	long a1 = a;
	long b0 = 0;
	long b1 = 1;

	while(a1 > 1)
	{
		long q = a0 / a1;
		long a2 = a0 - q*a1;
		long b2 = b0 - q*b1;

		a0 = a1;
		a1 = a2;
		b0 = b1;
		b1 = b2;
	}

	if(b1 < 0)
		b1 += m;
	assert(0 <= b1 && b1 < m);
	//assert(mulmod(b1, a, m) == 1);
	return b1;
}

unittest
{
	for(long m = 1; m <= 20; ++m)
		for(long x = 1; x < m; ++x)
			if(gcd(x,m) == 1)
				assert(x * invmod(x,m) % m == 1);
}

/**
 * congruence classes of the form [x] = x + nZ
 * (with 0 <= x < n by convention)
 */
struct IntMod
{
	long x;
	long n;

	this(long x, long n) pure nothrow
	{
		assert(n > 0);
		assert(0 <= x && x < n);
		this.x = x;
		this.n = n;
	}

	static IntMod make(long x, long n) pure nothrow
	{
		assert(n > 0);
		return IntMod((x%n+n)%n,n);
	}

	/** (additive) inverse */
	IntMod opUnary(string op)() const pure nothrow
		if(op == "-")
	{
		if(x == 0)
			return this;
		else
			return IntMod(n-x, n);
	}

	/** (multiplicative) inverse */
	IntMod inverse() const @property pure nothrow
	{
		return IntMod(invmod(x, n), n);
	}

	IntMod opBinary(string op)(long b) const pure nothrow
		if(op == "+" || op == "-" || op == "*" || op == "/")
	{
		return opBinary!op(make(b, n));
	}

	IntMod opBinary(string op)(IntMod b) const pure nothrow
		if(op == "+" || op == "-" || op == "*" || op == "/")
	{
		assert(n == b.n);

		final switch(op)
		{
			case "+": return IntMod(addmod(x, b.x, n), n);
			case "-": return IntMod(submod(x, b.x, n), n);
			case "*": return IntMod(mulmod(x, b.x, n), n);
			case "/": return IntMod(mulmod(x, invmod(b.x, n), n), n);
		}
	}

	IntMod opBinary(string op)(long e) const pure nothrow
		if(op == "^^")
	{
		return IntMod(powmod(x, e, n), n);
	}

	bool compatible(IntMod b) const pure nothrow
	{
		return (x - b.x) % gcd(n, b.n) == 0;
	}

	/** chinese remainder */
	IntMod opBinary(string op)(IntMod b) const pure
		if(op == "&")
	{
		// TODO: this is not overflow-safe
		long y, z;
		long d = euclid(n, b.n, y, z);
		if((x-b.x) % d != 0)
			throw new Exception(format("no solution: (%s %% %s) & (%s %% %s)", x, n, b.x, b.n));
		return make(x - (x - b.x)/d*y*n, n*b.n/d);
	}

	bool opEquals(int b) const pure nothrow
	{
		return opEquals(make(b, n));
	}

	bool opEquals(IntMod b) const pure nothrow
	{
		assert(n == b.n);
		return x == b.x;
	}

	string toString() const pure @property
	{
		return format("[%s]", x);
	}

	byte jacobi() const pure nothrow
	{
		assert(n % 2 != 0);
		return .jacobi(x, n);
	}
}

/**
 * modular square-root implemented using Cipolla's algorithm.
 * Only works for prime fields.
 */
IntMod sqrt(IntMod a)
{
	assert(isPrime(a.n));

	if(a.n == 2)
		return a;

	assert(jacobi(a.x, a.n) == 1);

	if(a.n % 4 == 3)
		return a ^^ ((a.n+1)/4);

	auto z = IntMod(0, a.n);
	while((z*z-a).jacobi() != -1)
		++z.x;
	auto omega = z*z-a;

	auto b = Quadratic!IntMod(z, IntMod(1,a.n), omega);
	b = b ^^ ((a.n+1)/2);
	assert(b.b == 0);
	return b.a;
}

unittest
{
	auto x = IntMod(3,17);
	auto y = IntMod(5,17);
	assert(x*y/y == 3);
}


//////////////////////////////////////////////////////////////////////
/// prime number generation and testing
//////////////////////////////////////////////////////////////////////

/** generate all prime numbers <= n using Eratosthenes */
Array!long calculatePrimes(long n)
{
	if(n < 0)
		n = 0;

	// excluding 2 and 3 as special cases, all primes have the form 6*k +- 1
	auto b5 = BitArray(n/6+1); // b5[k] represents 6*k+5
	auto b7 = BitArray(n/6+1); // b7[k] represents 6*k+7

	// limit for relevant prime divisors
	long limit = sqrti(max(6*(b7.length-1)+7, 6*(b5.length-1)+5));

	// mark all non-primes
	// NOTE: in older versions of DMD, making k an ulong triggers a strange bug:
	// issues.dlang.org/show_bug.cgi?id=13023
	for(long k = 0; k < n/6; ++k)
	{
		if(!b5[cast(size_t)k])
		{
			long p = 6*k+5;

			if(p > limit)
				break;

			for(long s = (p*(p+2)-5)/6; s < b5.length; s += p)
				b5[cast(size_t)s] = true;

			for(long s = (p*p-7)/6; s < b7.length; s += p)
				b7[cast(size_t)s] = true;
		}

		if(!b7[cast(size_t)k])
		{
			long p = 6*k+7;

			if(p > limit)
				break;

			for(long s = (p*(p+4)-5)/6; s < b5.length; s += p)
				b5[cast(size_t)s] = true;

			for(long s = (p*p-7)/6; s < b7.length; s += p)
				b7[cast(size_t)s] = true;
		}
	}

	// collect primes into array
	Array!long primes;
	primes.reserve(b5.count(false) + b7.count(false) + 2);
	primes.pushBack(2);
	primes.pushBack(3);
	for(long k = 0; k < n/6+1; ++k)
	{
		if(!b5[cast(size_t)k])
			primes.pushBack(6L*k+5);
		if(!b7[cast(size_t)k])
			primes.pushBack(6L*k+7);
	}

	// now we might have computed slightly more primes than requested,
	// so we remove them again (simpler than doing it right in the beginning)
	while(!primes.empty && primes[$-1] > n)
		primes.popBack;

	return primes;
}

/**
 * returns all primes in [a..b] using cached results from calculatePrimes
 */
immutable(long)[] primes(long a, long b)
{
	static immutable(long)[] cache;
	static long limit = 1;

	if(b > limit)
	{
		limit = max(b, 2*limit);
		cache = assumeUnique(calculatePrimes(limit).release);
	}

	return cache[].assumeSorted.upperBound(a-1).lowerBound(b+1).release;
}

/**
 * alias for primes(0, n)
 */
immutable(long)[] primes(long n)
{
	return primes(0, n);
}

/**
 * tests wether n is a strong probable prime to base a
 * conditions a >= 0 and n >= 3 odd
 * note that for a % n = -1,0,1 this returns always true
 */
bool isSPRP(long a, long n) pure nothrow
{
	assert(a >= 0);
	assert(n > 1 && n % 2 == 1);

	// the original definition assumes a < n-1, we use a straight-forward generalization
	a %= n;
	if(a == 0 || a == 1 || a == n-1)
		return true;

	long s = bsf(n-1);
	long d = (n-1) >> s;

	// now it is n-1 = d*2^s, t odd

	a = powmod(a, d, n);
	if(a == 1 || a == n-1)
		return true;

	while(--s)
	{
		a = mulmod(a, a, n);
		if(a == n-1)
			return true;
	}

	return false;
}

/** Test if n is prime. */
bool isPrime(long n) pure nothrow
{
	if(n < 53)	// includes the trivial n < 0 case
		switch(n)
		{
			case 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47:
				return true;
			default:
				return false;
		}

	foreach(p; [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47])
		if(n % p == 0)
			return false;

	if(n < 53*53)
		return true;

	return isPrimeMillerRabin(n);
}

/**
 * Deterministic version of Miller-Rabin primality test.
 * Use this after small factors have been tested.
 * see http://priv.ckp.pl/wizykowski/sprp.pdf for implementation details and
 * http://miller-rabin.appspot.com for the actual values used therein.
 */
bool isPrimeMillerRabin(long n) pure nothrow
{
	if(n < 291_831_L)
		return isSPRP(  126_401_071_349_994_536_L, n);

	if(n < 1_050_535_501_L)
		return isSPRP(          336_781_006_125_L, n)
		    && isSPRP(    9_639_812_373_923_155_L, n);

	if(n < 273_919_523_041_L)
		return isSPRP(                       15_L, n)
		    && isSPRP(            7_363_882_082_L, n)
		    && isSPRP(      992_620_450_144_556_L, n);

	if(n < 47_636_622_961_201_L)
		return isSPRP(                        2_L, n)
		    && isSPRP(                2_570_940_L, n)
		    && isSPRP(              211_991_001_L, n)
		    && isSPRP(            3_749_873_356_L, n);

	if(n < 3_770_579_582_154_547_L)
		return isSPRP(                        2_L, n)
		    && isSPRP(                  880_937_L, n)
			&& isSPRP(                2_570_940_L, n)
		    && isSPRP(              610_386_380_L, n)
		    && isSPRP(            4_130_785_767_L, n);

	if(n < 585_226_005_592_931_977_L)
		return isSPRP(                        2_L, n)
		    && isSPRP(      123_635_709_730_000_L, n)
		    && isSPRP(    9_233_062_284_813_009_L, n)
		    && isSPRP(   43_835_965_440_333_360_L, n)
		    && isSPRP(  761_179_012_939_631_437_L, n)
		    && isSPRP(1_263_739_024_124_850_375_L, n);

	if(true) // should work for all n < 2^64
		return isSPRP(                        2_L, n)
		    && isSPRP(                      325_L, n)
		    && isSPRP(                    9_375_L, n)
		    && isSPRP(                   28_178_L, n)
		    && isSPRP(                  450_775_L, n)
		    && isSPRP(                9_780_504_L, n)
		    && isSPRP(            1_795_265_022_L, n);
}

/**
 *  return lowest prime > n.
 *  not suitable for many continous primes, only for individual values.
 */
long nextPrime(long n) pure nothrow
{
	assert(n < 9223372036854775783L);	// largest 63 bit prime
	if(n < 2)
		return 2;

	n = (n+1) | 1;	// round to next odd prime
	while(!isPrime(n))
		n += 2;
	return n;
}

unittest
{
	assert(primes(3,19) == [3,5,7,11,13,17,19]);
	assert(isPrime(1000000007));
	assert(isPrime(9223372036854775783L)); // largest 63 bit prime
	assert(!isPrime(1000000007L*1000000009L));
	assert(equal(map!nextPrime(iota(0,21)), [2,2,3,5,5,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,23][]));

	//TODO: implement countPrimes which is actually fast
	//for(int n = 0; n < 1000; ++n)
	//	assert(countPrimes(n) == primes(n).length);
}

/** generate all Hilbert-primes <= n using a simple sieve (OEIS A057948) */
Array!long calculateHilbertPrimes(long n)
{
	if(n < 5)
		return Array!long.init;

	// initialize sieve array
	auto b = BitArray((n-1)/4+1); // b[k] represents 4*k+1
	b[0] = true;	// 1 is not a (Hilbert-)prime

	// cross out all non-primes
	for(long k = 0; k < b.length; ++k)
		if(!b[k])
		{
			long p = 4*k+1;
			long s = (p*p-1)/4;	// first index to cross out
			if(s >= b.length)
				break;
			for(; s < b.length; s += p)
				b[s] = true;
		}

	// collect primes into array
	Array!long primes;
	primes.reserve(b.count(false));
	for(long k = 0; k < b.length; ++k)
		if(!b[k])
			primes.pushBack(4*k+1);
	return primes;
}

unittest
{
	assert(equal(calculateHilbertPrimes(200)[], [5,9,13,17,21,29,33,37,41,49,53,57,61,69,73,77,89,93,97,101,109,113,121,129,133,137,141,149,157,161,173,177,181,193,197][]));
}


//////////////////////////////////////////////////////////////////////
/// integer factorization
//////////////////////////////////////////////////////////////////////

/**
 * Find a prime factor of a positive number. For small numbers we use a table
 * (automatically extended up to tableLimit). For larger numbers we do some
 * trial division and pollard rho alrogithm for the serious cases. Typically,
 * the smallest prime factor is returned, but that is not guaranteed. In
 * particular for large numbers with prime factors of similar size, a bigger
 * one might be found first.
 */
struct primeFactor
{
	static:

	// (smallest) prime factor of odd numbers. zero for prime numbers.
	private Array!ushort table;

	// limit up to which auto-extend the table
	enum long tableLimit = 1L<<28; // that is a 256MiB table at full size

	long currLimit() @property
	{
		return 2*table.length;
	}

	void makeTable(long limit)
	{
		if(limit <= currLimit)
			return;

		// change table to a larger integer type or do some tighter packing to increase this limit
		if(sqrti(limit) > ushort.max)
			assert(false);

		table.assign((limit+1)/2, 0);
		foreach(p; primes(3, sqrti(limit)))
			for(long k = p*p/2; k < table.length; k += p)
				if(table[k] == 0)
					table[k] = cast(ushort)p;
	}

	long opCall(long n)
	{
		assert(n > 1);

		// even numbers are easy
		if(n % 2 == 0)
			return 2;

		// already computed values
		if(n <= currLimit)
			if(table[n/2] == 0)
				return n;
			else
				return table[n/2];

		// auto extend table if n is not too large
		if(n <= tableLimit)
		{
			makeTable(min(max(n, 2*currLimit), tableLimit));
			if(table[n/2] == 0)
				return n;
			else
				return table[n/2];
		}

		// a little trial division
		foreach(long p; [3,5,7,11,13,17,19,23,29,31,37,41,43,47])
			if(n % p == 0)
				return p;

		// Actual factorizing algorithm.
		// by now we now that n is large and does not contain very small
		// factors. So we use the low-overhead version of isPrime
		long c = 1;
		while(!isPrimeMillerRabin(n))
			n = pollardRho(n, c++);
		return n;
	}
}

/**
 * range over all prime factors (with multiplicity) of a number.
 * See primeFactor(...) for algorithmic details.
 *  - Using this range does not invoke any memory allocation. So factoring a
 *    lot of small numbers should be quite fast
 *  - order of prime factors is typically ascending, but that is not guaranteed
 */
struct factors
{
	private long p = 1; // current prime factor
	private int r = 0; // multiplicity of current factor
	private long n; // remaining stuff to be factored (f already removed)

	@disable this();

	this(long n)
	{
		assert(n > 0);
		this.n = n;
		popFront();
	}

	bool empty() const @property
	{
		return p == 1;
	}

	Tuple!(long,int) front() const @property
	{
		return Tuple!(long,int)(p,r);
	}

	void popFront()
	{
		if(n == 1)
		{
			p = 1;
			r = 0;
			return;
		}

		p = primeFactor(n);
		assert(n % p == 0);
		r = 1;
		n /= p;
		while(n % p == 0)
		{
			++r;
			n /= p;
		}
	}
}

/**
 * prime factorization of a positive number. Actually just a thin wrapper
 * around Array!(Tuple!(long,int))(factors(n)). In most cases you should use
 * factors(n) directly. But this wrapper provides:
 *   - prime factors are guaranteed to be in ascending order
 *   - there is a nice (human readable) toString method
 */
struct Factorization
{
	Array!(Tuple!(long,int)) fs;
	alias fs this;

	this(long n)
	{
		assert(n > 0);

		foreach(f; factors(n))
			fs.pushBack(f);

		normalize();
		assert(multiply == n);
		//assert(allPrime);
	}

	/**
	 * puts the product in human readable form
	 */
	string toString() const @property
	{
		if(this.empty)
			return "1";
		string r;
		foreach(f; this)
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
		sort!"a[0] < b[0]"(this[]);
		foreach(i, f, ref bool rem; &this.prune)
			if(i+1 < this.length && f[0] == this[i+1][0])
			{
				this[i+1][1] += f[1];
				rem = true;
			}
	}

	/**
	 * multiply the factorization again. Mostly for checking.
	 */
	long multiply() pure nothrow const @property
	{
		long r = 1;
		foreach(f; this)
			for(int i = 0; i < f[1]; ++i)
				r *= f[0];
		return r;
	}

	/**
	 * checks if all factors are prime
	 */
	bool allPrime() pure nothrow const @property
	{
		foreach(f; this)
			if(!isPrime(f[0]))
				return false;
		return true;
	}
}

/**
 * find a factor of n using a simple Pollard rho algorithm.
 * returns either a proper factor of n (which is not necessarily prime), or n
 * if none was found. in the latter case try using a different value for c.
 * make sure n is composite before calling this, no check is performed (result
 * is correct but slow for prime n). For composite n the heuristic runtime is
 * O(n^(1/4)), but depends pseudo-randomly ony the parameter c.
 */
long pollardRho(long n, long c) pure nothrow
{
	assert(n > 0);
	assert(0 < c && c < n);

	// arbitrary start value. Could be randomized, the parameter c seems to be sufficient variation for now.
	long x = 0;
	long runLength = 1;

	while(true)
	{
		long y = x;

		for(long i = 0; i < runLength; ++i)
		{
			x = mulmod(x, x, n);
			x = addmod(x, c, n);
			long d = gcd(x-y, n);

			if(d != 1)
				return d;
		}

		runLength *= 2;
	}
}

unittest
{
	assert(Factorization(2*2*2*2*3*5*5*7)[] == [tuple(2,4),tuple(3,1),tuple(5,2),tuple(7,1)]);
	assert(Factorization(1000000007L*1000000009L)[] == [tuple(1000000007L,1), tuple(1000000009L,1)]);
}

/**
 * Range over all (positive) divisors of a positive number. Divisors are not
 * strictly ordered, but if x divides y, x will appear before y.
 * Not a real range with front/popFront, but only opApply. But no heap
 * allocations at all, which seems to make it slightly faster than putting
 * everything into an array and iterating over that.
 */
struct divisors
{
	long n;

	@disable this();

	this(long n)
	{
		assert(n > 0);
		this.n = n;
	}

	private static int apply(in int delegate(long) dg, long n, factors fs)
	{
		if(fs.empty)
			return dg(n);

		auto p = fs.front;
		fs.popFront;

		int r;
		for(; p[1] >= 0; n *= p[0], p[1]--)
			if((r = apply(dg,n,fs)) != 0)
				return r;
		return 0;
	}

	int opApply(in int delegate(long) dg) const
	{
		return apply(dg, 1, factors(n));
	}
}


//////////////////////////////////////////////////////////////////////
/// combinatoric functions
//////////////////////////////////////////////////////////////////////

/**
 * caclulate the factorial
 * n! = n * (n-1) * ... * 2 * 1
 */
long factorial(long n) pure nothrow
{
	assert(0 <= n && n <= 20);

	long r = 1;
	for(long i = 2; i <= n; ++i)
		r *= i;
	return r;
}

/**
 * calculate a binomial coefficient
 * (n over k) = n! / (k! * (n-k)!)
 */
long binomial(long n, long k) pure nothrow
{
	assert(n >= 0 && k >= 0);

	if(k > n)
		return 0;

	// use symmetry to minimize work
	if(n - k < k)
		k = n - k;

	long r = 1;
	for(long i = 1; i <= k; ++i)
	{
		assert(r <= long.max/(n+1-i));
		r *= (n+1-i);
		r /= i;
	}
	return r;
}

/**
 * calculate a binomial coefficient modulo a prime using Lucas's theorem
 */
long binomialMod(long n, long k, long p)
{
	assert(isPrime(p));
	assert(n >= 0 && k >= 0);
	if(k > n)
		return 0;

	long r = 1;
	long s = 1;

	while(k > 0)
	{
		long a = n % p;
		long b = k % p;
		n /= p;
		k /= p;

		if(b > a)
			return 0;

		if(a - b < b)
			b = a - b;

		for(long i = 1; i <= b; ++i)
		{
			r = mulmod(r, a + 1 - i, p);
			s = mulmod(s, i, p);
		}
	}

	return mulmod(r, invmod(s, p), p);
}

/**
 * calculate a multinomial coefficient
 * (∑a_i)! / ∏a_i!
 * for maximum efficiency, put the largest a_i into a[0]
 */
long multinomial(const(long)[] a...)
{
	if(a.length == 0)
		return 1;
	long r = 1;
	long k = a[0]+1;
	assert(a[0] >= 0);
	foreach(x; a[1..$])
	{
		assert(x >= 0);
		for(long i = 1; i <= x; ++i,++k)
		{
			if(r > long.max/k)
				throw new Exception("overflow in multinomial(...)");
			r *= k;
			r /= i;
		}
	}
	return r;
}


/**
 *  compute the n'th fibonacci number
 *  f(0) = 0, f(1) = 1, f(n+2) = f(n) + f(n+1)
 *  n must be in the range [-92,92], so that the result fits into a long
 */
long fibonacci(long n) pure nothrow
{
	assert(-92 <= n && n <= 92);

	if(n == 0)
		return 0;
	if(n < 0)
		if(n % 2 == 0)
			return -fibonacci(-n);
		else
			return fibonacci(-n);

	long a = 1, b = 0, c = 0, d = 1;
	if (n <= 0)
		return 0;
	for(n -= 1; n > 0; n /= 2)
	{
		if (n % 2 == 1)
		{
			long t = d*(b + a) + c*b;
			a = d*b + c*a;
			b = t;
		}
		long t = d*(2*c + d);
		c = c*c + d*d;
		d = t;
	}
	return a + b;
}

/**
 * fibonacci(n) % m
 */
long fibonacciMod(long n, long m) pure nothrow
{
	if(m == 2)
		return n%3 == 0 ? 0 : 1;

	assert(m >= 3);

	if(n == 0)
		return 0;
	if(n < 0)
		if(n % 2 == 0)
			return negmod(fibonacciMod(-n,m), m);
		else
			return fibonacciMod(-n,m);

	long a = 1, b = 0, c = 0, d = 1;
	if (n <= 0)
		return 0;
	for(n -= 1; n > 0; n /= 2)
	{
		if (n % 2 == 1)
		{
			long t = addmod(mulmod(d, addmod(b,a,m),m), mulmod(c,b,m), m);
			a = addmod(mulmod(d,b,m), mulmod(c,a,m), m);
			b = t;
		}
		long t = mulmod(d, addmod(mulmod(2,c,m), d, m), m);
		c = addmod(mulmod(c,c,m), mulmod(d,d,m), m);
		d = t;
	}
	return addmod(a,b,m);
}

/**
 * caclulate 1^k + 2^k + .. + n^k using Faulhaber's formula
 *
 * currently only implemented for k <= 3, because the general case needs
 * Bernoulli numbers which are not whole numbers so I need to think about a
 * beautiful way to write that (without rounding or actual rational arithmetic)
 */
long powerSum(long n, long k) pure nothrow
{
	assert(0 <= n);

	switch(k)
	{
		case 0: return n;
		case 1: return n*(n+1)/2;
		case 2: return n*(n+1)*(2*n+1)/6;
		case 3: return n*n*(n+1)*(n+1)/4;
		default: assert(false);
	}
}

/**
 * powerSum(n, k) % m.
 * takes O(k^2) time, but currently only works for prime m > k.
 */
long powerSumMod(long n, long k, long m) pure nothrow
{
	long sum = 0;
	long y = 0;
	for(long i = 1; i <= k+1; ++i)
	{
		y = addmod(y, powmod(i, k, m), m);
		long a = 1;
		long b = 1;

		for(long j = 0; j <= k+1; ++j)
			if(j != i)
			{
				a = mulmod(a, submod(n%m, j, m), m);
				b = mulmod(b, submod(i, j, m), m);
			}
		sum = addmod(sum, mulmod(mulmod(a, invmod(b,m), m), y, m), m);
	}
	return sum;
}

/**
 * calculate 1 + a + a^2 + ... + a^n = (a^(n+1) - 1) / (a - 1) mod m
 * without division which might be problematic if a-1 and m are not coprime
 * a = 0 is allowed (using 0^0 = 1)
 */
long geometricMod(long a, long n, long m) pure nothrow
{
	assert(m > 0);
	assert(0 <= a && a < m);
	assert(n >= 0);

	long factor = 1;
	long sum = 0;

	while(n > 0 && a != 0)
	{
		if(n % 2 == 0)
		{
			sum = addmod(sum, mulmod(factor, powmod(a, n%m, m), m), m);
			n--;
		}

		factor = mulmod(addmod(1,a,m), factor, m);
		a = mulmod(a, a, m);

		n /= 2;
	}

	return addmod(sum, factor, m);
}

/**
 * Count numbers in the range [1,n] that are multiples of one or more a[i].
 * Implemented using inclusion-exclusion principle. Takes up to O(n) time,
 * but may be significantly faster if the a[i] are few or big. Assumes that
 * the a[i] are coprime, positive and in ascending order.
 */
long countMultiples(const(long)[] a, long n)
{
	assert(n >= 0);
	assert(a.length == 0 || a[0] > 0);
	for(size_t i = 1; i < a.length; ++i)
		assert(a[i-1] < a[i]);
	/*for(size_t i = 0; i < a.length; ++i)
		for(size_t j = i+1; j < a.length; ++j)
			assert(gcd(a[i],a[j]) == 1);*/

	static long f(const(long)[] a, long n)
	{
		long count = n;
		foreach(i, x; a)
		{
			if(x > n)
				break;
			count -= f(a[i+1..$], n/x);
		}
		return count;
	}

	return n - f(a,n);
}

unittest
{
	assert(factorial(0) == 1);
	assert(factorial(5) == 120);

	assert(binomial(0,0) == 1);
	assert(binomial(4,2) == 6);
	for(int n = 0; n < 61; ++n)
		for(int k = 0; k < 61; ++k)
			foreach(p; [2,3,5,7,11,13,17,19])
				assert(binomial(n,k)%p == binomialMod(n,k,p));
	for(int n = 1; n < 61; ++n)
		for(int k = 1; k < 61; ++k)
			assert(binomial(n,k) == binomial(n-1,k) + binomial(n-1,k-1));

	assert(fibonacci(0) == 0);
	assert(fibonacci(1) == 1);
	for(int i = -50; i < 50; ++i)
	{
		assert(fibonacci(i) == fibonacci(i-1) + fibonacci(i-2));
		for(int m = 2; m < 20; ++m)
			assert((fibonacci(i)%m+m)%m == fibonacciMod(i,m));
	}

	assert(powerSum(5, 0) == 1 + 1 + 1 + 1 + 1);
	assert(powerSum(5, 1) == 1 + 2 + 3 + 4 + 5);
	assert(powerSum(5, 2) == 1*1 + 2*2 + 3*3 + 4*4 + 5*5);
	assert(powerSum(5, 3) == 1*1*1 + 2*2*2 + 3*3*3 + 4*4*4 + 5*5*5);

	foreach(p; primes(5,100))
		foreach(k; 0..4)
			foreach(n; 1..p)
				assert(powerSum(n,k)%p == powerSumMod(n,k,p));

	for(int a = 0; a < 15; ++a)
		for(long n = 0; n < 15; ++n)
		{
			long x = geometricMod(a, n, long.max);
			long y;
			if(a == 0)
				y = 1;
			else if(a == 1)
				y = n+1;
			else
				y = (a^^(n+1)-1)/(a-1);

			assert(x == y);
		}

	assert(countMultiples([2,3,5,7], 10000) == 7715);
}


//////////////////////////////////////////////////////////////////////
/// rounded real functions
//////////////////////////////////////////////////////////////////////

/** returns floor(log_b(n)) */
int logi(long n, long b) pure nothrow
{
	assert(n > 0);
	assert(b > 1);

	int r = 0;
	while(n >= b)
	{
		++r;
		n /= b;
	}
	return r;
}

/** returns floor(sqrt(a)) */
long sqrti(long a) pure nothrow
{
	assert(a >= 0);
	auto r = cast(long)sqrt(cast(real)a);
	assert(r*r <= a && a < (r+1)*(r+1)); // not sure about rounding, so check...
	return r;
}

/** return floor(cbrt(a)) */
long cbrti(long a) /*pure*/ nothrow
{
	assert(a >= 0);
	auto r = cast(long)cbrt(cast(real)a);
	assert(r*r*r <= a && a < (r+1)*(r+1)*(r+1));
	return r;
}


//////////////////////////////////////////////////////////////////////
/// number theoretic functions
//////////////////////////////////////////////////////////////////////

/** returns largest k such that p^k divides n */
int powerOf(long n, long p) pure nothrow
{
	assert(n > 0);
	assert(p > 1);

	int r = 0;
	while(n % p == 0)
	{
		++r;
		n /= p;
	}
	return r;
}

struct MultiplicativeFunction(alias fun, alias mult = "a*b", long neutral = 1)
{
	static:

	enum tableLimit = 1L<<25; // that is a 256MiB table at full size
	alias f = binaryFun!(fun, "p", "e");
	alias mul = binaryFun!(mult, "a", "b");

	private Array!long table;

	/** make a table of all values n <= limit */
	void makeTable(long limit)
	{
		if(limit < cast(long)table.length)
			return;

		table.assign(limit+1, neutral);
		foreach(p; primes(limit))
			for(long q = p, e = 1; q <= limit; q *= p, ++e)
			{
				long a = f(p, e);
				for(long n = 1; n <= limit/q; ++n)
					if(n % p != 0)
						table[n*q] = mul(table[n*q], a);
			}
	}

	/** compute f(n) */
	long opCall(long n)
	{
		assert(n > 0);

		// already computed value
		if(n < table.length)
			return table[n];

		// automatically extend table
		if(n <= tableLimit)
		{
			makeTable(min(tableLimit, max(n, 2*table.length)));
			assert(n < table.length);
			return table[n];
		}

		// exceed table limit -> compute by explicit factoring
		long r = neutral;
		foreach(i; factors(n))
			r = mul(r, f(i[0], i[1]));
		return r;
	}
}

/** number of divisors / sigma_0 */
alias tau = MultiplicativeFunction!"e+1";

/** sum of divisors / sigma_1 */
alias sigma = MultiplicativeFunction!"(p^^(e+1) - 1) / (p-1)";

/** sum of square of divisors / sigma_2 */
alias sigma2 = MultiplicativeFunction!"(p^^(2*e+2) - 1) / (p*p-1)";

/** Euler's totient function */
alias phi = MultiplicativeFunction!"p^^(e-1)*(p-1)";

/** Carmichael function / reduced totient function */
alias carmichael = MultiplicativeFunction!("(p==2 && e > 2) ? 2L^^(e-2) : p^^(e-1)*(p-1)", lcm);

/** radical function */
alias rad = MultiplicativeFunction!"p";

/** Möbius function */
alias mu = MultiplicativeFunction!("e==1?-1:0");

/** number of distinct prime factors */
alias omega = MultiplicativeFunction!("1", "a+b", 0);

/** number of (possibly not distinct) prime factors */
alias Omega = MultiplicativeFunction!("e", "a+b", 0);

/** gcd(n,n'). It is equal to 1 if and only if n is square free */
alias gcdDerivative = MultiplicativeFunction!("e%p==0 ? p^^e : p^^(e-1)");

/** arithmetic derivative */
long derivative(long n)
{
	if(n == 0)
		return 0;
	if(n < 0)
		return -derivative(-n);

	long r = 0;
	foreach(f; factors(n))
		r += n / f[0] * f[1];
	return r;
}

unittest
{
	import std.algorithm : equal;
	import std.range : iota;

	assert(equal(map!tau(iota(1,21)), [1,2,2,3,2,4,2,4,3,4,2,6,2,4,4,5,2,6,2,6][]));
	assert(equal(map!sigma(iota(1,21)), [1,3,4,7,6,12,8,15,13,18,12,28,14,24,24,31,18,39,20,42][]));
	assert(equal(map!sigma2(iota(1,21)), [1,5,10,21,26,50,50,85,91,130,122,210,170,250,260,341,290,455,362,546][]));
	assert(equal(map!phi(iota(1,21)), [1,1,2,2,4,2,6,4,6,4,10,4,12,6,8,8,16,6,18,8][]));
	assert(equal(map!carmichael(iota(1,21)), [1,1,2,2,4,2,6,2,6,4,10,2,12,6,4,4,16,6,18,4][]));
	assert(equal(map!rad(iota(1,21)), [1,2,3,2,5,6,7,2,3,10,11,6,13,14,15,2,17,6,19,10][]));
	assert(equal(map!mu(iota(1,21)), [1,-1,-1,0,-1,1,-1,0,0,1,-1,0,-1,1,1,0,-1,0,-1,0][]));
	assert(equal(map!omega(iota(1,21)), [0,1,1,1,1,2,1,1,1,2,1,2,1,2,2,1,1,2,1,2][]));
	assert(equal(map!Omega(iota(1,21)), [0,1,1,2,1,2,1,3,2,2,1,3,1,2,2,4,1,3,1,3][]));
	assert(equal(map!gcdDerivative(iota(1,21)), [1,1,1,4,1,1,1,4,3,1,1,4,1,1,1,16,1,3,1,4][]));
	assert(equal(map!derivative(iota(1,21)), [0,1,1,4,1,5,1,12,6,7,1,16,1,9,8,32,1,21,1,24][]));
}

/** Jacobi symbol (a/n) defined for any odd integer n */
byte jacobi(long a, long n) pure nothrow
{
	assert(n > 0 &&  n % 2 == 1);

	a %= n;
	if(a < 0)
		a += n;
	if(n == 1)
		return 1;

	byte r = 1;
	while(true)
	{
		if(a == 0)
			return 0;

		if(a == 1)
			return r;

		if ((a & 1) == 0)
		{
			if(n % 8 == 3 || n % 8 == 5)
				r = -r;

			a /= 2;
		}
		else
		{
			if(a % 4 == 3 && n % 4 == 3)
				r = -r;

			swap(a, n);
			a %= n;
		}
	}
}

/**
 * calculate greatest common divisor of a and b.
 * Sign of a and b is ignored, result is always >= 0.
 * Convention: gcd(0,x) = abs(x) = gcd(x,0).
 */
long gcd(long a, long b) pure nothrow
{
	while(true)
	{
		if(a == 0)
			return b>=0?b:-b;
		b %= a;

		if(b == 0)
			return a>=0?a:-a;
		a %= b;
	}
}

/**
 * extended euclid algorithm
 * returns gcd(a,b) = x*a + y*b
 * result >= 0
 * limitations: a,b != 0
 */
long euclid(long a, long b, ref long x, ref long y) pure nothrow
{
	assert(a != 0 && b != 0);

	long a0 = a, x0 = 1, y0 = 0;
	long a1 = b, x1 = 0, y1 = 1;

	while(a1 != 0)
	{
		long q = a0/a1;

		long a2 = a0 - q*a1;
		long x2 = x0 - q*x1;
		long y2 = y0 - q*y1;

		a0 = a1, a1 = a2;
		x0 = x1, x1 = x2;
		y0 = y1, y1 = y2;
	}

	if(a0 < 0)
	{
		x = -x0;
		y = -y0;
		return -a0;
	}
	else
	{
		x = x0;
		y = y0;
		return a0;
	}
}

/**
 * calculate least common multiple of a and b.
 * Sign of a and b is ignored, result is always >= 0.
 * Convention: lcm(0,x) = 0 = lcm(x,0).
 */
long lcm(long a, long b) pure nothrow
{
	if(a == 0 || b == 0)
		return 0;
	else
		return a/gcd(a,b)*b;
}

/**
 * count numbers in the range [1,n] that are coprime to x. Generalization of
 * Euler's totient function, implemented by explicitly factorizing x.
 */
long countCoprimes(long x, long n)
{
	static Array!long ps;
	ps.resize(0);
	foreach(p; factors(x))
		ps.pushBack(p[0]);
	return n - countMultiples(ps[], n);
}

/** determine if n is a perfect square (OEIS A000290) */
bool isSquare(long n) pure nothrow
{
	if(n < 0)
		return false;
	auto r = sqrti(n);
	return r*r == n;
}

/** determine if n is a perfect cube (OEIS A000578) */
bool isCube(long n) /*pure*/ nothrow
{
	assert(n >= 0);
	auto r = cbrti(n);
	return r*r*r == n;
}

/**
 * Determine if n is square-free, i.e. not divisable by a square other than 1.
 * Takes O(n^(1/3+ϵ)) time. (OEIS A005117)
 */
bool isSquareFree(long n)
{
	assert(n > 0);
	foreach(p; primes(cbrti(n)))
	{
		if(p*p*p > n)
			break;
		if(n % p == 0)
		{
			n /= p;
			if(n % p == 0)
				return false;
		}
	}
	return n == 1 || !isSquare(n);
}

/**
 * Count square-free numbers <= n. O(n^(1/2+ϵ)). (OEIS A013928, different offset)
 */
long countSquareFree(long n)
{
	long f(long n, const long[] primes)
	{
		long r = n;
		foreach(i, p; primes)
		{
			if(n/p/p == 0)
				break;
			r -= f(n/p/p, primes[i+1..$]);
		}
		return r;
	}

	if(n < 1)
		return 0;
	return f(n, primes(sqrti(n)));
}

unittest
{
	assert(countSquareFree(100000000) == 60792694);
	assert(count!isSquareFree(iota(1,10000)) == 6083);
	assert(equal(filter!isSquareFree(iota(1,100)), [1,2,3,5,6,7,10,11,13,14,15,17,19,21,22,23,26,29,30,31,33,34,35,37,38,39,41,42,43,46,47,51,53,55,57,58,59,61,62,65,66,67,69,70,71,73,74,77,78,79,82,83,85,86,87,89,91,93,94,95,97][]));
}

/**
 * Find the smallest primitive root modulo n. A root exists whenever n is
 * 1, 2 or 4 or is of the form p^i or 2*p^i for p an odd prime.
 */
long primitiveRoot(long n)
{
	if(n <= 1)
		throw new Exception("invalid modulus");
	if(n == 2)
		return 1;

	auto p = phi(n);
	auto fs = Factorization(p);

	outer: for(long x = 2; x < n; ++x)
	{
		if(gcd(x,n) != 1)
			continue outer;

		foreach(f; fs)
			if(powmod(x, p/f[0], n) == 1)
				continue outer;

		return x;
	}

	throw new Exception("no primitive root found");
}

/** check if x is a primitive root mod n */
bool isPrimitiveRoot(long x, long n)
{
	assert(n > 1);
	assert(0 <= x && x < n);

	if(n == 2)
		return x % n != 0;

	if(gcd(x,n) != 1)
		return false;

	auto p = phi(n);
	foreach(f; factors(p))
		if(powmod(x, p/f[0], n) == 1)
			return false;

	return true;
}

/**
 * Find all d'th roots of 1 modulo m.
 * Currently only implemented if there exists a primitive root modulo n.
 */
Array!long unityRoots(long d, long m)
{
	auto n = phi(m);	// order of the multiplicative group (Z/mZ)*
	auto x = primitiveRoot(m);	// generator of that group (assuming it is cyclic)

	Array!long r;
	for(long k = 0; k < gcd(d,n); ++k)
		r.pushBack(powmod(x, n/gcd(d,n)*k, m));
	return r;
}

/** Caclulate f(f(f(...f(x)...))) using Brent's cycle finding method */
T functionPower(alias f, T)(long n, T x0)
{
	assert(0 <= n);

	T x = x0;
	T y = x;
	long safe = 0;

	for(long i = 0; i < n; ++i, x = unaryFun!f(x))
	{
		if(i != 0 && x == y) // found cycle ?
		{
			// advance by whole cycles
			long len = i - safe;
			assert(len > 0);
			long c = (n-i)/len;
			i += c*len;

			// compute last (incomplete) cycle and return
			for(; i < n; ++i)
				x = unaryFun!f(x);
			return x;
		}

		if(i >= 2*safe) // new safe point
		{
			safe = i;
			y = x;
		}
	}

	return x; // only if no cycle was found
}

unittest
{
	assert(functionPower!"(a+1)%5"(10,0) == 0);
	assert(functionPower!"(a+1)%5"(123456789123456789L,0) == 4);
}
