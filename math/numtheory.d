module math.numtheory;

/**
 * Number Theory and combinatorical functions for 64-bit integers and smaller
 *
 * For similar functionality using arbitrary large numbers, see math.integer.
 * Also, there is no special consideration for architectures without a native
 * 64-bit integer type, so performance might be bad on such CPUs, even in
 * cases where 32-bit computions would suffice.
 */

import jive.array;
import jive.bitarray;
private import std.math : log, sqrt;
private import core.bitop : bsf;
private import std.typecons;
private import std.algorithm;
private import std.range;
private import std.exception;
private import std.format;
private import std.functional : unaryFun;


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
	for(; b != 0; b >>= 1, a = mulmod(a, a, m))
		if(b & 1)
			r = mulmod(r, a, m);
	return r;
}

/**
 * calculate the modular inverse a^-1 % m
 * conditions: m > 0 and 0 <= a,result < m
 * triggers division by zero if gcd(x, m) != 1
 */
long invmod(long a, long m) pure nothrow
{
	assert(m > 0);
	assert(0 <= a && a < m);

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
	return b1;
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

	IntMod opBinary(string op)(int b) const pure nothrow
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

	IntMod opBinary(string op)(int e) const pure nothrow
		if(op == "^^")
	{
		return IntMod(powmod(x, e, n), n);
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

	string toString() const pure nothrow @property
	{
		return "["~to!string(x)~"]";
	}
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

/** generate all prime numbers below n using a simple sieve */
Array!long calculatePrimesBelow(long n)
{
	if(n < 0)
		n = 0;

	if(n > size_t.max) // can only trigger on a 32 bit system
		throw new Exception("can not calculate that many primes");

	// (exclusive) limit for relevant prime divisors
	long limit = cast(long)sqrt(cast(double)n) + 1;
	assert(limit*limit > n); // not sure about floating-point rounding for very large n ...

	// excluding 2 and 3 as special cases, all primes have the form 5*k +- 1
	auto b5 = BitArray(cast(size_t)n/6); // b5[k] represents 6*k+5
	auto b7 = BitArray(cast(size_t)n/6); // b7[k] represents 6*k+7

	// mark all non-primes
	// NOTE: in older versions of DMD, making k an ulong triggers a strange bug:
	// issues.dlang.org/show_bug.cgi?id=13023
	for(long k = 0; k < n/6; ++k)
	{
		if(!b5[cast(size_t)k])
		{
			long p = 6*k+5;

			if(p >= limit)
				break;

			for(long s = (p*(p+2)-5)/6; s < b5.length; s += p)
				b5[cast(size_t)s] = true;

			for(long s = (p*p-7)/6; s < b7.length; s += p)
				b7[cast(size_t)s] = true;
		}

		if(!b7[cast(size_t)k])
		{
			long p = 6*k+7;

			if(p >= limit)
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
	for(long k = 0; k < n/6; ++k)
	{
		if(!b5[cast(size_t)k])
			primes.pushBack(6L*k+5);
		if(!b7[cast(size_t)k])
			primes.pushBack(6L*k+7);
	}

	// now we might have computed slightly more primes than requested,
	// so we remove them again (simpler than doing it right in the beginning)
	while(!primes.empty && primes[$-1] >= n)
		primes.popBack;

	return primes;
}

/**
 * returns all primes in [a..b) using cached results from calculatePrimesBelow
 */
immutable(long)[] primesBetween(long a, long b)
{
	static immutable(long)[] cache;
	static long limit = 2;

	if(b > limit)
	{
		limit = max(b, limit+limit/2);
		cache = assumeUnique(cast(long[])calculatePrimesBelow(limit));
	}

	return cache[].assumeSorted.upperBound(a-1).lowerBound(b).release;
}

/**
 * alias for primesBetween(0, n)
 */
immutable(long)[] primesBelow(long n)
{
	return primesBetween(0, n);
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

/**
 * test if n is prime
 *
 * see http://priv.ckp.pl/wizykowski/sprp.pdf for implementation details and
 * http://miller-rabin.appspot.com for the actual values used therein.
 */
bool isPrime(long n) pure nothrow
{
	if(n < 2)
		return false;

	if(n % 2 == 0)
		return n == 2;

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

unittest
{
	assert(primesBetween(3,20) == [3,5,7,11,13,17,19]);
	assert(isPrime(1000000007));
	assert(isPrime(9223372036854775783L)); // largest 63 bit prime
	assert(!isPrime(1000000007L*1000000009L));
}


//////////////////////////////////////////////////////////////////////
/// integer factorization
//////////////////////////////////////////////////////////////////////

/**
 * convenience wrapper around Array!(Tuple!(long,int)) for integer factorization
 */
struct Factorization
{
	Array!(Tuple!(long,int)) factors;
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
	long multiply() pure nothrow const @property
	{
		long r = 1;
		foreach(f; factors)
			for(int i = 0; i < f[1]; ++i)
				r *= f[0];
		return r;
	}

	/**
	 * checks if all factors are prime
	 */
	bool allPrime() pure nothrow const @property
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
Factorization factor(long n)
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
			long d = f[i][0];
			for(int c = 1; d == f[i][0]; ++c)
				d = findFactor(f[i][0], 0, c);

			f[i][0] /= d;
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
long findFactor(long n, long x0, long c) pure nothrow
{
	assert(n > 0);
	assert(0 < c && c < n);

	long x = x0; // arbitrary start value. Actual randomization might be good...
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
	assert(factor(2*2*2*2*3*5*5*7)[] == [tuple(2,4),tuple(3,1),tuple(5,2),tuple(7,1)]);
	assert(factor(1000000007L*1000000009L)[] == [tuple(1000000007L,1), tuple(1000000009L,1)]);
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
long binomial(long n, long k)
{
	assert(0 <= n);
	assert(0 <= k && k <= n);

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
 * caclulate 1^k + 2^k + .. + n^k using Faulhaber's formula
 *
 * currently only implemented for k <= 3, because the general case needs
 * Bernoulli numbers which are not whole numbers so I need to think about a
 * beautiful way to write that (without rounding or actual rational arithmetic)
 */
long powerSum(long k, long n)
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

unittest
{
	assert(factorial(0) == 1);
	assert(factorial(5) == 120);
	assert(binomial(0,0) == 1);
	assert(binomial(4,2) == 6);

	assert(powerSum(0, 5) == 1 + 1 + 1 + 1 + 1);
	assert(powerSum(1, 5) == 1 + 2 + 3 + 4 + 5);
	assert(powerSum(2, 5) == 1*1 + 2*2 + 3*3 + 4*4 + 5*5);
	assert(powerSum(3, 5) == 1*1*1 + 2*2*2 + 3*3*3 + 4*4*4 + 5*5*5);
}


//////////////////////////////////////////////////////////////////////
/// number theoretic functions
//////////////////////////////////////////////////////////////////////


// TODO: templated function for general multiplicative functions

/**
 * Euler's totient function.
 */
long phi(long n)
{
	foreach(p; factor(n))
		n = n/p[0]*(p[0]-1);
	return n;
}

//////////////////////////////////////////////////////////////////////
/// other stuff not fitting in the above categories
//////////////////////////////////////////////////////////////////////

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
 * Find the smallest primitive root modulo n. A root exists whenever n is
 * 1, 2 or 4 or is of the form p^i or 2*p^i for p an odd prime.
 */
long primitiveRoot(long n)
{
	import std.stdio;

	if(n <= 1)
		throw new Exception("invalid modulus");
	if(n == 2)
		return 1;

	auto p = phi(n);
	auto fs = factor(p);

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
