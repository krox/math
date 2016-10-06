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
private import std.functional : unaryFun, binaryFun;
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

	string toString() const pure nothrow @property
	{
		return "["~to!string(x)~"]";
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

	// excluding 2 and 3 as special cases, all primes have the form 6*k +- 1
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
		cache = assumeUnique(calculatePrimesBelow(limit).release);
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

Array!long divisors(long n)
{
    // factor the number
    auto fs = factor(n);

    // number of divisors (equal to tau(this.multiply))
    long count = 1;
    foreach(f; fs.factors)
        count *= f[1] + 1;

    Array!long d;
    d.reserve(count);
    d.pushBack(1);

    foreach(f; fs.factors)
    {
        long oldCount = d.length;
        for(long i = 0; i < f[1]*oldCount; ++i)
            d.pushBack(f[0] * d[$-oldCount]);
    }

    assert(d.length == count);
    sort(d[]);
    return d;
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
	assert(n < int.max);	// need to use "mulmod" and such if you want this

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
			r = (r * (a + 1 - i)) % p;
			s = (s * i) % p;
		}
	}

	return (r * invmod(s,p)) % p;
}

/**
 * caclulate 1^k + 2^k + .. + n^k using Faulhaber's formula
 *
 * currently only implemented for k <= 3, because the general case needs
 * Bernoulli numbers which are not whole numbers so I need to think about a
 * beautiful way to write that (without rounding or actual rational arithmetic)
 */
long powerSum(long k, long n) pure nothrow
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

	for(int n = 0; n < 61; ++n)
		for(int k = 0; k < 61; ++k)
			foreach(p; [2,3,5,7,11,13,17,19])
				assert(binomial(n,k)%p == binomialMod(n,k,p));

	for(int 1 = 0; n < 61; ++n)
		for(int k = 1; k < 61; ++k)
			assert(binomial(n,k) == binomial(n-1,k) + binomial(n-1,k-1));

	assert(powerSum(0, 5) == 1 + 1 + 1 + 1 + 1);
	assert(powerSum(1, 5) == 1 + 2 + 3 + 4 + 5);
	assert(powerSum(2, 5) == 1*1 + 2*2 + 3*3 + 4*4 + 5*5);
	assert(powerSum(3, 5) == 1*1*1 + 2*2*2 + 3*3*3 + 4*4*4 + 5*5*5);
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

	enum tableLimit = 1_000_000L;
	alias f = binaryFun!(fun, "p", "e");
	alias mul = binaryFun!(mult, "a", "b");

	private Array!long table;

	/** make a table of all values n <= limit */
	void makeTable(long limit)
	{
		if(limit < cast(long)table.length)
			return;

		table.assign(limit+1, neutral);
		foreach(p; primesBelow(limit+1))
			for(long x = p; x < table.length; x += p)
				table[x] = mul(table[x], f(p, powerOf(x, p)));
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
		foreach(i; factor(n))
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
	foreach(f; factor(n))
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
 * determine if a is square-free (i.e. not divisable by a square other than 1)
 */
bool isSquareFree(long a)
{
	foreach(p; primesBelow(sqrti(a)+1))
		if(a % (p*p) == 0)
			return false;

	return true;
}

/** determine if a is a perfect square */
bool isSquare(long a) pure nothrow
{
	if(a < 0)
		return false;
	auto r = sqrti(a);
	return r*r == a;
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
