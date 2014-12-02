module math.numtheory;

/**
 * Number Theory functions for 64-bit integers and smaller
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


/** generate all prime numbers below n using a simple sieve */
Array!long primesBelow(long n)
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
 * Calculate base^exp % mod using binary exponentiation.
 * Sign of mod is ignored and result is always >= 0.
 */
long powmod(long base, long exp, long mod) //pure
{
	if(mod < 0)
		mod = -mod;

	if(mod > uint.max)
		throw new Exception("exponentiation with large modulus not supported (TODO)");

	base %= mod;
	if(base < 0)
		base += mod;

	if(exp < 0)
	{
		base = modinv(base, mod);
		exp = -exp;
	}

	// signed -> unsigned cast is safe now, and this way, slightly
	// larger values for mod can be allowed: (m-1)^2 has to fit into ulong
	ulong b = base;
	ulong e = exp;
	ulong m = mod;
	ulong r = 1;

	while(e)
	{
		if(e&1)
			r = (r*b)%m;
		e >>= 1;
		b = (b*b) % m;
	}

	return r; // cast back to signed should be fine now
}

/**
 * Compute inverse of x modulo mod.
 * Sign of mod is ignored, result is always > 0.
 * Throws if gcd(x,mod) != 1.
 */
long modinv(long x, long mod) //pure
{
	// cool trick if mod is known to be a prime
	//return powmod(x,mod-2,mod);

	if(mod < 0)
		mod = -mod;
	x %= mod;

	// invariant: a = y * x + _ * mod
	long a_old = mod;
	long a = x;
	long y_old = 0;
	long y = 1;

	while(a != 0)
	{
		long q = a_old/a;

		long t = a;
		a = a_old%a;
		a_old = t;

		t = y;
		y = y_old - q*y;
		y_old = t;
	}

	if(a_old != 1)
		throw new Exception("number not invertible (gcd(x,mod) != 1)");

	return y_old>=0?y_old:y_old+mod;
}

/**
 * Compute greatest common divisor of a and b.
 * Sign of a and b is ignored, result is always >= 0.
 * Convention: gcd(0,x) = x = gcd(x,0).
 */
long gcd(long a, long b)
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

/** tests wether n is a strong probable prime to base a */
bool isSPRP(long a, long n)
{
	if(n > uint.max) // binary exponentiation wont work (also bsf on 32 bit system)
		throw new Exception("overflow in isSPRP()");
	if(a <= 0)
		throw new Exception("a-SPRP is only defined for a>0");

	if(n < 2)
		return false;

	// the original definition assumes a < n-1, we use a straight-forward generalization
	a %= n;
	if(a == 0 || a == 1 || a == n-1)
		return true;

	long s = bsf(cast(size_t)(n-1));
	long d = (n-1) >> s;

	// now it is n-1 = t*2^s, t odd

	a = powmod(a, d, n);
	if(a == 1 || a == n-1)
		return true;

	while(--s)
	{
		a = (cast(ulong)a*cast(ulong)a)%cast(ulong)n;
		if(a == n-1)
			return true;
	}
	return false;
}

/**
 * test wether a number n < 3_071_837_692_357_849 is prime
 * see http://priv.ckp.pl/wizykowski/sprp.pdf for details
 */
bool isPrime(long n) //pure
{
	if(n<2)
		return false;
	if(n%2==0)
		return n == 2;

	// TODO: better optimization for small n
	// i.e.: table-lookup and/or more/smarter trial-division
	for(long d = 3; d < 20; d += 2)
	{
		if(d*d > n)
			return true;
		if(n%d == 0)
			return false;
	}

	// NOTE: only use bases < 2^32. isSPRP() and powmod cant handle higher values
	if(n < 227_132_641)
		return isSPRP(          660, n)
		    && isSPRP(   56_928_287, n);

	if(n < 105_936_894_253)
		return isSPRP(2, n)
		    && isSPRP(1_005_905_886, n)
		    && isSPRP(1_340_600_841, n);

	if(n < 31_858_317_218_647)
		return isSPRP(            2, n)
		    && isSPRP(      642_735, n)
		    && isSPRP(  553_174_392, n)
		    && isSPRP(3_046_413_974, n);

	if(n < 3_071_837_692_357_849)
		return isSPRP(            2, n)
		    && isSPRP(       75_088, n)
		    && isSPRP(      642_735, n)
		    && isSPRP(  203_659_041, n)
		    && isSPRP(3_613_982_119, n);


	throw new Exception("number to high for prime testing");
}
