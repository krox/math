module math.prime;

/**
 * This module is only for 'small' prime numbers, i.e. numbers below 2^60 or so.
 * For larger numbers see math.integer, which builds upon the GMP library.
 */
import jive.array;
import jive.bitarray;
private import std.math : log;
private import core.bitop : bsf;

/** generate all prime numbers below n using a simple sieve */
Array!ulong primesBelow(ulong n)
{
	auto b = BitArray(cast(size_t)n/2);
	Array!ulong primes;
	if(n <= 10)
		primes.reserve(4);
	else
		primes.reserve(cast(size_t)(n/log(n)*(1+1/log(n)+2.51/log(n)/log(n))));

	if(n <= 2)
		return primes;
	primes.pushBack(2);
	if(n <= 3)
		return primes;

	size_t k = 1;
	for( ; /*k < n/2*/; ++k)
		if(!b[k])
		{
			ulong p = 2*k+1;

			if(p*p > n)
				break;

			primes.pushBack(p);

			for(ulong s = p*p/2; s < n/2; s += p)
				b[cast(size_t)s] = true;
		}


	for( ; k < n/2; ++k)
		if(!b[k])
			primes.pushBack(2*k+1);

	return primes;
}

/**
 * calculate b^e % m using binary exponentiation
 * NOTE: only works for m < 2^32
 */
ulong powmod(ulong base, ulong exp, ulong mod) pure
{
	base %= mod;
	ulong r = 1;
	while(exp)
	{
		if(exp&1)
			r = (r*base) % mod;
		exp >>= 1;
		base = (base*base) % mod;
	}
	return r;
}

ulong modinv(ulong x, ulong mod)
{
	return powmod(x,mod-2,mod);	// only works for mod prime (TODO)
}

long gcd(long a, long b)
{
	if(b == 0)
		return a;
	else
		return gcd(b, a%b);
}

/**
 * tests wether n is a strong probable prime to base a
 */
bool isSPRP(ulong a, ulong n)
{
	// the original definition assumes a < n-1, we use a straight-forward generalization
	a %= n;
	if(a == 0 || a == 1 || a == -1)
		return true;

	long s = bsf(cast(size_t)(n-1));
	long d = (n-1) >> s;

	// now it is n-1 = t*2^s, t odd

	a = powmod(a, d, n);
	if(a == 1 || a == n-1)
		return true;

	while(--s)
	{
		a = (a*a)%n;
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
