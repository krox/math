module math.factor_integer;

private import std.algorithm;
private import std.typetuple;
private import std.typecons;

private import jive.array;

private import math.integer;
private import math.prime;


Array!(Tuple!(Integer,int)) factor(Integer n)
{
	if(n <= 0)
		throw new Exception("can only factor postive integers");

	Array!(Tuple!(Integer,int)) r;

	if(n.isPrime)
	{
		r.pushBack(tuple(n,1));
		return r;
	}

	auto z = new GmpInteger;
	__gmpz_set(z.ptr, n.ptr);

	foreach(fac; removeLowFactors(z))
		r.pushBack(tuple(Integer(fac[0]),fac[1]));

	Array!(Tuple!(GmpInteger, int)) stack;
	if(__gmpz_cmp_si(z.ptr, 1) != 0)
		stack.pushBack(tuple(z, 1));

	while(!stack.empty)
	{
		auto k = stack.popBack;
		if(__gmpz_probab_prime_p(k[0].ptr, 25))
			r.pushBack(tuple(Integer(cast(immutable)k[0]), k[1]));
		else
		{
			auto t = findFactor(k[0]);
			auto cnt = cast(int)__gmpz_remove(k[0].ptr, k[0].ptr, t.ptr);
			assert(cnt >= 1);
			if(__gmpz_cmp_si(k[0].ptr, 1) != 0) // remaining is 1 if k was a power of t
				stack.pushBack(k);
			stack.pushBack(tuple(t, k[1]*cnt));
		}
	}

	sort(r[]);
	return r;
}

private static
{
	enum primorialLimit = 5000;
	GmpInteger primorial;
	Array!ulong primesLow;

	static this()
	{
		primorial = new GmpInteger;
		__gmpz_primorial_ui(primorial.ptr, primorialLimit);
		primesLow = primesBelow(primorialLimit);
	}
}

private Array!ulong factorLowSquareFree(GmpInteger n)
{
	Array!ulong r;

	foreach(p; primesLow)
	{
		if(__gmpz_cmp_si(n.ptr, 1) == 0)
			break;

		if(__gmpz_divisible_ui_p(n.ptr, p))
		{
			__gmpz_divexact_ui(n.ptr, n.ptr, p);
			r.pushBack(p);
		}
	}

	assert(__gmpz_cmp_si(n.ptr, 1) == 0);
	return r;
}

private Array!(Tuple!(ulong,int)) removeLowFactors(GmpInteger n)
{
	Array!(Tuple!(ulong,int)) r;

	scope g = new GmpInteger;
	scope w = new GmpInteger;
	scope a = new GmpInteger;

	__gmpz_gcd(g.ptr, n.ptr, primorial.ptr);

	for(int i = 1; __gmpz_cmp_si(g.ptr, 1) != 0; ++i)
	{
		__gmpz_divexact(n.ptr, n.ptr, g.ptr);
		__gmpz_gcd(w.ptr, g.ptr, n.ptr);
		__gmpz_divexact(a.ptr, g.ptr, w.ptr);

		foreach(p; factorLowSquareFree(a))
			r.pushBack(tuple(p, i));

		__gmpz_swap(w.ptr, g.ptr);
	}

	return r;
}

private GmpInteger findFactor(const GmpInteger n)
{
	auto inc = new GmpInteger;
	__gmpz_sqrt(inc.ptr, n.ptr);
	if(__gmpz_perfect_square_p(n.ptr))
		return inc;

	scope a = new GmpInteger;
	__gmpz_mul(a.ptr, inc.ptr, inc.ptr);
	__gmpz_sub(a.ptr, a.ptr, n.ptr);

	__gmpz_mul_ui(inc.ptr, inc.ptr, 2);
	__gmpz_add_ui(inc.ptr, inc.ptr, 1);

	while(true)
	{
		__gmpz_add(a.ptr, a.ptr, inc.ptr);
		__gmpz_add_ui(inc.ptr, inc.ptr, 2);

		if(__gmpz_perfect_square_p(a.ptr))
		{
			__gmpz_fdiv_q_ui(inc.ptr, inc.ptr, 2); // will round down
			__gmpz_sqrt(a.ptr, a.ptr);
			__gmpz_sub(inc.ptr, inc.ptr, a.ptr);
			return inc;
		}
	}
}
