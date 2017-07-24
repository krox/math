module math.factor;

/**
 * Factorizing algorithms. Currently only for polynomials over finite fields.
 * Might add integer polynomials, or integers itself later.
 */

private import std.typecons;
private import std.typetuple;
private import std.algorithm;
private import std.exception;

private import jive.array;

import math.integer;
import math.finitefield;
import math.polynomial;

/**
 * Split a polynomial over a finite field into irreducible factors.
 * Implemented using the randomized Cantor-Zassenhaus algorithm.
 */
Array!(Tuple!(Polynomial!FFE,int)) factor(Polynomial!FFE f)
{
	if(f.degree < 0 || f.leading != 1)
		throw new Exception("can not factorize non-normalized polynomial");

	auto field = f.coeffs[0].field;
	assert(field !is null);

	alias Factor = Tuple!(Polynomial!FFE,int);
	Array!Factor sff; // square-free factors

	if(f.degree == 0) // factorization of 1 is the empty product
		return sff;

	auto c = gcd(f, f.derivative).normalize; // NOTE: the derivative can be 0 (in which case c = f)
	auto w = f/c;

	for(int i = 1; w.degree != 0; ++i)
	{
		auto y = gcd(w, c).normalize;
		if(y.degree != w.degree)
			sff.pushBack(tuple(w/y, i));
		w = y;
		c = c/y;
	}

	if(c.degree != 0)
		foreach(fac; factor(pthRoot(c)))
			sff.pushBack(tuple(fac[0], fac[1]*field.p));

	Array!Factor r;
	foreach(fac; sff)
		foreach(p; factorSquareFree(fac[0]))
			r.pushBack(tuple(p, fac[1]));

	sort(r[]);
	return r;
}

/**
 * Split a square-free polynomial over a finite field into irreducible factors.
 */
Array!(Polynomial!FFE) factorSquareFree(Polynomial!FFE f)
{
	if(f.degree < 0 || f.leading != 1)
		throw new Exception("can not factorize non-normalized polynomial");

	auto field = f.coeffs[0].field;

	if(!f.isSquareFree)
		throw new Exception("polynomial is not actually square-free");

	Array!(Polynomial!FFE) r;

	for(int i = 1; f.degree > 0; ++i)
	{
		auto x = Polynomial!FFE([FFE(0,field), FFE(1,field)]);
		auto g = x.powmod(Integer(field.q)^^i, f) - x;

		auto fac = gcd(f, g).normalize; // g = gcd(x^q^i - x, f) = all factors of degree i in f
		if(fac.degree <= 0)
			continue;

		r.pushBack(factorEqualDegree(fac, i)[]);
		f = f / fac;
	}

	return r;
}

/**
 * Split a square-free polynomial over a finite field with only factors of degree d into irreducible factors.
 * NOTE: if f is not if this form, this function will get stuck.
 */
Array!(Polynomial!FFE) factorEqualDegree(Polynomial!FFE f, int d)
{
	if(f.degree < 0 || f.leading != 1)
		throw new Exception("can not factorize non-normalized polynomial");

	assert(f.degree % d == 0);

	auto field = f.coeffs[0].field;
	if(field.p % 2 == 0)
		throw new Exception("factorization over fields of even characteristic not implemented (yet)");
	auto e = Integer(field.q)^^d/2; // NOTE: this rounds down

	Array!(Polynomial!FFE) r;
	if(f.degree == 0)
		return r;

	r.reserve(f.degree/d); // without this, appending while iterating would be illegal (and a very hard to find bug)

	r.pushBack(f);

	while(r.length < f.degree/d)
	{
		auto g = Polynomial!FFE.random(f.degree-1, field);
		g = g.powmod(e, f)-Polynomial!FFE(1);
		if(g.degree == 0) // early out for useless g
			continue;

		foreach(ref p; r[])
			if(p.degree > d)
			{
				auto h = gcd(p, g);
				if(h.degree == p.degree || h.degree == 0)
					continue;

				h = h.normalize();
				p = p/h;
				r.pushBack(h);
			}
	}

	return r;
}

/**
 * p-th root of a polynomial over F_p^n
 */
Polynomial!FFE pthRoot(Polynomial!FFE f)
{
	if(f.degree < 0)
		return f;

	auto field = f.coeffs[0].field;
	for(int i = 0; i <= f.degree; ++i)
		if(i%field.p != 0 && f.coeffs[i] != 0)
			throw new Exception("polynomial "~f.toString~" does not have a pth root");

	auto r = new FFE[f.degree / field.p + 1];
	for(int i = 0; i < r.length; ++i)
		r[i] = f.coeffs[i*field.p] ^^ (field.q/field.p);
	return Polynomial!FFE(assumeUnique(r));
}
