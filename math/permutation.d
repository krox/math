module math.permutation;

private import std.algorithm : min, max;
private import std.conv : to;
private import std.exception : assumeUnique;
private import std.random : uniform;

private import jive.bitarray;

private import math.integer;

struct Permutation
{
	immutable(int)[] coeffs;

	this(immutable(int)[] coeffs)
	{
		while(coeffs.length && coeffs[$-1] == coeffs.length-1)
			coeffs = coeffs[0..$-1];
		this.coeffs = coeffs;
	}

	static Permutation random(int d) @property
	{
		assert(d >= 0);
		auto r = new int[d];
		for(int i = 0; i < d; ++i)
		{
			int j = uniform(0,i+1);
			r[i] = r[j];
			r[j] = i;
		}
		return Permutation(assumeUnique(r));
	}

	int degree() const @property
	{
		return cast(int)coeffs.length;
	}

	string toString() const @property
	{
		if(degree == 0)
			return "()";

		auto done = BitArray(degree);
		string s;
		for(int i = 0; i < degree; ++i)
			if(!done[i] && coeffs[i] != i)
			{
				done[i] = true;
				s ~= "(" ~ to!string(i);
				for(int j = coeffs[i]; j != i; j = coeffs[j])
				{
					done[j] = true;
					s ~= ","~to!string(j);
				}
				s ~= ")";
			}
		return s;
	}

	Integer order() const @property
	{
		Integer r = 1;

		auto done = BitArray(degree);
		for(int i = 0; i < degree; ++i)
			if(!done[i] && coeffs[i] != i)
			{
				done[i] = true;
				int l = 1;
				for(int j = coeffs[i]; j != i; j = coeffs[j])
				{
					done[j] = true;
					++l;
				}
				r = lcm(r, l);
			}

		return r;
	}

	int opCall(int i) const
	{
		if(i < 0 || i >= degree)
			return i;
		return coeffs[i];
	}

	Permutation opBinary(string op)(Permutation b) const
		if(op == "*")
	{
		int d = max(degree, b.degree);
		auto r = new int[d];

		for(int i = 0; i < d; ++i)
			r[i] = opCall(b(i));

		return Permutation(assumeUnique(r));
	}

	Permutation inverse() const @property
	{
		auto r = new int[degree];

		for(int i = 0; i < degree; ++i)
			r[opCall(i)] = i;
		return Permutation(assumeUnique(r));
	}
}
