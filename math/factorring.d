module math.factorring;

import std.conv : format;
import std.algorithm : move;
import math.integer;


struct Coset(Ring)
{
	Ring val;
	Ring mod;

	this(Ring val, Ring mod)
	{
		this.val = val % mod;
		this.mod = mod;
	}

	this(int val, Ring mod)
	{
		this(Ring(val), mod);
	}

	this(int val, int mod)
	{
		this(Ring(val), Ring(mod));
	}

	static Coset random(Ring mod)
	{
		static if(is(Ring == Integer))
		{
			return Coset(Ring.random(mod), mod); // TODO: use fast no-check constructor
		}
		else assert(false, "TODO");
	}

	string toString() const @property
	{
		return "["~val.toString~"]";
	}

	/** return 1/this */
	Coset inverse() const @property
	{
		static if(is(Ring == Integer))
		{
			return Coset(val.inverseMod(mod), mod);
		}
		else
		{
			assert(false, "TODO");
		}
	}

	Coset opUnary(string op)() const
		if(op == "-")
	{
		return Coset(-val, mod);
	}

	Coset opBinary(string op, T)(T rhs) const
		if(is(T == int) || is(T == Integer) || is(T == Ring))
	{
		     static if(op == "+") return Coset(val + rhs, mod);
		else static if(op == "-") return Coset(val - rhs, mod);
		else static if(op == "*") return Coset(val * rhs, mod);
		else static assert(false, "binary assign '"~op~"' is not defined");
	}

	Coset opBinary(string op)(Coset rhs) const
	{
		static if(op == "/")
			return this * rhs.inverse;
		else
		{
			assert(this.mod == rhs.mod);
			return opBinary!op(rhs.val);
		}
	}

	bool opEquals(Coset r) const
	{
		assert(mod == r.mod);
		return val == r.val;
	}

	bool opEquals(T)(T r) const
		if(!is(T : Coset))
	{
		return opEquals(Coset(r, mod));
	}
}

/** faster special case for when the modulus is a small compile-time integer */
struct IntCoset(ulong mod)
{
	static assert(mod <= uint.max);

	ulong val; // actually always in [0..mod), and mod has to fit in 32 bit

	this(long val)
	{
		val %= mod;
		val += mod;
		val %= mod;
		this.val = val;
	}

	string toString() const @property
	{
		return format("[%s]", val);
	}

	/** return 1/this */
	IntCoset inverse() const @property
	{
		assert(false, "TODO");
	}

	IntCoset opUnary(string op)() const
		if(op == "-")
	{
		return IntCoset(-val);
	}

	IntCoset opBinary(string op)(IntCoset rhs) const
	{
		     static if(op == "+") return IntCoset(val + rhs.val);
		else static if(op == "-") return IntCoset(val - rhs.val);
		else static if(op == "*") return IntCoset(val * rhs.val);
		else static if(op == "/") return IntCoset(val * rhs.inverse.val);
		else static assert(false, "binary assign '"~op~"' is not defined");
	}

	bool opEquals(IntCoset r) const
	{
		return val == r.val;
	}
}
