module math.floating;

private import std.string : toStringz;
//private import std.format;
private import std.conv : to;
private import std.typecons : Rebindable;
private import std.exception : assumeUnique;
private import std.random : unpredictableSeed;
private import std.algorithm : canFind;
//private import jive.array;
private import math.integer;
private import bindings.gmp;
private import bindings.mpfr;

/**
 * Arbitrary precision floating point type using the GNU MPFR library.
 */
struct Floating
{
	/**
	 * All calculations will be done in this precision.
	 */
	static mpfr_prec_t prec = 53; // this emulates IEEE double precision

	//////////////////////////////////////////////////////////////////////
	/// constructors
	//////////////////////////////////////////////////////////////////////

	static enum nan = Floating.init;

	this(immutable MpfrFloat f) pure nothrow
	{
		this.f = f;
	}

	/** constructor for given value */
	this(int v) nothrow
	{
		auto r = new MpfrFloat(prec);
		mpfr_set_si(r, v, MPFR_RNDN);
		this(cast(immutable)r);
	}

	/** ditto */
	this(double v) nothrow
	{
		auto r = new MpfrFloat(prec);
		mpfr_set_d(r, v, MPFR_RNDN);
		this(cast(immutable)r);
	}

	/** ditto */
	this(string v)
	{
		auto r = new MpfrFloat(prec);
		if(mpfr_set_str(r, toStringz(v), 0, MPFR_RNDN) != 0)
			throw new Exception("bad string given to MPFR");
		this(cast(immutable)r);
	}

	/** random number in [0,1) */
	static Floating random()
	{
		auto r = new MpfrFloat(prec);
		mpfr_urandomb(r, &rand);
		return Floating(cast(immutable)r);
	}


	//////////////////////////////////////////////////////////////////////
	/// conversion
	//////////////////////////////////////////////////////////////////////

	/** returns value as (decimal) string */
	string toString() const @property
	{
		char* buf;
		mpfr_asprintf(&buf, "%.30Rg", ptr);
		string r = to!string(buf);
		mpfr_free_str(buf);
		return r;
	}

	/** returns value as int, asserting it is small enough */
	int opCast(T)() const pure
		if(is(T == int))
	{
		assert(int.min <= this && this <= int.max);
		return cast(int)mpfr_get_si(ptr, MPFR_RNDN);
	}

	double opCast(T)() const pure
		if(is(T == double))
	{
		return mpfr_get_d(ptr, MPFR_RNDN);
	}


	//////////////////////////////////////////////////////////////////////
	/// size metric and comparisons
	//////////////////////////////////////////////////////////////////////

	bool isNaN() const pure nothrow @property
	{
		return f is null || mpfr_nan_p(this);
	}

	bool opEquals(double b) const pure nothrow
	{
		return mpfr_cmp_d(this, b) == 0;
	}

	bool opEquals(Integer b) const pure nothrow
	{
		return mpfr_cmp_z(this, b) == 0;
	}

	bool opEquals(Floating b) const pure nothrow
	{
		return mpfr_cmp(this, b) == 0;
	}

	int opCmp(double b) const pure nothrow
	{
		return mpfr_cmp_d(this, b);
	}

	int opCmp(Integer b) const pure nothrow
	{
		return mpfr_cmp_z(this, b);
	}

	int opCmp(Floating b) const pure nothrow
	{
		return mpfr_cmp(this, b);
	}


	//////////////////////////////////////////////////////////////////////
	/// arithmetic operations
	//////////////////////////////////////////////////////////////////////

	Floating opUnary(string op)() const
	{
		auto r = new MpfrFloat(prec);

		     static if(op == "+") mpfr_set(r, this, MPFR_RNDN);
		else static if(op == "-") mpfr_neg(r, this, MPFR_RNRN);
		else static assert(false);

		return Floating(cast(immutable)r);
	}

	Floating opBinary(string op)(Floating b) const
	{
		auto r = new MpfrFloat(prec);

		     static if(op == "+") mpfr_add(r, this, b, MPFR_RNDN);
		else static if(op == "-") mpfr_sub(r, this, b, MPFR_RNDN);
		else static if(op == "*") mpfr_mul(r, this, b, MPFR_RNDN);
		else static if(op == "/") mpfr_div(r, this, b, MPFR_RNDN);
		else static if(op == "^^") mpfr_pow(r, this, b, MPFR_RNDN);
		else static assert(false);

		return Floating(cast(immutable)r);
	}

	Floating opBinary(string op)(double b) const
	{
		auto r = new MpfrFloat(prec);

		     static if(op == "+") mpfr_add_d(r, this, b, MPFR_RNDN);
		else static if(op == "-") mpfr_sub_d(r, this, b, MPFR_RNDN);
		else static if(op == "*") mpfr_mul_d(r, this, b, MPFR_RNDN);
		else static if(op == "/") mpfr_div_d(r, this, b, MPFR_RNDN);
		else static assert(false);

		return Floating(cast(immutable)r);
	}

	Floating opBinaryRight(string op)(double b) const
	{
		auto r = new MpfrFloat(prec);

		     static if(op == "+") mpfr_add_d(r, this, b, MPFR_RNDN);
		else static if(op == "-") mpfr_d_sub(r, b, this, MPFR_RNDN);
		else static if(op == "*") mpfr_mul_d(r, this, b, MPFR_RNDN);
		else static if(op == "/") mpfr_d_div(r, b, this, MPFR_RNDN);
		else static assert(false);

		return Floating(cast(immutable)r);
	}

	/** opOpAssign for convenience */
	Floating opOpAssign(string op, S)(S b)
	{
		this = this.opBinary!op(b);
		return this;
	}

	//////////////////////////////////////////////////////////////////////
	/// mathematical functions
	//////////////////////////////////////////////////////////////////////

	private enum unaryFuns = [
		"abs", "round", "trunc", "ceil", "floor", "frac", // rounding
		"sin", "cos", "tan", "cot", "sec", "csc", // basic trigonometry
		"sinh", "cosh", "tanh", "coth", "sech", "csch", // hyperbolic trigonometry
		"asin", "acos", "atan", "asinh", "acosh", "atanh", // inverse trigonometry
		"sqrt", "cbrt", "sqr", "rec_sqrt", // roots and squares
		"log", "log2", "log10", "log1p", // logarithms
		"exp", "exp2", "exp10", "expm1", // exponentials
		"eint", "li2", "erf", "erfc", "gamma", "lngamma", // special functions
		"digamma", "zeta", "j0", "j1", "y0", "y1", "ai", // ditto
		];

	/** unary functions */
	Floating opDispatch(string op)() const
		if(canFind(unaryFuns, op))
	{
		auto r = new MpfrFloat(prec);
		static if(canFind(["round", "trunc", "ceil", "floor"], op))
			mixin("mpfr_"~op~"(r, this);");
		else
			mixin("mpfr_"~op~"(r, this, MPFR_RNDN);");
		return Floating(cast(immutable)r);
	}

	private enum binaryFuns = [
		"pow", "agm", "reldiff", "fmod", "atan2", "hypot", "min", "max", "dim",
		];

	/** binary functions */
	Floating opDispatch(string op)(Floating b) const
		if(canFind(binaryFuns, op))
	{
		auto r = new MpfrFloat(prec);
		mixin("mpfr_"~op~"(r, this, b, MPFR_RNDN);");
		return Floating(cast(immutable)r);
	}

	alias sqAbs = opDispatch!"sqr";

	//////////////////////////////////////////////////////////////////////
	/// internals
	//////////////////////////////////////////////////////////////////////

	Rebindable!(immutable(MpfrFloat)) f;
	alias f this;

	private static mp_randstate_t rand;

	static this()
	{
		__gmp_randinit_default(&rand);
		__gmp_randseed_ui(&rand, unpredictableSeed);
	}
}

/**
 * convenience wrapper for a (mutable) MPFR float.
 */
final class MpfrFloat
{
	mpfr_t f;

	inout(mpfr_t)* ptr() inout pure nothrow @property
	{
		return &f;
	}

	alias ptr this;

	this(mpfr_prec_t prec) pure nothrow
	{
		mpfr_init2(this, prec);
	}

	~this() pure nothrow
	{
		mpfr_clear(this);
	}
}
