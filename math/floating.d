module math.floating;

private import std.string : toStringz;
private import std.format;
private import std.conv : to;
private import std.exception : assumeUnique;
private import std.random : unpredictableSeed;
private import std.typecons : Rebindable;
private import math.mpfr;
private import math.gmp;

/**
 * Arbitrary precision floating point type using the GNU MPFR library.
 *
 * Copy-on-write semantics, i.e. every result is newly allocated. This means
 * it should be very compatible with float/double/real but kinda inefficient.
 * May change in future when I figure out a good way to do expression templates
 * or something...
 */
struct Floating
{
	static mpfr_prec_t precision;

	//////////////////////////////////////////////////////////////////////
	/// constructors
	//////////////////////////////////////////////////////////////////////

	this(immutable MpfrFloat z)
	{
		this.z = z;
	}

	/** constructor for given value */
	this(int v)
	{
		auto r = new MpfrFloat(precision);
		mpfr_set_si(r.ptr, v, MPFR_RNDN);
		this(cast(immutable)r);
	}

	/** ditto */
	this(double v)
	{
		auto r = new MpfrFloat(precision);
		mpfr_set_d(r.ptr, v, MPFR_RNDN);
		this(cast(immutable)r);
	}

	/** ditto */
	this(string v)
	{
		// TODO: throw exception on bad strings
		auto r = new MpfrFloat(precision);
		mpfr_set_str(r.ptr, toStringz(v), 0, MPFR_RNDN);
		this(cast(immutable)r);
	}

	/** random in [0,1) */
	static Floating random()
	{
		auto r = new MpfrFloat(precision);
		mpfr_urandomb(r.ptr, &rand);
		return Floating(cast(immutable)r);
	}


	//////////////////////////////////////////////////////////////////////
	/// conversion
	//////////////////////////////////////////////////////////////////////

	/** returns value as (decimal) string */
	string toString() const @property
	{
		char* buf;
		mpfr_asprintf(&buf, "%Rg", ptr);
		string r = to!string(buf);
		mpfr_free_str(buf);
		return r;
	}

	/** returns (rounded) value as double */
	double opCast(T)() const
		if(is(T == double))
	{
		return mpfr_get_d(ptr, MPFR_RNDN);
	}


	//////////////////////////////////////////////////////////////////////
	/// size metric and comparisons
	//////////////////////////////////////////////////////////////////////

	bool isNan() const @property
	{
		return mpfr_nan_p(ptr) != 0;
	}

	bool opEquals(T)(T b) const
	{
		     static if(is(T == int))      return mpfr_cmp_si(ptr, b) == 0;
		else static if(is(T == double))   return mpfr_cmp_d(ptr, b) == 0;
		else static if(is(T == Floating)) return mpfr_equal_p(ptr, b.ptr) != 0;
		else static assert(false);
	}

	int opCmp(T)(T b) const
	{
		     static if(is(T == int))      return mpfr_cmp_si(ptr, b);
		else static if(is(T == double))   return mpfr_cmp_d(ptr, b);
		else static if(is(T == Floating)) return mpfr_cmp(ptr, b.ptr);
		else static assert(false);
	}


	//////////////////////////////////////////////////////////////////////
	/// arithmetic operations
	//////////////////////////////////////////////////////////////////////

	Integer opUnary(string op)() const
		if(op == "-")
	{
		auto r = new MpfrFloat(precision);
		mpfr_neg(r.ptr, ptr, MPFR_RNDN);
		return Floating(cast(immutable)r);
	}

	Floating opBinary(string op)(Floating b) const
	{
		auto r = new MpfrFloat(precision);
		mpfrBinary!op(r.ptr, ptr, b.ptr);
		return Floating(cast(immutable)r);
	}

	Floating opBinary(string op)(int b) const
	{
		auto r = new MpfrFloat(precision);
		mpfrBinary!op(r.ptr, ptr, b);
		return Floating(cast(immutable)r);
	}

	Floating opBinaryRight(string op)(int a) const
	{
		auto r = new MpfrFloat(precision);
		mpfrBinary!op(r.ptr, a, ptr);
		return Floating(cast(immutable)r);
	}


	//////////////////////////////////////////////////////////////////////
	/// internals
	//////////////////////////////////////////////////////////////////////

	Rebindable!(immutable(MpfrFloat)) z;

	immutable(mpfr_t)* ptr() const @property
	{
		return &z.z;
	}

	private static mp_randstate_t rand;

	static this()
	{
		__gmp_randinit_default(&rand);
		__gmp_randseed_ui(&rand, unpredictableSeed);
	}
}
