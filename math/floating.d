module math.floating;

private import std.string : toStringz;
private import std.format;
private import std.conv : to;
private import std.exception : assumeUnique;
private import std.random : unpredictableSeed;
private import std.typecons : Rebindable;
private import std.traits;
private import math.mpfr;
private import math.gmp;

/**
 * Arbitrary precision floating point type using the GNU MPFR library.
 */
struct Floating
{
	static mpfr_prec_t precision = 64;

	//////////////////////////////////////////////////////////////////////
	/// constructors / destructors
	//////////////////////////////////////////////////////////////////////

	this(this)
	{
		if(z.mpfr_d)
		{
			mpfr_t old = z;
			this.z = this.z.init;
			mpfr_init2(ptr, mpfr_get_prec(&old));
			mpfr_set(ptr, &old, MPFR_RNDN);
		}
	}

	~this()
	{
		if(z.mpfr_d)
			mpfr_clear(ptr);
	}

	private this(mpfr_t z)
	{
		this.z = z;
	}

	/** constructor for given value */
	this(A)(auto ref A a, mpfr_prec_t prec = precision)
	{
		mpfr_init2(ptr, prec);
		opAssign(a);
	}

	auto opAssign(T)(auto ref T expr)
	{
		if(!z.mpfr_d)
			mpfr_init2(ptr, precision);
		expr.toExpr.eval(ptr);
		return this.toExpr();
	}

	/** random in [0,1) */
	static Floating random(mpfr_prec_t prec = precision)
	{
		mpfr_t z;
		mpfr_init2(&z, prec);
		mpfr_urandomb(&z, &rand);
		return Floating(z);
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
	/// misc
	//////////////////////////////////////////////////////////////////////

	bool isNan() const @property
	{
		return mpfr_nan_p(ptr) != 0;
	}


	//////////////////////////////////////////////////////////////////////
	/// arithmetic operations
	//////////////////////////////////////////////////////////////////////

	auto toExpr() const @property
	{
		return ValueExpr!(const(mpfr_t)*)(ptr);
	}

	alias toExpr this;

	//////////////////////////////////////////////////////////////////////
	/// internals
	//////////////////////////////////////////////////////////////////////

	mpfr_t z;

	inout(mpfr_t)* ptr() inout @property
	{
		return &z;
	}

	private static mp_randstate_t rand;

	static this()
	{
		__gmp_randinit_default(&rand);
		__gmp_randseed_ui(&rand, unpredictableSeed);
	}
}

private string unaryDecl(string name)
{
	return "auto "~name~"(A)(auto ref A a)"
	"{ return UnaryExpr!(\""~name~"\", typeof(a.toExpr()))(a.toExpr());}";
}

// rounding
mixin(unaryDecl("abs"));
mixin(unaryDecl("round"));
mixin(unaryDecl("trunc"));
mixin(unaryDecl("ceil"));
mixin(unaryDecl("floor"));
mixin(unaryDecl("frac"));

// trigonometry
mixin(unaryDecl("sin"));
mixin(unaryDecl("cos"));
mixin(unaryDecl("tan"));
mixin(unaryDecl("cot"));
mixin(unaryDecl("sec"));
mixin(unaryDecl("csc"));

// hyperbolic trigonometry
mixin(unaryDecl("sinh"));
mixin(unaryDecl("cosh"));
mixin(unaryDecl("tanh"));
mixin(unaryDecl("coth"));
mixin(unaryDecl("sech"));
mixin(unaryDecl("csch"));

// inverse trigonometry
mixin(unaryDecl("asin"));
mixin(unaryDecl("acos"));
mixin(unaryDecl("atan"));
mixin(unaryDecl("asinh"));
mixin(unaryDecl("acosh"));
mixin(unaryDecl("atanh"));

// roots and squares
mixin(unaryDecl("sqrt"));
mixin(unaryDecl("cbrt"));
mixin(unaryDecl("sqr"));
mixin(unaryDecl("rec_sqrt"));

// exponential and logarithm
mixin(unaryDecl("log"));
mixin(unaryDecl("log2"));
mixin(unaryDecl("log10"));
mixin(unaryDecl("log1p"));
mixin(unaryDecl("exp"));
mixin(unaryDecl("exp2"));
mixin(unaryDecl("exp10"));
mixin(unaryDecl("expm1"));

// other special functions
mixin(unaryDecl("eint"));
mixin(unaryDecl("li2"));
mixin(unaryDecl("erf"));
mixin(unaryDecl("erfc"));
mixin(unaryDecl("gamma"));
mixin(unaryDecl("lngamma"));
mixin(unaryDecl("lgamma"));
mixin(unaryDecl("digamma"));
mixin(unaryDecl("zeta"));
mixin(unaryDecl("j0"));
mixin(unaryDecl("j1"));
mixin(unaryDecl("y0"));
mixin(unaryDecl("y1"));
mixin(unaryDecl("ai"));


private string binaryDecl(string name)
{
	return "auto "~name~"(A, B)(auto ref A a, auto ref B b)"
	"{ return BinaryExpr!(\""~name~"\", typeof(a.toExpr()), typeof(b.toExpr()))"
	"(a.toExpr(), b.toExpr());}";
}

mixin(binaryDecl("pow"));
mixin(binaryDecl("agm"));
mixin(binaryDecl("reldiff"));
mixin(binaryDecl("fmod"));
mixin(binaryDecl("atan2"));
mixin(binaryDecl("hypot")); // sqrt(a*a + b*b)
mixin(binaryDecl("min"));
mixin(binaryDecl("max"));
mixin(binaryDecl("dim")); // abs(a - b)

private int compare(A, B)(auto ref A _a, auto ref B _b)
{
	static if(is(A T : ValueExpr!T))
		alias a = _a;
	else
		Floating a = _a;

	static if(is(B U : ValueExpr!U))
		alias b = _b;
	else
		Floating b = _b;

	return mpfrCompare(a.val, b.val);
}

/** "base class" for all expressions */
template Expr()
{
	auto toExpr()
	{
		return this;
	}

	enum this_is_a_mpfr_expression = true;

	auto opUnary(string op)()
	{
		return UnaryExpr!(op,typeof(this))(this);
	}

	auto opBinary(string op, B)(auto ref B b)
	{
		return BinaryExpr!(op, typeof(this), typeof(b.toExpr()))(this, b.toExpr());
	}

	bool opEquals(T)(auto ref T b) const
	{
		return compare(this, b.toExpr()) == 0;
	}

	int opCmp(T)(auto ref T b) const
	{
		return compare(this, b.toExpr());
	}
}

/** a (read-only) mpfr value or a builtin type */
struct ValueExpr(T)
{
	mixin Expr;

	T val;

	void eval(mpfr_t* dest)
	{
		mpfrUnary!"set"(dest, val);
	}
}

auto toExpr(A)(A a)
	if(isFloatingPoint!A || isSigned!A || isUnsigned!A)
{
	return ValueExpr!A(a);
}

/** unary expression, including one parameter functions */
struct UnaryExpr(string op, A)
{
	mixin Expr;

	A a;

	void eval(mpfr_t* dest)
	{
		static if(is(A T : ValueExpr!T))
		{
			mpfrUnary!op(dest, a.val);
		}
		else
		{
			a.eval(dest);
			mpfrUnary!op(dest, dest);
		}
	}
}

/** binary expression */
struct BinaryExpr(string op, A, B)
{
	mixin Expr;

	A a;
	B b;

	void eval(mpfr_t* dest)
	{
		// Value <-> Value
		static if(is(A T : ValueExpr!T) && is(B U : ValueExpr!U))
		{
			mpfrBinary!op(dest, a.val, b.val);
		}

		// Value <-> Expression
		else static if(is(A V : ValueExpr!V))
		{
			static if(is(V : const(mpfr_t)*))
			{
				if(a.val is dest)
				{
					mpfr_t tmp;
					mpfr_init2(&tmp, mpfr_get_prec(dest));
					b.eval(&tmp);
					mpfrBinary!op(dest, a.val, &tmp);
					mpfr_clear(&tmp);
				}
				else
				{
					b.eval(dest);
					mpfrBinary!op(dest, a.val, dest);
				}
			}
			else
			{
				b.eval(dest);
				mpfrBinary!op(dest, a.val, dest);
			}
		}

		// Expression <-> Value
		else static if(is(B V : ValueExpr!V))
		{
			static if(is(V : const(mpfr_t)*))
			{
				if(b.val is dest)
				{
					mpfr_t tmp;
					mpfr_init2(&tmp, mpfr_get_prec(dest));
					a.eval(&tmp);
					mpfrBinary!op(dest, &tmp, b.val);
					mpfr_clear(&tmp);
				}
				else
				{
					a.eval(dest);
					mpfrBinary!op(dest, dest, b.val);
				}
			}
			else
			{
				a.eval(dest);
				mpfrBinary!op(dest, dest, b.val);
			}
		}

		// Expression <-> Expression
		else
		{
			a.eval(dest);

			mpfr_t tmp;
			mpfr_init2(&tmp, mpfr_get_prec(dest));
			b.eval(&tmp);

			mpfrBinary!op(dest, dest, &tmp);
			mpfr_clear(&tmp);
		}
	}
}
