module math.mpfr;

/**
 * binding for the GNU MPFR library.
 */

public import core.stdc.config : c_long, c_ulong;
private import std.traits;
private import math.gmp;


/**
 * overloaded template version of mpfr_*(dest, a, b, rnd)
 * a and a can be any combination of mpfr_t* / int / uint / double / ..
 */
int mpfrUnary(string op, A)(mpfr_t* dest, const A a, mpfr_rnd_t rnd = MPFR_RNDN)
{
	     static if(op == "-") enum opName = "neg";
	else static if(op == "+") enum opName = "set";
	else enum opName = op;

	// mpfr
	static if(is(A : const(mpfr_t)*))
		enum funName = opName;

	// builtin
	else static if(isFloatingPoint!A)
		enum funName = opName~"_d";
	else static if(isUnsigned!A)
		enum funName = opName~"_ui";
	else static if(isSigned!A)
		enum funName = opName~"_si";

	else static assert(false);

	return mixin("mpfr_"~funName)(dest, a, rnd);
}

/** ditto */
int mpfrBinary(string op, A, B)(mpfr_t* dest, const A a, const B b, mpfr_rnd_t rnd = MPFR_RNDN)
{
	     static if(op == "+") enum opName = "add";
	else static if(op == "-") enum opName = "sub";
	else static if(op == "*") enum opName = "mul";
	else static if(op == "/") enum opName = "div";
	else enum opName = op;

	// mpfr <-> mpfr
	static if(is(A : const(mpfr_t)*) && is(B : const(mpfr_t)*))
		enum funName = opName;

	// mpfr <-> builtin
	else static if(is(A : const(mpfr_t)*) && isFloatingPoint!B)
		enum funName = opName~"_d";
	else static if(is(A : const(mpfr_t)*) && isUnsigned!B)
		enum funName = opName~"_ui";
	else static if(is(A : const(mpfr_t)*) && isSigned!B)
		enum funName = opName~"_si";

	// builtin <-> mpfr
	else static if(isFloatingPoint!A && is(B : const(mpfr_t)*))
		enum funName = "d_"~opName;
	else static if(isUnsigned!A && is(B : const(mpfr_t)*))
		enum funName = "ui_"~opName;
	else static if(isSigned!A && is(B : const(mpfr_t)*))
		enum funName = "si_"~opName;

	else static assert(false);

	return mixin("mpfr_"~funName)(dest, a, b, rnd);
}

/** ditto */
int mpfrCompare(A, B)(const A a, const B b, mpfr_rnd_t rnd = MPFR_RNDN)
{
	// mpfr <-> mpfr
	static if(is(A : const(mpfr_t)*) && is(B : const(mpfr_t)*))
		return mpfr_cmp(a, b);

	// mpfr <-> builtin
	else static if(is(A : const(mpfr_t)*) && isFloatingPoint!B)
		return mpfr_cmp_d(a, b);
	else static if(is(A : const(mpfr_t)*) && isUnsigned!B)
		return mpfr_cmp_ui(a, b);
	else static if(is(A : const(mpfr_t)*) && isSigned!B)
		return mpfr_cmp_si(a, b);

	// builtin <-> mpfr
	else static if(isFloatingPoint!A && is(B : const(mpfr_t)*))
		return -mpfr_cmp_d(b, a);
	else static if(isUnsigned!A && is(B : const(mpfr_t)*))
		return -mpfr_cmp_ui(b, a);
	else static if(isSigned!A && is(B : const(mpfr_t)*))
		return -mpfr_cmp_si(b, a);

	else static assert(false);
}

/**
 * Functions not present in MPFR because they are redunant by commutativity.
 * But having them simplifies some (template) code, so here they go.
 */
int mpfr_si_add (mpfr_t* dst, c_long a, const(mpfr_t)* b, mpfr_rnd_t rnd)
	{ return mpfr_add_si(dst, b, a, rnd); }
int mpfr_si_mul (mpfr_t* dst, c_long a, const(mpfr_t)* b, mpfr_rnd_t rnd)
	{ return mpfr_mul_si(dst, b, a, rnd); }
int mpfr_ui_add (mpfr_t* dst, c_ulong a, const(mpfr_t)* b, mpfr_rnd_t rnd)
	{ return mpfr_add_ui(dst, b, a, rnd); }
int mpfr_ui_mul (mpfr_t* dst, c_ulong a, const(mpfr_t)* b, mpfr_rnd_t rnd)
	{ return mpfr_mul_ui(dst, b, a, rnd); }
int mpfr_d_add (mpfr_t* dst, double a, const(mpfr_t)* b, mpfr_rnd_t rnd)
	{ return mpfr_add_d(dst, b, a, rnd); }
int mpfr_d_mul (mpfr_t* dst, double a, const(mpfr_t)* b, mpfr_rnd_t rnd)
	{ return mpfr_mul_d(dst, b, a, rnd); }

extern(C):

/* Define MPFR version number */
enum MPFR_VERSION_MAJOR = 3;
enum MPFR_VERSION_MINOR = 1;
enum MPFR_VERSION_PATCHLEVEL = 2;
enum MPFR_VERSION_STRING = "3.1.2-p11";


/* Definition of rounding modes (DON'T USE MPFR_RNDNA!).
   Warning! Changing the contents of this enum should be seen as an
   interface change since the old and the new types are not compatible
   (the integer type compatible with the enumerated type can even change,
   see ISO C99, 6.7.2.2#4), and in Makefile.am, AGE should be set to 0.

   MPFR_RNDU must appear just before MPFR_RNDD (see
   MPFR_IS_RNDUTEST_OR_RNDDNOTTEST in mpfr-impl.h).

   MPFR_RNDF has been added, though not implemented yet, in order to avoid
   to break the ABI once faithful rounding gets implemented.

   If you change the order of the rounding modes, please update the routines
   in texceptions.c which assume 0=RNDN, 1=RNDZ, 2=RNDU, 3=RNDD, 4=RNDA.
*/
enum
{
    MPFR_RNDN=0,  /* round to nearest, with ties to even */
    MPFR_RNDZ,    /* round toward zero */
    MPFR_RNDU,    /* round toward +Inf */
    MPFR_RNDD,    /* round toward -Inf */
    MPFR_RNDA,    /* round away from zero */
    MPFR_RNDF,    /* faithful rounding (not implemented yet) */
    MPFR_RNDNA=-1 /* round to nearest, with ties away from zero (mpfr_round) */
}

alias mpfr_rnd_t = int;
alias mpfr_exp_t = c_long; // might be c_int on certain hardware (see mpfr.h _MPFR_EXP/PREC_FORMAT)
alias mpfr_prec_t = c_long; // ditto
alias mpfr_sign_t = int;

static assert(mpfr_prec_t.sizeof <= mpfr_exp_t.sizeof);

/* Definition of the main structure */
struct mpfr_t
{
  mpfr_prec_t  _mpfr_prec;
  mpfr_sign_t  _mpfr_sign;
  mpfr_exp_t   _mpfr_exp;

  //mp_limb_t* _mpfr_d;
  size_t mpfr_d;
  // having this pointer in a non-pointer variable is fine because it is never actually used in D code.
  // In particular it does never point to GC-memory.
}

/* DON'T USE THIS! (For MPFR-public macros only, see below.)
   The mpfr_sgn macro uses the fact that __MPFR_EXP_NAN and __MPFR_EXP_ZERO
   are the smallest values. */
enum __MPFR_EXP_MAX = mpfr_exp_t.max;
enum __MPFR_EXP_NAN = 1 - __MPFR_EXP_MAX;
enum __MPFR_EXP_ZERO = 0 - __MPFR_EXP_MAX;
enum __MPFR_EXP_INF  = 2 - __MPFR_EXP_MAX;



/* For those who need a direct and fast access to the sign field.
   However it is not in the API, thus use it at your own risk: it might
   not be supported, or change name, in further versions!
   Unfortunately, it must be defined here (instead of MPFR's internal
   header file mpfr-impl.h) because it is used by some macros below.
*/
//#define MPFR_SIGN(x) ( (x)->_mpfr_sign)

/* Stack interface */
enum
{
  MPFR_NAN_KIND = 0,
  MPFR_INF_KIND = 1,
  MPFR_ZERO_KIND = 2,
  MPFR_REGULAR_KIND = 3
}


immutable(char)* mpfr_get_version();
immutable(char)* mpfr_get_patches();
int mpfr_buildopt_tls_p();
int mpfr_buildopt_decimal_p();
int mpfr_buildopt_gmpinternals_p();
immutable(char)* mpfr_buildopt_tune_case();

mpfr_exp_t mpfr_get_emin     ();
int        mpfr_set_emin     (mpfr_exp_t);
mpfr_exp_t mpfr_get_emin_min ();
mpfr_exp_t mpfr_get_emin_max ();
mpfr_exp_t mpfr_get_emax     ();
int        mpfr_set_emax     (mpfr_exp_t);
mpfr_exp_t mpfr_get_emax_min ();
mpfr_exp_t mpfr_get_emax_max ();

void mpfr_set_default_rounding_mode (mpfr_rnd_t);
mpfr_rnd_t mpfr_get_default_rounding_mode ();
const(char) * mpfr_print_rnd_mode (mpfr_rnd_t);

void mpfr_clear_flags ();
void mpfr_clear_underflow ();
void mpfr_clear_overflow ();
void mpfr_clear_divby0 ();
void mpfr_clear_nanflag ();
void mpfr_clear_inexflag ();
void mpfr_clear_erangeflag ();

void mpfr_set_underflow ();
void mpfr_set_overflow ();
void mpfr_set_divby0 ();
void mpfr_set_nanflag ();
void mpfr_set_inexflag ();
void mpfr_set_erangeflag ();

int mpfr_underflow_p ();
int mpfr_overflow_p ();
int mpfr_divby0_p ();
int mpfr_nanflag_p ();
int mpfr_inexflag_p ();
int mpfr_erangeflag_p ();

int mpfr_check_range (mpfr_t*, int, mpfr_rnd_t);

void mpfr_init2 (mpfr_t*, mpfr_prec_t);
void mpfr_init (mpfr_t*);
void mpfr_clear (mpfr_t*);

/*void mpfr_inits2 (mpfr_prec_t, mpfr_t*, ...);
void mpfr_inits (mpfr_t*, ...);
void mpfr_clears (mpfr_t*, ...);*/

int mpfr_prec_round (mpfr_t*, mpfr_prec_t, mpfr_rnd_t);
int mpfr_can_round (const(mpfr_t)*, mpfr_exp_t, mpfr_rnd_t, mpfr_rnd_t, mpfr_prec_t);
mpfr_prec_t mpfr_min_prec (const(mpfr_t)*);

mpfr_exp_t mpfr_get_exp (const(mpfr_t)* x) { return x._mpfr_exp; }
int mpfr_set_exp (mpfr_t*, mpfr_exp_t);
mpfr_prec_t mpfr_get_prec (const(mpfr_t)* x) { return x._mpfr_prec; }
void mpfr_set_prec (mpfr_t*, mpfr_prec_t);
void mpfr_set_prec_raw (mpfr_t*, mpfr_prec_t);
void mpfr_set_default_prec (mpfr_prec_t);
mpfr_prec_t mpfr_get_default_prec ();

int mpfr_set_d (mpfr_t*, double, mpfr_rnd_t);
int mpfr_set_flt (mpfr_t*, float, mpfr_rnd_t);
/*#ifdef MPFR_WANT_DECIMAL_FLOATS
int mpfr_set_decimal64 (mpfr_t*, _Decimal64, mpfr_rnd_t);
#endif*/
//int mpfr_set_ld (mpfr_t*, long double, mpfr_rnd_t);
int mpfr_set_z (mpfr_t*, const(mpz_t)*, mpfr_rnd_t);
int mpfr_set_z_2exp (mpfr_t*, const(mpz_t)*, mpfr_exp_t, mpfr_rnd_t);
void mpfr_set_nan (mpfr_t*);
void mpfr_set_inf (mpfr_t*, int);
void mpfr_set_zero (mpfr_t*, int);
//int mpfr_set_f (mpfr_t*, mpf_srcptr, mpfr_rnd_t);
//int mpfr_get_f (mpf_ptr, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_set_si (mpfr_t*, c_long, mpfr_rnd_t);
int mpfr_set_ui (mpfr_t*, c_ulong, mpfr_rnd_t);
int mpfr_set_si_2exp (mpfr_t*, c_long, mpfr_exp_t, mpfr_rnd_t);
int mpfr_set_ui_2exp (mpfr_t*,c_ulong,mpfr_exp_t,mpfr_rnd_t);
//int mpfr_set_q (mpfr_t*, mpq_srcptr, mpfr_rnd_t);
int mpfr_set_str (mpfr_t*, const char *, int, mpfr_rnd_t);
int mpfr_init_set_str (mpfr_t*, const char *, int, mpfr_rnd_t);
int mpfr_set4 (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t, int);
int mpfr_abs (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_set (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_neg (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_signbit (const(mpfr_t)*);
int mpfr_setsign (mpfr_t*, const(mpfr_t)*, int, mpfr_rnd_t);
int mpfr_copysign (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);

mpfr_exp_t mpfr_get_z_2exp (const(mpz_t)*, const(mpfr_t)*);
float mpfr_get_flt (const(mpfr_t)*, mpfr_rnd_t);
double mpfr_get_d (const(mpfr_t)*, mpfr_rnd_t);
/*#ifdef MPFR_WANT_DECIMAL_FLOATS
_Decimal64 mpfr_get_decimal64 (const(mpfr_t)*, mpfr_rnd_t);
#endif*/
//long double mpfr_get_ld (const(mpfr_t)*, mpfr_rnd_t);
double mpfr_get_d1 (const(mpfr_t)*);
double mpfr_get_d_2exp (c_long*, const(mpfr_t)*, mpfr_rnd_t);
//long double mpfr_get_ld_2exp (c_long*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_frexp (mpfr_exp_t*, mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
c_long mpfr_get_si (const(mpfr_t)*, mpfr_rnd_t);
c_ulong mpfr_get_ui (const(mpfr_t)*, mpfr_rnd_t);
char*mpfr_get_str (char*, mpfr_exp_t*, int, size_t, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_get_z (const(mpz_t)* z, const(mpfr_t)* f, mpfr_rnd_t);

void mpfr_free_str (char *);

int mpfr_urandom (mpfr_t*, mp_randstate_t*, mpfr_rnd_t);
int mpfr_grandom (mpfr_t*, mpfr_t*, mp_randstate_t*, mpfr_rnd_t);
int mpfr_urandomb (mpfr_t*, mp_randstate_t*);

void mpfr_nextabove (mpfr_t*);
void mpfr_nextbelow (mpfr_t*);
void mpfr_nexttoward (mpfr_t*, const(mpfr_t)*);

// actually these take "..." as last argument, but I dont think I want to use that
int mpfr_printf (const char*, const(mpfr_t)*);
int mpfr_asprintf (char**, const char*, const(mpfr_t)*);
int mpfr_sprintf (char*, const char*, const(mpfr_t)*);
int mpfr_snprintf (char*, size_t, const char*, const(mpfr_t)*);

int mpfr_pow (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_pow_si (mpfr_t*, const(mpfr_t)*, c_long, mpfr_rnd_t);
int mpfr_pow_ui (mpfr_t*, const(mpfr_t)*, c_ulong, mpfr_rnd_t);
int mpfr_ui_pow_ui (mpfr_t*, c_ulong, c_ulong, mpfr_rnd_t);
int mpfr_ui_pow (mpfr_t*, c_ulong, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_pow_z (mpfr_t*, const(mpfr_t)*, const(mpz_t)*, mpfr_rnd_t);

int mpfr_sqrt (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_sqrt_ui (mpfr_t*, c_ulong, mpfr_rnd_t);
int mpfr_rec_sqrt (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_add (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_sub (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_mul (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_div (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_add_ui (mpfr_t*, const(mpfr_t)*, c_ulong, mpfr_rnd_t);
int mpfr_sub_ui (mpfr_t*, const(mpfr_t)*, c_ulong, mpfr_rnd_t);
int mpfr_ui_sub (mpfr_t*, c_ulong, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_mul_ui (mpfr_t*, const(mpfr_t)*, c_ulong, mpfr_rnd_t);
int mpfr_div_ui (mpfr_t*, const(mpfr_t)*, c_ulong, mpfr_rnd_t);
int mpfr_ui_div (mpfr_t*, c_ulong, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_add_si (mpfr_t*, const(mpfr_t)*, c_long, mpfr_rnd_t);
int mpfr_sub_si (mpfr_t*, const(mpfr_t)*, c_long, mpfr_rnd_t);
int mpfr_si_sub (mpfr_t*, c_long,  const(mpfr_t)*, mpfr_rnd_t);
int mpfr_mul_si (mpfr_t*, const(mpfr_t)*, c_long, mpfr_rnd_t);
int mpfr_div_si (mpfr_t*, const(mpfr_t)*, c_long, mpfr_rnd_t);
int mpfr_si_div (mpfr_t*, c_long, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_add_d (mpfr_t*, const(mpfr_t)*, double, mpfr_rnd_t);
int mpfr_sub_d (mpfr_t*, const(mpfr_t)*, double, mpfr_rnd_t);
int mpfr_d_sub (mpfr_t*, double, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_mul_d (mpfr_t*, const(mpfr_t)*, double, mpfr_rnd_t);
int mpfr_div_d (mpfr_t*, const(mpfr_t)*, double, mpfr_rnd_t);
int mpfr_d_div (mpfr_t*, double, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_sqr (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);

int mpfr_const_pi (mpfr_t*, mpfr_rnd_t);
int mpfr_const_log2 (mpfr_t*, mpfr_rnd_t);
int mpfr_const_euler (mpfr_t*, mpfr_rnd_t);
int mpfr_const_catalan (mpfr_t*, mpfr_rnd_t);

int mpfr_agm (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_log   (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_log2  (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_log10 (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_log1p (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_exp   (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_exp2  (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_exp10 (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_expm1 (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_eint  (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_li2   (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);

int mpfr_cmp    (const(mpfr_t)*, const(mpfr_t)*);
int mpfr_cmp3   (const(mpfr_t)*, const(mpfr_t)*, int);
int mpfr_cmp_d  (const(mpfr_t)*, double);
//int mpfr_cmp_ld (const(mpfr_t)*, long double);
int mpfr_cmpabs (const(mpfr_t)*, const(mpfr_t)*);
int mpfr_cmp_ui (const(mpfr_t)*, c_ulong);
int mpfr_cmp_si (const(mpfr_t)*, c_long);
int mpfr_cmp_ui_2exp (const(mpfr_t)*, c_ulong, mpfr_exp_t);
int mpfr_cmp_si_2exp (const(mpfr_t)*, c_long, mpfr_exp_t);
void mpfr_reldiff (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_eq  (const(mpfr_t)*, const(mpfr_t)*, c_ulong);
int mpfr_sgn (const(mpfr_t)*);
/*#define mpfr_sgn(_x)                                               \
  ( (_x)->_mpfr_exp < __MPFR_EXP_INF ?                              \
   (mpfr_nan_p (_x) ? mpfr_set_erangeflag () : (mpfr_void) 0), 0 : \
   MPFR_SIGN (_x) )*/

int mpfr_mul_2exp (mpfr_t*, const(mpfr_t)*, c_ulong, mpfr_rnd_t);
int mpfr_div_2exp (mpfr_t*, const(mpfr_t)*, c_ulong, mpfr_rnd_t);
int mpfr_mul_2ui  (mpfr_t*, const(mpfr_t)*, c_ulong, mpfr_rnd_t);
int mpfr_div_2ui  (mpfr_t*, const(mpfr_t)*, c_ulong, mpfr_rnd_t);
int mpfr_mul_2si  (mpfr_t*, const(mpfr_t)*, c_long, mpfr_rnd_t);
int mpfr_div_2si  (mpfr_t*, const(mpfr_t)*, c_long, mpfr_rnd_t);

int mpfr_rint  (mpfr_t*,const(mpfr_t)*, mpfr_rnd_t);
int mpfr_round (mpfr_t* x, const(mpfr_t)* y) { return mpfr_rint(x, y, MPFR_RNDNA); }
int mpfr_trunc (mpfr_t* x, const(mpfr_t)* y) { return mpfr_rint(x, y, MPFR_RNDZ); }
int mpfr_ceil  (mpfr_t* x, const(mpfr_t)* y) { return mpfr_rint(x, y, MPFR_RNDU); }
int mpfr_floor (mpfr_t* x, const(mpfr_t)* y) { return mpfr_rint(x, y, MPFR_RNDD); }
int mpfr_rint_round (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_rint_trunc (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_rint_ceil  (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_rint_floor (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_frac (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_modf (mpfr_t*, mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_remquo (mpfr_t*, c_long*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_remainder (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fmod (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_fits_uc_long_p (const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fits_sc_long_p (const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fits_uint_p    (const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fits_sint_p    (const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fits_ushort_p  (const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fits_sshort_p  (const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fits_uintmax_p (const(mpfr_t)*,mpfr_rnd_t);
int mpfr_fits_intmax_p  (const(mpfr_t)*, mpfr_rnd_t);

void mpfr_extract (const(mpz_t)*, const(mpfr_t)*, uint);
void mpfr_swap (mpfr_t*, mpfr_t*);
void mpfr_dump (const(mpfr_t)*);

int mpfr_nan_p (const(mpfr_t)* x) { return x._mpfr_exp == __MPFR_EXP_NAN; }
int mpfr_inf_p (const(mpfr_t)* x) { return x._mpfr_exp == __MPFR_EXP_INF; }
int mpfr_number_p (const(mpfr_t)*);
int mpfr_integer_p (const(mpfr_t)*);
int mpfr_zero_p (const(mpfr_t)* x) { return x._mpfr_exp == __MPFR_EXP_ZERO; }
int mpfr_regular_p (const(mpfr_t)* x) { return x._mpfr_exp > __MPFR_EXP_INF; }

int mpfr_greater_p (const(mpfr_t)*, const(mpfr_t)*);
int mpfr_greaterequal_p (const(mpfr_t)*, const(mpfr_t)*);
int mpfr_less_p (const(mpfr_t)*, const(mpfr_t)*);
int mpfr_lessequal_p (const(mpfr_t)*, const(mpfr_t)*);
int mpfr_lessgreater_p (const(mpfr_t)*,const(mpfr_t)*);
int mpfr_equal_p (const(mpfr_t)*, const(mpfr_t)*);
int mpfr_unordered_p (const(mpfr_t)*, const(mpfr_t)*);

int mpfr_atanh (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_acosh (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_asinh (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_cosh (mpfr_t*,const(mpfr_t)*, mpfr_rnd_t);
int mpfr_sinh (mpfr_t*,const(mpfr_t)*, mpfr_rnd_t);
int mpfr_tanh (mpfr_t*,const(mpfr_t)*, mpfr_rnd_t);
int mpfr_sinh_cosh (mpfr_t*, mpfr_t*,  const(mpfr_t)*, mpfr_rnd_t);

int mpfr_sech (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_csch (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_coth (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);

int mpfr_acos (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_asin (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_atan (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_sin (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_sin_cos (mpfr_t*, mpfr_t*,  const(mpfr_t)*, mpfr_rnd_t);
int mpfr_cos (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_tan (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_atan2 (mpfr_t*,const(mpfr_t)*,const(mpfr_t)*, mpfr_rnd_t);
int mpfr_sec (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_csc (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_cot (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);

int mpfr_hypot (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_erf (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_erfc (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_cbrt (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_root (mpfr_t*,const(mpfr_t)*,c_ulong,mpfr_rnd_t);
int mpfr_gamma (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_lngamma (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_lgamma (mpfr_t*,int*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_digamma (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_zeta (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_zeta_ui (mpfr_t*,c_ulong,mpfr_rnd_t);
int mpfr_fac_ui (mpfr_t*, c_ulong, mpfr_rnd_t);
int mpfr_j0 (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_j1 (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_jn (mpfr_t*, c_long, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_y0 (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_y1 (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_yn (mpfr_t*, c_long, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_ai (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_min (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_max (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_dim (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_mul_z (mpfr_t*, const(mpfr_t)*, const(mpz_t)*, mpfr_rnd_t);
int mpfr_div_z (mpfr_t*, const(mpfr_t)*, const(mpz_t)*, mpfr_rnd_t);
int mpfr_add_z (mpfr_t*, const(mpfr_t)*, const(mpz_t)*, mpfr_rnd_t);
int mpfr_sub_z (mpfr_t*, const(mpfr_t)*, const(mpz_t)*, mpfr_rnd_t);
int mpfr_z_sub (mpfr_t*, const(mpz_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_cmp_z (const(mpfr_t)*, const(mpz_t)*);

//int mpfr_mul_q (mpfr_t*, const(mpfr_t)*, mpq_srcptr, mpfr_rnd_t);
//int mpfr_div_q (mpfr_t*, const(mpfr_t)*, mpq_srcptr, mpfr_rnd_t);
//int mpfr_add_q (mpfr_t*, const(mpfr_t)*, mpq_srcptr, mpfr_rnd_t);
//int mpfr_sub_q (mpfr_t*, const(mpfr_t)*, mpq_srcptr, mpfr_rnd_t);
//int mpfr_cmp_q (const(mpfr_t)*, mpq_srcptr);

//int mpfr_cmp_f (const(mpfr_t)*, mpf_srcptr);

int mpfr_fma (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fms (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_sum (mpfr_t*, const(mpfr_t*)*, c_ulong, mpfr_rnd_t);

void mpfr_free_cache ();

int  mpfr_subnormalize (mpfr_t*, int, mpfr_rnd_t);

int  mpfr_strtofr (mpfr_t*, const char *,   char **, int, mpfr_rnd_t);


/*#define mpfr_cmp_ui(b,i) mpfr_cmp_ui_2exp(b,i,0)
#define mpfr_cmp_si(b,i) mpfr_cmp_si_2exp(b,i,0)
#define mpfr_set(a,b,r)  mpfr_set4(a,b,r,MPFR_SIGN(b) )
#define mpfr_abs(a,b,r)  mpfr_set4(a,b,r,1)
#define mpfr_copysign(a,b,c,r) mpfr_set4(a,b,r,MPFR_SIGN(c) )
#define mpfr_setsign(a,b,s,r) mpfr_set4(a,b,r,(s) ? -1 : 1)
#define mpfr_signbit(x)  (MPFR_SIGN(x) < 0)
#define mpfr_cmp(b, c)   mpfr_cmp3(b, c, 1)
#define mpfr_mul_2exp(y,x,n,r) mpfr_mul_2ui( (y),(x),(n),(r) )
#define mpfr_div_2exp(y,x,n,r) mpfr_div_2ui( (y),(x),(n),(r) )*/
