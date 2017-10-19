module bindings.mpfr;

/**
 * Binding for the GNU MPFR library, version 3.
 */

private import core.stdc.config : c_long, c_ulong, c_long_double;
private import bindings.gmp;

extern(C):

// FIXME: not sure if these are correct on all platforms.
alias mpfr_prec_t = c_long;
alias mpfr_uprec_t = c_ulong;
alias mpfr_exp_t = c_long;
alias mpfr_uexp_t = c_ulong;
alias mpfr_sign_t = int;

enum mpfr_prec_t MPFR_PREC_MIN = 2;
enum mpfr_prec_t MPFR_PREC_MAX = mpfr_prec_t.max;
enum mpfr_exp_t MPFR_EMAX_DEFAULT = (1<<30)-1;
enum mpfr_exp_t MPFR_EMIN_DEFAULT = -MPFR_EMAX_DEFAULT;

enum MPFR_RNDN = 0;
enum MPFR_RNDZ = 1;
enum MPFR_RNDU = 2;
enum MPFR_RNDD = 3;
enum MPFR_RNDA = 4;
enum MPFR_RNDF = 5;
enum MPFR_RNDNA = -1;
alias mpfr_rnd_t = int;

enum MPFR_NAN_KIND = 0;
enum MPFR_INF_KIND = 1;
enum MPFR_ZERO_KIND = 2;
enum MPFR_REGULAR_KIND = 3;
alias mpfr_kind_t = int;

struct mpfr_t
{
	mpfr_prec_t _mpfr_prec;
	mpfr_sign_t _mpfr_sign;
	mpfr_exp_t  _mpfr_exp;
	mp_limb_t*  _mpfr_d;
}

immutable(char)* mpfr_get_version ();
immutable(char)* mpfr_get_patches ();
int mpfr_buildopt_tls_p          ();
int mpfr_buildopt_decimal_p      ();
int mpfr_buildopt_gmpinternals_p ();
immutable(char)* mpfr_buildopt_tune_case ();

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
const(char)* mpfr_print_rnd_mode (mpfr_rnd_t);

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

void mpfr_init  (mpfr_t*);
void mpfr_inits  (mpfr_t*, ...);
void mpfr_set_default_prec (mpfr_prec_t);
mpfr_prec_t mpfr_get_default_prec ();
int mpfr_init_set_str (mpfr_t*, const(char)*, int, mpfr_rnd_t);

int mpfr_check_range (mpfr_t*, int, mpfr_rnd_t);
void mpfr_free_cache ();

pure nothrow
{

void mpfr_init2 (mpfr_t*, mpfr_prec_t);
void mpfr_clear (mpfr_t*);

void mpfr_inits2 (mpfr_prec_t, mpfr_t*, ...);
void mpfr_clears (mpfr_t*, ...);

int mpfr_prec_round (mpfr_t*, mpfr_prec_t, mpfr_rnd_t);
int mpfr_can_round (const(mpfr_t)*, mpfr_exp_t, mpfr_rnd_t, mpfr_rnd_t, mpfr_prec_t);
mpfr_prec_t mpfr_min_prec (const(mpfr_t)*);

mpfr_exp_t mpfr_get_exp (const(mpfr_t)*);
int mpfr_set_exp (mpfr_t*, mpfr_exp_t);
mpfr_prec_t mpfr_get_prec (const(mpfr_t)*);
void mpfr_set_prec (mpfr_t*, mpfr_prec_t);
void mpfr_set_prec_raw (mpfr_t*, mpfr_prec_t);

int mpfr_set_d (mpfr_t*, double, mpfr_rnd_t);
int mpfr_set_flt (mpfr_t*, float, mpfr_rnd_t);

int mpfr_set_ld (mpfr_t*, c_long_double, mpfr_rnd_t);
int mpfr_set_z (mpfr_t*, const(mpz_t)*, mpfr_rnd_t);
int mpfr_set_z_2exp (mpfr_t*, const(mpz_t)*, mpfr_exp_t, mpfr_rnd_t);
void mpfr_set_nan (mpfr_t*);
void mpfr_set_inf (mpfr_t*, int);
void mpfr_set_zero (mpfr_t*, int);
int mpfr_set_si (mpfr_t*, c_long, mpfr_rnd_t);
int mpfr_set_ui (mpfr_t*, c_ulong, mpfr_rnd_t);
int mpfr_set_si_2exp (mpfr_t*, c_long, mpfr_exp_t, mpfr_rnd_t);
int mpfr_set_ui_2exp (mpfr_t*,c_ulong,mpfr_exp_t,mpfr_rnd_t);
//int mpfr_set_q (mpfr_t*, mpq_srcptr, mpfr_rnd_t);
int mpfr_set_str (mpfr_t*, const(char)*, int, mpfr_rnd_t);
int mpfr_set4 (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t, int);
int mpfr_abs (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_set (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_neg (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_signbit (const(mpfr_t)*);
int mpfr_setsign (mpfr_t*, const(mpfr_t)*, int, mpfr_rnd_t);
int mpfr_copysign (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);

mpfr_exp_t mpfr_get_z_2exp (mpz_t*, const(mpfr_t)*);
float mpfr_get_flt (const(mpfr_t)*, mpfr_rnd_t);
double mpfr_get_d (const(mpfr_t)*, mpfr_rnd_t);

c_long_double mpfr_get_ld (const(mpfr_t)*, mpfr_rnd_t);
double mpfr_get_d1 (const(mpfr_t)*);
double mpfr_get_d_2exp (c_long*, const(mpfr_t)*, mpfr_rnd_t);
c_long_double mpfr_get_ld_2exp (c_long*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_frexp (mpfr_exp_t*, mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
c_long mpfr_get_si (const(mpfr_t)*, mpfr_rnd_t);
c_ulong mpfr_get_ui (const(mpfr_t)*, mpfr_rnd_t);
char*mpfr_get_str (char*, mpfr_exp_t*, int, size_t, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_get_z (mpz_t* z, const(mpfr_t)* f, mpfr_rnd_t);

void mpfr_free_str (char *);

int mpfr_urandom (mpfr_t*, mp_randstate_t*, mpfr_rnd_t);
int mpfr_grandom (mpfr_t*, mpfr_t*, mp_randstate_t*, mpfr_rnd_t);
int mpfr_urandomb (mpfr_t*, mp_randstate_t*);

void mpfr_nextabove (mpfr_t*);
void mpfr_nextbelow (mpfr_t*);
void mpfr_nexttoward (mpfr_t*, const(mpfr_t)*);

int mpfr_printf (const char*, ...);
int mpfr_asprintf (char**, const char*, ...);
int mpfr_sprintf (char*, const char*, ...);
int mpfr_snprintf (char*, size_t, const char*, ...);

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
int mpfr_si_sub (mpfr_t*, c_long, const(mpfr_t)*, mpfr_rnd_t);
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

int mpfr_log (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_log2 (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_log10 (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_log1p (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_exp (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_exp2 (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_exp10 (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_expm1 (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_eint (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_li2 (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);

int mpfr_cmp  (const(mpfr_t)*, const(mpfr_t)*);
int mpfr_cmp3 (const(mpfr_t)*, const(mpfr_t)*, int);
int mpfr_cmp_d (const(mpfr_t)*, double);
int mpfr_cmp_ld (const(mpfr_t)*, c_long_double);
int mpfr_cmpabs (const(mpfr_t)*, const(mpfr_t)*);
int mpfr_cmp_ui (const(mpfr_t)*, c_ulong);
int mpfr_cmp_si (const(mpfr_t)*, c_long);
int mpfr_cmp_ui_2exp (const(mpfr_t)*, c_ulong, mpfr_exp_t);
int mpfr_cmp_si_2exp (const(mpfr_t)*, c_long, mpfr_exp_t);
void mpfr_reldiff (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_eq (const(mpfr_t)*, const(mpfr_t)*, c_ulong);
int mpfr_sgn (const(mpfr_t)*);

int mpfr_mul_2exp (mpfr_t*, const(mpfr_t)*, c_ulong, mpfr_rnd_t);
int mpfr_div_2exp (mpfr_t*, const(mpfr_t)*, c_ulong, mpfr_rnd_t);
int mpfr_mul_2ui (mpfr_t*, const(mpfr_t)*, c_ulong, mpfr_rnd_t);
int mpfr_div_2ui (mpfr_t*, const(mpfr_t)*, c_ulong, mpfr_rnd_t);
int mpfr_mul_2si (mpfr_t*, const(mpfr_t)*, c_long, mpfr_rnd_t);
int mpfr_div_2si (mpfr_t*, const(mpfr_t)*, c_long, mpfr_rnd_t);

int mpfr_rint (mpfr_t*,const(mpfr_t)*, mpfr_rnd_t);
int mpfr_round (mpfr_t*, const(mpfr_t)*);
int mpfr_trunc (mpfr_t*, const(mpfr_t)*);
int mpfr_ceil (mpfr_t*, const(mpfr_t)*);
int mpfr_floor (mpfr_t*, const(mpfr_t)*);
int mpfr_rint_round (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_rint_trunc (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_rint_ceil (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_rint_floor (mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_frac (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_modf (mpfr_t*, mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_remquo (mpfr_t*, c_long*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_remainder (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fmod (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_fits_ulong_p (const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fits_slong_p (const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fits_uint_p (const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fits_sint_p (const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fits_ushort_p (const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fits_sshort_p (const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fits_uintmax_p (const(mpfr_t)*,mpfr_rnd_t);
int mpfr_fits_intmax_p (const(mpfr_t)*, mpfr_rnd_t);

void mpfr_extract (mpz_t*, const(mpfr_t)*, uint);
void mpfr_swap (mpfr_t*, mpfr_t*);
void mpfr_dump (const(mpfr_t)*);

int mpfr_nan_p (const(mpfr_t)*);
int mpfr_inf_p (const(mpfr_t)*);
int mpfr_number_p (const(mpfr_t)*);
int integer_p (const(mpfr_t)*);
int mpfr_zero_p (const(mpfr_t)*);
int mpfr_regular_p (const(mpfr_t)*);

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
int mpfr_sinh_cosh (mpfr_t*, mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);

int mpfr_sech (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_csch (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_coth (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);

int mpfr_acos (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_asin (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_atan (mpfr_t*,const(mpfr_t)*,mpfr_rnd_t);
int mpfr_sin (mpfr_t*, const(mpfr_t)*,mpfr_rnd_t);
int mpfr_sin_cos (mpfr_t*, mpfr_t*, const(mpfr_t)*, mpfr_rnd_t);
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

/*
int mpfr_mul_q (mpfr_t*, const(mpfr_t)*, mpq_srcptr, mpfr_rnd_t);
int mpfr_div_q (mpfr_t*, const(mpfr_t)*, mpq_srcptr, mpfr_rnd_t);
int mpfr_add_q (mpfr_t*, const(mpfr_t)*, mpq_srcptr, mpfr_rnd_t);
int mpfr_sub_q (mpfr_t*, const(mpfr_t)*, mpq_srcptr, mpfr_rnd_t);
int mpfr_cmp_q (const(mpfr_t)*, mpq_srcptr);
*/

int mpfr_fma (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_fms (mpfr_t*, const(mpfr_t)*, const(mpfr_t)*, const(mpfr_t)*, mpfr_rnd_t);
int mpfr_sum (mpfr_t*, const(mpfr_t*)*, c_ulong, mpfr_rnd_t);

int  mpfr_subnormalize (mpfr_t*, int, mpfr_rnd_t);

int  mpfr_strtofr (mpfr_t*, const(char)*, char **, int, mpfr_rnd_t);

size_t mpfr_custom_get_size   (mpfr_prec_t);
void   mpfr_custom_init    (void *, mpfr_prec_t);
void * mpfr_custom_get_significand (const(mpfr_t)*);
mpfr_exp_t mpfr_custom_get_exp  (const(mpfr_t)*);
void   mpfr_custom_move       (mpfr_t*, void *);
void   mpfr_custom_init_set   (mpfr_t*, int, mpfr_exp_t, mpfr_prec_t, void *);
int    mpfr_custom_get_kind   (const(mpfr_t)*);

}
