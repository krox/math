module math.gmp;

/**
 * binding for the GNU Multiple Precision Arithmetic Library.
 */

public import core.stdc.config : c_long, c_ulong;

/**
 * convenience wrapper for a (mutable) GMP integer.
 */
final class GmpInteger
{
	mpz_t z;

	inout(mpz_t)* ptr() inout pure nothrow @property
	{
		return &z;
	}

	this() pure nothrow
	{
		__gmpz_init(&z);
	}

	~this() pure nothrow
	{
		__gmpz_clear(&z);
	}
}


extern(C):

alias c_ulong mp_bitcnt_t;

extern immutable int __gmp_bits_per_limb;
extern immutable int __gmp_errno;
extern immutable char* __gmp_version;

private alias size_t limb;

static this()
{
	assert(limb.sizeof*8 == __gmp_bits_per_limb, "wrong gmp limb size");
}

struct mpz_t
{
	int _mp_alloc;
	int _mp_size;
	/*limb * */ size_t _mp_d;
	// having this pointer in a non-pointer variable is fine because it is never actually used in D code.
	// In particular it does never point to GC-memory.
}


size_t __gmpz_sizeinbase (const mpz_t* op , int base ) pure nothrow;

void __gmpz_init(mpz_t* x ) pure nothrow;
void __gmpz_init2(mpz_t* x ,size_t n) pure nothrow;
void __gmpz_clear(mpz_t* x) pure nothrow;
void __gmpz_realloc2(mpz_t* x, mp_bitcnt_t n) pure nothrow;	// set to 0 if it does not fit

void __gmpz_set    (mpz_t* rop, const mpz_t* op) pure nothrow;
void __gmpz_set_ui (mpz_t* rop, c_ulong op) pure nothrow;
void __gmpz_set_si (mpz_t* rop, c_long op) pure nothrow;
void __gmpz_set_d  (mpz_t* rop, double op) pure nothrow;
int  __gmpz_set_str(mpz_t* rop , const char * str , int base ) pure nothrow;	// white space allowed/ignored, if base=0, 0x/0b&/0 are recognized. returns 0 if entire string is valid number

void __gmpz_swap   (mpz_t* rop1 , mpz_t* rop2 ) pure nothrow;

void __gmpz_init_set    (mpz_t* rop, const mpz_t* op ) pure nothrow;
void __gmpz_init_set_ui (mpz_t* rop, c_ulong op) pure nothrow;
void __gmpz_init_set_si (mpz_t* rop, c_long op) pure nothrow;
void __gmpz_init_set_d  (mpz_t* rop, double op) pure nothrow;
int  __gmpz_init_set_str(mpz_t* rop, const char * str, int base) pure nothrow;

c_ulong __gmpz_get_ui (const mpz_t* op) pure nothrow;
c_long  __gmpz_get_si (const mpz_t* op) pure nothrow;
double  __gmpz_get_d  (const mpz_t* op) pure nothrow;
double  __gmpz_get_d_2exp (c_long * exp, const mpz_t* op) pure nothrow;
char*   __gmpz_get_str(char* str, int base, const mpz_t* op) pure nothrow;	// str==null or buffer of size __gmpz_sizeinbase (op, base ) + 2

void __gmpz_add       (mpz_t* rop, const mpz_t* op1, const mpz_t* op2) pure nothrow;
void __gmpz_add_ui    (mpz_t* rop, const mpz_t* op1,      c_ulong op2) pure nothrow;
void __gmpz_sub       (mpz_t* rop, const mpz_t* op1, const mpz_t* op2) pure nothrow;
void __gmpz_sub_ui    (mpz_t* rop, const mpz_t* op1,      c_ulong op2) pure nothrow;
void __gmpz_ui_sub    (mpz_t* rop,      c_ulong op1, const mpz_t* op2) pure nothrow;
void __gmpz_mul       (mpz_t* rop, const mpz_t* op1, const mpz_t* op2) pure nothrow;
void __gmpz_mul_si    (mpz_t* rop, const mpz_t* op1,       c_long op2) pure nothrow;
void __gmpz_mul_ui    (mpz_t* rop, const mpz_t* op1,      c_ulong op2) pure nothrow;
void __gmpz_addmul    (mpz_t* rop, const mpz_t* op1, const mpz_t* op2) pure nothrow;
void __gmpz_addmul_ui (mpz_t* rop, const mpz_t* op1,      c_ulong op2) pure nothrow;
void __gmpz_submul    (mpz_t* rop, const mpz_t* op1, const mpz_t* op2) pure nothrow;
void __gmpz_submul_ui (mpz_t* rop, const mpz_t* op1,      c_ulong op2) pure nothrow;

void __gmpz_mul_2exp (mpz_t* rop , const mpz_t* op1 , size_t op2 ) pure nothrow;	// rop = op1*2^op2
void __gmpz_neg (mpz_t* rop , const mpz_t* op ) pure nothrow;
void __gmpz_abs (mpz_t* rop , const mpz_t* op ) pure nothrow;

void __gmpz_cdiv_q (mpz_t* q , const mpz_t* n , const mpz_t* d ) pure nothrow;
void __gmpz_cdiv_r (mpz_t* r , const mpz_t* n , const mpz_t* d ) pure nothrow;
void __gmpz_cdiv_qr (mpz_t* q , mpz_t* r , const mpz_t* n , const mpz_t* d ) pure nothrow;

c_ulong __gmpz_cdiv_q_ui (mpz_t* q , const mpz_t* n, c_ulong d ) pure nothrow;
c_ulong __gmpz_cdiv_r_ui (mpz_t* r , const mpz_t* n, c_ulong d ) pure nothrow;
c_ulong __gmpz_cdiv_qr_ui (mpz_t* q , mpz_t* r , const mpz_t* n, c_ulong d ) pure nothrow;
c_ulong __gmpz_cdiv_ui (const mpz_t* n , c_ulong d ) pure nothrow;
void __gmpz_cdiv_q_2exp (mpz_t* q , const mpz_t* n , size_t b ) pure nothrow;
void __gmpz_cdiv_r_2exp (mpz_t* r , const mpz_t* n , size_t b ) pure nothrow;
void __gmpz_fdiv_q (mpz_t* q , const mpz_t* n , const mpz_t* d ) pure nothrow;
void __gmpz_fdiv_r (mpz_t* r , const mpz_t* n , const mpz_t* d ) pure nothrow;
void __gmpz_fdiv_qr (mpz_t* q , mpz_t* r , const mpz_t* n , const mpz_t* d ) pure nothrow;
c_ulong __gmpz_fdiv_q_ui (mpz_t* q , const mpz_t* n ,c_ulong d ) pure nothrow;
c_ulong __gmpz_fdiv_r_ui (mpz_t* r , const mpz_t* n ,c_ulong d ) pure nothrow;
c_ulong __gmpz_fdiv_qr_ui (mpz_t* q , mpz_t* r , const mpz_t* n ,c_ulong d ) pure nothrow;
c_ulong __gmpz_fdiv_ui (const mpz_t* n , c_ulong d ) pure nothrow;
void __gmpz_fdiv_q_2exp (mpz_t* q , const mpz_t* n , size_t b ) pure nothrow;
void __gmpz_fdiv_r_2exp (mpz_t* r , const mpz_t* n , size_t b ) pure nothrow;
void __gmpz_tdiv_q (mpz_t* q , const mpz_t* n , const mpz_t* d ) pure nothrow;
void __gmpz_tdiv_r (mpz_t* r , const mpz_t* n , const mpz_t* d ) pure nothrow;
void __gmpz_tdiv_qr (mpz_t* q , mpz_t* r , const mpz_t* n , const mpz_t* d ) pure nothrow;
c_ulong __gmpz_tdiv_q_ui (mpz_t* q , const mpz_t* n ,c_ulong d ) pure nothrow;
c_ulong __gmpz_tdiv_r_ui (mpz_t* r , const mpz_t* n ,c_ulong d ) pure nothrow;
c_ulong __gmpz_tdiv_qr_ui (mpz_t* q , mpz_t* r , const mpz_t* n ,c_ulong d ) pure nothrow;
c_ulong __gmpz_tdiv_ui (const mpz_t* n , c_ulong d ) pure nothrow;
void __gmpz_tdiv_q_2exp (mpz_t* q , const mpz_t* n , size_t b ) pure nothrow;
void __gmpz_tdiv_r_2exp (mpz_t* r , const mpz_t* n , size_t b ) pure nothrow;

void __gmpz_mod (mpz_t* r , const mpz_t* n , const mpz_t* d ) pure nothrow;
c_ulong __gmpz_mod_ui (mpz_t* r , const mpz_t* n , c_ulong d ) pure nothrow;

void __gmpz_divexact (mpz_t* q , const mpz_t* n , const mpz_t* d ) pure nothrow;
void __gmpz_divexact_ui (mpz_t* q , const mpz_t* n , c_ulong d ) pure nothrow;

int __gmpz_divisible_p (const mpz_t* n , const mpz_t* d ) pure nothrow;
int __gmpz_divisible_ui_p (const mpz_t* n , c_ulong d ) pure nothrow;
int __gmpz_divisible_2exp_p (const mpz_t* n , size_t b ) pure nothrow;

int __gmpz_congruent_p (const mpz_t* n , const mpz_t* c , const mpz_t* d ) pure nothrow;
int __gmpz_congruent_ui_p (const mpz_t* n , c_ulong c , c_ulong d ) pure nothrow;
int __gmpz_congruent_2exp_p (const mpz_t* n , const mpz_t* c , size_t b ) pure nothrow;

void __gmpz_powm (mpz_t* rop , const mpz_t* base , const mpz_t* exp , const mpz_t* mod ) pure nothrow;
void __gmpz_powm_ui (mpz_t* rop , const mpz_t* base , c_ulong exp , const mpz_t* mod ) pure nothrow;
void __gmpz_powm_sec (mpz_t* rop , const mpz_t* base , const mpz_t* exp , const mpz_t* mod ) pure nothrow;
void __gmpz_pow_ui (mpz_t* rop , const mpz_t* base , c_ulong exp ) pure nothrow;
void __gmpz_ui_pow_ui (mpz_t* rop , c_ulong base , c_ulong exp) pure nothrow;

int __gmpz_root (mpz_t* rop , const mpz_t* op , c_ulong n ) pure nothrow;	// returns nonzero if exact

void __gmpz_rootrem (mpz_t* root , mpz_t* rem , const mpz_t* u , c_ulong n ) pure nothrow;
void __gmpz_sqrt (mpz_t* rop , const mpz_t* op ) pure nothrow;

void __gmpz_sqrtrem (mpz_t* rop1 , mpz_t* rop2 , const mpz_t* op ) pure nothrow;

int __gmpz_perfect_power_p (const mpz_t* op ) pure nothrow;	// nonzero if perfect power
int __gmpz_perfect_square_p (const mpz_t* op ) pure nothrow;

int __gmpz_probab_prime_p (const mpz_t* n , int reps ) nothrow;	// 2=prime, 1=prob.prime, 0=composite, reps~25

void __gmpz_nextprime (mpz_t* rop , const mpz_t* op ) pure nothrow;
void __gmpz_gcd (mpz_t* rop , const mpz_t* op1 , const mpz_t* op2 ) pure nothrow;
c_ulong __gmpz_gcd_ui (mpz_t* rop , const mpz_t* op1 , c_ulong op2 ) pure nothrow;
void __gmpz_gcdext (mpz_t* g , mpz_t* s , mpz_t* t , const mpz_t* a , const mpz_t* b ) pure nothrow;
void __gmpz_lcm (mpz_t* rop , const mpz_t* op1 , const mpz_t* op2 ) pure nothrow;
void __gmpz_lcm_ui (mpz_t* rop , const mpz_t* op1 , c_ulong op2 ) pure nothrow;
int __gmpz_invert (mpz_t* rop , const mpz_t* op1 , const mpz_t* op2 ) pure nothrow;
int __gmpz_jacobi (const mpz_t* a , const mpz_t* b ) pure nothrow;
int __gmpz_legendre (const mpz_t* a , const mpz_t* p ) pure nothrow;

int __gmpz_kronecker (const mpz_t* a , const mpz_t* b ) pure nothrow;
int __gmpz_kronecker_si (const mpz_t* a , c_long b ) pure nothrow;
int __gmpz_kronecker_ui (const mpz_t* a , c_ulong b ) pure nothrow;
int __gmpz_si_kronecker (c_long a , const mpz_t* b ) pure nothrow;
int __gmpz_ui_kronecker (c_ulong a , const mpz_t* b ) pure nothrow;

mp_bitcnt_t __gmpz_remove (mpz_t* rop , const mpz_t* op , const mpz_t* f ) pure nothrow;

void __gmpz_fac_ui (mpz_t* rop , c_ulong n ) pure nothrow;
void __gmpz_2fac_ui (mpz_t* rop , c_ulong n ) pure nothrow;
void __gmpz_mfac_uiui (mpz_t* rop , c_ulong n , c_ulong m ) pure nothrow;
void __gmpz_primorial_ui (mpz_t* rop , c_ulong n ) pure nothrow;
void __gmpz_bin_ui (mpz_t* rop , const mpz_t* n , c_ulong k ) pure nothrow;
void __gmpz_bin_uiui (mpz_t* rop , c_ulong n , c_ulong k ) pure nothrow;
void __gmpz_fib_ui (mpz_t* fn , c_ulong n ) pure nothrow;
void __gmpz_fib2_ui (mpz_t* fn , mpz_t* fnsub1 , c_ulong n ) pure nothrow;
void __gmpz_lucnum_ui (mpz_t* ln , c_ulong n ) pure nothrow;
void __gmpz_lucnum2_ui (mpz_t* ln , mpz_t* lnsub1 , c_ulong n ) pure nothrow;

int __gmpz_cmp (const mpz_t* op1 , const mpz_t* op2 ) pure nothrow;
int __gmpz_cmp_d (const mpz_t* op1 , double op2 ) pure nothrow;
int __gmpz_cmp_si (const mpz_t* op1 , c_long op2 ) pure nothrow;
int __gmpz_cmp_ui (const mpz_t* op1 , c_ulong op2 ) pure nothrow;

int __gmpz_cmpabs (const mpz_t* op1 , const mpz_t* op2 ) pure nothrow;
int __gmpz_cmpabs_d (const mpz_t* op1 , double op2 ) pure nothrow;
int __gmpz_cmpabs_ui (const mpz_t* op1 , c_ulong op2 ) pure nothrow;

int __gmpz_sgn (const mpz_t* op ) pure nothrow;

void __gmpz_and (mpz_t* rop , const mpz_t* op1 , const mpz_t* op2 ) pure nothrow;
void __gmpz_ior (mpz_t* rop , const mpz_t* op1 , const mpz_t* op2 ) pure nothrow;
void __gmpz_xor (mpz_t* rop , const mpz_t* op1 , const mpz_t* op2 ) pure nothrow;
void __gmpz_com (mpz_t* rop , const mpz_t* op ) pure nothrow;

mp_bitcnt_t __gmpz_popcount (const mpz_t* op ) pure nothrow; // returns bitcnt.max on negative
mp_bitcnt_t __gmpz_hamdist (const mpz_t* op1 , const mpz_t* op2 ) pure nothrow; //returns bitcnt.max if different signs
mp_bitcnt_t __gmpz_scan0 (const mpz_t* op , size_t starting_bit ) pure nothrow;
mp_bitcnt_t __gmpz_scan1 (const mpz_t* op , size_t starting_bit ) pure nothrow;
void __gmpz_setbit (mpz_t* rop , size_t bit_index ) pure nothrow;
void __gmpz_clrbit (mpz_t* rop , size_t bit_index ) pure nothrow;
void __gmpz_combit (mpz_t* rop , size_t bit_index ) pure nothrow;
int __gmpz_tstbit (const mpz_t* op , size_t bit_index ) pure nothrow;

struct mp_randstate_t
{
  mpz_t _mp_seed;	  /* _mp_d member points to state of the generator. */
  int   _mp_alg;  /* Currently unused. */
  void* _mp_algdata; /* Pointer to function pointers structure.  */
}

void __gmp_randinit_default (mp_randstate_t*);
void __gmp_randseed_ui (mp_randstate_t*, c_ulong);
void __gmp_randclear (mp_randstate_t*);
c_ulong __gmp_urandomb_ui (mp_randstate_t*, c_ulong b); // generates b bits
c_ulong __gmp_urandomm_ui (mp_randstate_t*, c_ulong n); // generates <= n
void __gmpz_urandomb (mpz_t* rop, mp_randstate_t* state, mp_bitcnt_t b); // generates b bits
void __gmpz_urandomm (mpz_t* rop, mp_randstate_t* state, const mpz_t* n); // generates <= n
