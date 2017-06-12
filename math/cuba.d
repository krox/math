/*
	cuba.h
		Prototypes for the Cuba library
		this file is part of Cuba
		last modified 28 Nov 14 th
*/

module math.cuba;
extern(C):

alias integrand_t = int function(const(int)* ndim, const(double)* x, const(int)* ncomp, double* f, void* userdata);
alias peakfinder_t = void function(const(int)* ndim, const(double)* b, int* n, double* x, void* userdata);

void llVegas(int ndim, int ncomp,
  integrand_t integrand, void* userdata, long nvec,
  double epsrel, double epsabs,
  int flags, int seed,
  long mineval, long maxeval,
  long nstart, long nincrease,
  long nbatch,
  int gridno, const(char)* statefile, void* spin,
  long* neval, int* fail,
  double* integral, double* error, double* prob);

void llSuave(int ndim, int ncomp,
  integrand_t integrand, void* userdata, long nvec,
  double epsrel, double epsabs,
  int flags, int seed,
  long mineval, long maxeval,
  long nnew, long nmin,
  double flatness, const(char)* statefile, void* spin,
  int *nregions, long *neval, int *fail,
  double* integral, double* error, double* prob);

void llDivonne(int ndim, int ncomp,
  integrand_t integrand, void* userdata, long nvec,
  double epsrel, double epsabs,
  int flags, int seed,
  long mineval, long maxeval,
  int key1, int key2, int key3, int maxpass,
  double border, double maxchisq, double mindeviation,
  long ngiven, int ldxgiven, double* xgiven,
  long nextra, peakfinder_t peakfinder,
  const(char)* statefile, void* spin,
  int *nregions, long *neval, int *fail,
  double* integral, double* error, double* prob);

void llCuhre(int ndim, int ncomp,
  integrand_t integrand, void *userdata, long nvec,
  double epsrel, double epsabs,
  int flags,
  long mineval, long maxeval,
  int key,
  const(char)* statefile, void *spin,
  int *nregions, long *neval, int *fail,
  double* integral, double* error, double* prob);

void cubafork(void *pspin);
void cubawait(void *pspin);

void cubacores(int n, int p);
void cubaaccel(int n, int p);

void cubainit(void function(), void *arg);
void cubaexit(void function(), void *arg);
