
#include <stdio.h>
#include <stdarg.h>
#include <string.h> /* memset */
#include "misc.h"
#include <stdlib.h>
#include <math.h>

#ifdef MEX_COMPILE_FLAG
#include <mex.h>
#endif


#ifdef __cplusplus
extern "C" {
#endif

int verbose = 0;

void debug_printf(const char *fmt, ... )
{
    if( verbose ){
        va_list ap;
        va_start(ap, fmt);

        char buff[512] = "";
        
        vsprintf(buff, fmt, ap);

        int openok = 1;
#ifdef MEX_COMPILE_FLAG
        mexPrintf(buff);
#else
        fprintf(stdout, buff);
#endif
        va_end(ap);
    }
}


void my_assert( int exp, const char *err){
    if( !exp ){
#ifdef MEX_COMPILE_FLAG
        mexErrMsgTxt(err);
#else
        fprintf(stderr, err);
        exit(0);
#endif
    }
}

void md_calc_strides(unsigned int D, long str[], const long dim[], size_t size)
{
	long old = size;

    unsigned int i;
	for (i = 0; i < D; i++) {
		str[i] = (1 == dim[i]) ? 0 : old;
		old *= dim[i];
	}
}


/* product of D-length array dims[] */
long
md_calc_size(const unsigned int D, const long *dims)
{
    if( D <= 0 ){
        my_assert(0, "D <= 0");
    }
    return (D == 1) ? dims[0] : dims[D-1]*md_calc_size(D-1, dims);
}

/* convert subscript to column major linear index */
long
sub2ind(const unsigned int D, const long *strides, const long *sub)
{
    long res = 0;
    unsigned int i;
    for( i = 0 ; i < D ; i++ ){
        res += sub[i]*strides[i];
    }
    return res;
}

/* mul2 */
void mul2(long N, double *dst, const int *src1,  const double *src2)
{
    long i;
    for( i =0 ; i < N ; i++ ){
        dst[i] = src1[i]*src2[i];
    }
}

/* mul */
void mul(long N, double *dst, const double *src1,  const double *src2)
{
    long i;
    for( i =0 ; i < N ; i++ ){
        dst[i] = src1[i]*src2[i];
    }
}

/* divide an array by denom */
void
smpy(long N, double beta, double *dst, const double *src)
{
    long i;
    for( i = 0 ; i < N ; i++ ){
        dst[i] = beta*src[i];
    }
}

/* dst[i] = src[i]/sumd(src) */
double normalize( const long N, double *dst, const double *src){
    double tot = sumd(N, src);
    if( tot > 0 ){
        smpy(N, 1.0/tot, dst, src);
        return tot;
    }else{
        return 0;
    }

}

/* hard thresholding */
void hardThreshold(const long N, double *x, const double tau){
    int i;
    for( i = 0 ; i < N ; i++ ){
        if( fabs(x[i]) < tau ){
            x[i] = 0;
        }
    }
}



/* sum */
int sumi(const long N, const int *src)
{
    int res = 0.0f;
    long i;
    for( i = 0 ; i < N ; i++ ){
        res += src[i];
    }
    return res;
}

/* sum */
double sumd(const long N, const double *src)
{
    double res = 0.0f;
    long i;
    for( i = 0 ; i < N ; i++ ){
        res += src[i];
    }
    return res;
}


/* convert linear index to MD subscript */
void
ind2sub(const unsigned int D, const long *dims, long *sub, const long ind)
{
    
    unsigned int i;
    if( D == 1 ){
        sub[0] = ind;
    }else if( D == 2){
        sub[1] = ind / dims[0];
        sub[0] = ind - sub[1]*dims[0];
    }else if( D == 3){
        sub[2] = ind / (dims[0]*dims[1]);
        sub[1] = (ind - sub[2]*dims[0]*dims[1]) / dims[0];
        sub[0] = (ind - sub[2]*dims[0]*dims[1] - sub[1]*dims[0]);
    }
}


void* xmalloc(size_t s)
{
	void* p = malloc(s);

	if (NULL == p)
        my_assert(0, "Could not allocate memory");

	return p;
}

void rpermute(const long n, long *a) {
    long k;
    for (k = 0; k < n; k++)
        a[k] = k;
    for (k = n-1; k > 0; k--) {
        long j = rand() % (k+1);
        long temp = a[j];
        a[j] = a[k];
        a[k] = temp;
    }
}

/* randomly permutes the elements of the array perm of length n */
void randperm( const long n, long perm[])
{
    long *ind = xmalloc(n*sizeof(long));
    long *tmp = xmalloc(n*sizeof(long));
    long k;
    rpermute(n, ind);
    for( k = 0 ; k < n ; k++ ){
        tmp[k] = perm[ind[k]];
    }
    for( k = 0 ; k < n ; k++ ){
        perm[k] = tmp[k];
    }
    free(ind);
    free(tmp);
}
                    
/* returns if samples is in bounds */
int in_bounds( const long D, long *sample, const long *dims ){
    if( D == 1 ){
        return (sample[0] >= 0 && sample[0] < dims[0]);
    }else{
        return (sample[0] >= 0 && sample[0] < dims[0])
                  && in_bounds(D-1, &sample[1], &dims[1]);
    }
}

/* 1 if C[0 .. N-1] is all zeros */
int arr_is_zero( const unsigned long N, const double *C)
{
    unsigned long i;
    int res = 1;
    for( i = 0 ; i < N ; i++ ){
        if( C[i] != 0.0 ){
            res = 0;
            break;
        }
    }
    return res;
}

/* mod function */
long mod(const long x, const long y){
    long r = x % y;
    return r < 0 ? r + y : r;
}


#ifdef __cplusplus
}
#endif

