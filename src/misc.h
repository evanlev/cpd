
extern int verbose;

void debug_printf(const char *fmt, ... );
void my_assert( int exp, const char *err);
void md_calc_strides(unsigned int D, long str[], const long dim[], size_t size);
long md_calc_size(const unsigned int D, const long *dims);
long sub2ind(const unsigned int D, const long *strides, const long *sub);

void mul(long N, double *dst, const double *src1,  const double *src2);
void mul2(long N, double *dst, const int *src1,  const double *src2);
double sumd(const long N, const double *src);
int sumi(const long N, const int *src);
void ind2sub(const unsigned int D, const long *dims, long *sub, const long ind);
void hardThreshold(const long N, double *x, const double tau);
double normalize( const long N, double *dst, const double *src);
void smpy(long N, double beta, double *dst, const double *src);

void* xmalloc(size_t s);

void rpermute(const long n, long *a);
void randperm( const long n, long perm[]);
int in_bounds( const long D, long *sample, const long *dims );
int arr_is_zero( const unsigned long N, const double *C);
long mod( const long t, const long nt);

