
#define KDIMS 2
#define DIMS 3
#define Y_DIM 0u
#define Z_DIM 1u
#define T_DIM 2u
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) >= (Y) ? (X) : (Y))
#define SCALE_MINDISTANCE_K 0.9
#define SCALE_MINDISTANCE_T 0.9
#define THRESH_MINDISTANCE_K 0.25
#define THRESH_MINDISTANCE_T 0.1

/* vd pdf options - polynomial is 3+ */ 
enum pdfs {GAUSSIAN, EXPONENTIAL, UNIFORM, POLYNOMIAL};

/* min distance shape options */
enum shape_opt {CROSS, L1_BALL, ELLIPSOID, CONES, PLANES_AND_CONES};

/* min distance relaxation options */
#define DIST_K 0
#define DIST_KT 1



struct pattern_s{
    double *masks;
    long dims[DIMS];
    long strides[DIMS];
    long nptsKT;
    long nptsK;
    long *numSamplest;
    long numSamples;
    int isPeriodicInK;
    int isPeriodicInT;
};

struct distanceCriteria{
    /* min distance variables */
    double *dky;
    double *dt;
    int shapeOpt;

    double alpha;

    /* options for relaxing min distance */
    int distRelaxationOpt; /* 1 = k-space and time, 0 = k-space only */
    double dky_scale;
    double dt_scale;
    double dky_thresh;
    double dt_thresh;
    
};

struct samplingConstraints{
    /* min distance */
    struct distanceCriteria *dist;

    /* max samples */
    long maxSamplesPerPhase;
    long maxTotalSamples;

    /* defines k-space subdomain */
    const double *feasiblePoints;
};

struct pattern_s *init_data_str( const int *dims, const int isPeriodicInK, const int isPeriodicInT)
struct samplingConstraints*
init_constraints(const long maxSamplesPerPhase, const long maxTotalSamples, 
                const double *feasiblePoints)

struct distanceCriteria*
init_dst(const long *dims, const int shapeOpt, const int distRelaxationOpt, const double alpha)


