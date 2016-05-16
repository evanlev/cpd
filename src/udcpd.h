#ifndef UDCPD_H
#define UDCPD_H 1

#ifdef __cplusplus
extern "C" {
#endif



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


/* min distance shape options */
enum shape_opt {CROSS, L1_BALL, L2_BALL, CONES, PLANE_AND_CONES};

/* min distance relaxation options */
#define DIST_K 0
#define DIST_KT 1

struct pattern_s{
    int *masks;
    int *maskTot;
    long dims[DIMS];
    long pat_strs[DIMS];
    long nptsKT;
    long nptsK;
    long *numSamplest;
    long numSamples;
    int isPeriodicInK;
};

struct distanceCriteria{
    /* min distance variables */
    double dky_min;
    double dt_min;
    enum shape_opt shapeOpt;

    double alpha;

    /* options for relaxing min distance */
    double dky_min_scale;
    double dt_min_scale;
    double dky_min_thresh;
    double dt_min_thresh;
    
};

struct samplingConstraints{
    /* min distance */
    struct distanceCriteria *dist;

    /* max samples */
    const long *maxSamplesPerPhase;
    long maxTotalSamples;

    /* defines k-space subdomain */
    const int *feasiblePoints;
};

extern struct pattern_s *
init_data_str( const long *dims , const int isPeriodicInK);

void genUDCPD(const long *dims, 
              int *pattern,
                    const int *feasiblePoints, 
                    const double FOVRatio, 
                    const double C, 
                    const long shapeOpt,
                    const double ky_mindist);

#ifdef __cplusplus
}
#endif

#endif // UDCPD_H

