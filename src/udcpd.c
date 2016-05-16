#include <string.h> /* memset */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef MEX_COMPILE_FLAG
#include <mex.h>
#endif
#include "misc.h"
#include "udcpd.h"

#ifdef __cplusplus
extern "C" {
#endif

struct pattern_s *
init_data_str( const long *dims , const int isPeriodicInK)
{
    struct pattern_s *data = (struct pattern_s *) xmalloc(sizeof(struct pattern_s));
    
    memcpy(data->dims, dims, DIMS*sizeof(long));
    data->numSamples = 0;
    data->nptsKT = md_calc_size(DIMS, data->dims);
    data->nptsK  = md_calc_size(KDIMS, data->dims);
    md_calc_strides(DIMS, data->pat_strs, dims, 1);
    data->masks = xmalloc(data->nptsKT*sizeof(int)); 
    data->maskTot = xmalloc(data->nptsK*sizeof(int)); 
    
    data->numSamplest = (long *) xmalloc(dims[T_DIM]*sizeof(long));
    memset(data->numSamplest, 0, dims[T_DIM]*sizeof(long));
    memset(data->masks, 0, data->nptsKT*sizeof(int));
    
    data->isPeriodicInK = isPeriodicInK;
    return data;
}

static struct samplingConstraints*
init_constraints(const long *maxSamplesPerPhase, const long maxTotalSamples, 
                const int *feasiblePoints)
{
    struct samplingConstraints *constraints = (struct samplingConstraints *) xmalloc(sizeof(struct samplingConstraints));
    /*
    constraints->dist = NULL;
    */
    constraints->maxSamplesPerPhase = maxSamplesPerPhase;
    constraints->maxTotalSamples = maxTotalSamples;
    constraints->feasiblePoints = feasiblePoints;
    return constraints;
}

static struct distanceCriteria*
init_dst(const int shapeOpt, const double FOVRatio)
{
    struct distanceCriteria *dist = (struct distanceCriteria *) xmalloc(sizeof(struct distanceCriteria));
    
    /* initialize variables for min distance */
    dist->dky_min = 0;
    dist->dt_min  = 0;
    dist->alpha = FOVRatio;
    dist->shapeOpt = shapeOpt;

    /* initialize options for relaxing min distance */
    dist->dky_min_scale  = SCALE_MINDISTANCE_K;
    dist->dt_min_scale   = SCALE_MINDISTANCE_T;
    dist->dky_min_thresh = THRESH_MINDISTANCE_K;
    dist->dt_min_thresh  = THRESH_MINDISTANCE_T;
    return dist;
}

/* Multiply an ellipsoid by 0 */
static void
zeroOutKTNeighborsf( double *cpdf, const long sample[],
                  const struct distanceCriteria *dist, const long *dims, 
                  const long *pat_strs,
                  const int isPeriodicInK){
   
    long sampleInd = sub2ind(DIMS, pat_strs, sample);
    cpdf[sampleInd] = 0;
    
    double dky_minsq = pow(dist->dky_min,2);
    /* zero out neighbors in k-t space */
    /* loop over a box covering the ellipsoid (max deltat, deltak) and decide to zero */
    int ikt[DIMS];
    for( ikt[Y_DIM] = 0 ; ikt[Y_DIM] < dist->dky_min ; ikt[Y_DIM]++  ){
        for( ikt[Z_DIM] = 0 ; ikt[Z_DIM] < dist->dky_min/dist->alpha ; ikt[Z_DIM]++ ){
            for( ikt[T_DIM] = 0 ; ikt[T_DIM] < dims[T_DIM]  ; ikt[T_DIM]++ ){
                double deltaz = (double)ikt[Z_DIM] * dist->alpha;
                double deltay = (double)ikt[Y_DIM];
                double deltat = (double)ikt[T_DIM];
                int isZero = 0;
                if( ikt[Z_DIM] == 0 && ikt[Y_DIM] == 0 ){
                    isZero = 1;
                }else{
                    switch( dist->shapeOpt ){
                        case CROSS:
                        {
                            /* Cross shape */
                            double normDistl2k = (pow(deltaz,2) + pow(deltay,2)) / 
                                                    dky_minsq;
                            if( (normDistl2k == 0.0) || ( ikt[T_DIM] == 0 && normDistl2k < 1) ){
                                isZero = 1;
                            }
                            break;
                        }
                        case L1_BALL:
                        {
                            /* l1 ball */
                            double normDistl1 = (deltaz + deltay)/dist->dky_min + deltat/dist->dt_min;
                            if( normDistl1 <= 1 ){
                                isZero = 1;

                            }
                        }
                        case L2_BALL:
                        {                            
                            double normDistl2t = deltat*deltat / pow(dist->dt_min,2);
                            double normDistl2k = (pow(deltaz,2) + pow(deltay,2)) / 
                                                    dky_minsq;
                            double normDistl2  = normDistl2k + normDistl2t;
                            if( normDistl2 <= 1 ){
                                isZero = 1;
                            }
                            break;
                        }
                        case CONES:
                        {
                            /* Cones shape - ellipse in ky,kz, line in t with line through ky,kz = 0 */
                            double normDistl2k = pow(deltaz,2) + pow(deltay,2);
                            if( normDistl2k + dist->dky_min/dist->dt_min*fabs(ikt[T_DIM]) < dky_minsq ){
                                isZero = 1;
                            }
                            break;
                        }
                        case PLANE_AND_CONES:
                        {
                            /* Union of plane at t= 0 and cone at y,z,t = 0 for radial */
                            double normDistl2k = (pow(deltaz,2) + pow(deltay,2)) / 
                                                    dky_minsq;
                            if( ikt[T_DIM] == 0 || normDistl2k <= 1 - fabs(ikt[T_DIM])/dist->dt_min ){
                                isZero = 1;
                            }
                            break;
                        }
                    }
                }
                /* Shape has symmetry such that all octants are the same */
                if( isZero ){
                    long neighborSample[3];
                    int tsgn, ysgn, zsgn;
                    for( tsgn = -1 ; tsgn <= 1 ; tsgn += 2 )
                    for( ysgn = -1 ; ysgn <= 1 ; ysgn += 2 )
                    for( zsgn = -1 ; zsgn <= 1 ; zsgn += 2 ){
                        if( isPeriodicInK ){
                            neighborSample[Y_DIM] = mod( sample[Y_DIM] + ysgn*ikt[Y_DIM], dims[Y_DIM] );
                            neighborSample[Z_DIM] = mod( sample[Z_DIM] + zsgn*ikt[Z_DIM], dims[Z_DIM] );
                            neighborSample[T_DIM] = mod( sample[T_DIM] + tsgn*ikt[T_DIM], dims[T_DIM]);
                        }else{
                            neighborSample[Y_DIM] = sample[Y_DIM] + ysgn*ikt[Y_DIM];
                            neighborSample[Z_DIM] = sample[Z_DIM] + zsgn*ikt[Z_DIM];
                            neighborSample[T_DIM] = sample[T_DIM] + tsgn*ikt[T_DIM];
                        }
                        if( in_bounds(DIMS, neighborSample, dims) ){
                            cpdf[sub2ind(DIMS, pat_strs, neighborSample)] = 0;
                        }
                    }
                }
            } /* end for t = 0 ; t < */
        } /* end for z = 0 ; */
    } /* end for y = 0 ; y < ... */
}

void
addSamplesAtMinimumDistance( struct pattern_s *data, 
                             struct samplingConstraints *constraints){
    
    long sample[DIMS];
    long N = data->nptsKT;
    long nptsK = data->nptsK;

    /* Initialize conditional PDF. Note: does not sum to 1 */
    double* cpdf = (double *) xmalloc(sizeof(double)*N);
    for( long i = 0 ; i < N ; i++ )
        cpdf[i] = 1;
    
    /* Condition on feasible points */
    for( long t = 0 ; t < data->dims[T_DIM] ; t++ ){
        mul2(nptsK, &cpdf[t*nptsK], constraints->feasiblePoints, &cpdf[t*nptsK]);
    }
    
    /* Condition on points satisfying minimum "distance" */
    for( long i = 0 ; i < N ; i++ ){
        if( data->masks[i] ){
            ind2sub(DIMS, data->dims, sample, i);
            zeroOutKTNeighborsf( cpdf, sample, constraints->dist, data->dims, 
                                  data->pat_strs, data->isPeriodicInK);
        }
    }

    /* Set up random queue */
    long *randQ = xmalloc(N*sizeof(long));
    long lenRandQ = 0; /* length of the random queue */
    for( long i = 0 ; i < N ; i++ ){
        if( cpdf[i] > 0 ){
            randQ[lenRandQ++] = i;
        }
    }
    if( lenRandQ == 0 ){
        debug_printf("rand Q empty, breaking...\n");
        free(cpdf);
        return;
    }
    randperm(lenRandQ, randQ);

    for( long i = 0 ; i < lenRandQ && data->numSamples < constraints->maxTotalSamples ; i++ ){
        long sampleInd = randQ[i];
        long sample[DIMS];
        if( sampleInd < 0 || sampleInd > N ){
            my_assert(0, "sample out of range.");
        }
        ind2sub(DIMS, data->dims, sample, sampleInd);
        if( data->numSamplest[sample[T_DIM]] < constraints->maxSamplesPerPhase[sample[T_DIM]] 
                && data->numSamples < constraints->maxTotalSamples 
                && cpdf[sampleInd] > 0 ){
            
            /* add randQ[i] to current mask */
            data->masks[sampleInd] = 1;
            data->maskTot[sampleInd % data->nptsK] = 1;
            data->numSamplest[sample[T_DIM]]++;
            data->numSamples++;
            
            zeroOutKTNeighborsf(cpdf, sample, constraints->dist, data->dims, 
                                data->pat_strs, data->isPeriodicInK);
       }
    } /* end randQ loop */
    free(randQ);   
    free(cpdf);
}


/* 
 * return 1/0 = success/failure, failure if distance constraints are fully relaxed
 */
static int
relax_distance_constraints( struct distanceCriteria *dist)
{
    
    if( dist->dt_min == 0 && dist->dky_min == 0 ){
        return 0;
    }
    dist->dt_min  *= dist->dt_min_scale;
    if( dist->dt_min < dist->dt_min_thresh ){
        dist->dt_min = 0;
    }
    
    /* relax k-space distance */
    dist->dky_min *= dist->dky_min_scale;
    if( dist->dky_min < dist->dky_min_thresh ){
        dist->dky_min = 0;
    }
    return 1;
}

void genUDCPD(const long *dims, 
              int *pattern,
                    const int *feasiblePoints, 
              const double FOVRatio, 
              const double C, 
              const long shapeOpt,
              const double mindist_scale){
    
    struct pattern_s *data = init_data_str(dims, 0);

    /* max samples per phase = ny*nz / # phases */
    long maxTotalSamples = (long) sumi(data->nptsK, feasiblePoints);
    long maxSamplesPerPhase[dims[T_DIM]];
    for( long i = 0 ; i < dims[T_DIM] ; i++ ){
        maxSamplesPerPhase[i] = maxTotalSamples/dims[T_DIM];
    }
    /* + remainder */
    if( maxTotalSamples % dims[T_DIM]){
        for( long i = 0 ; i <  maxTotalSamples % dims[T_DIM] ; i++ ){
            maxSamplesPerPhase[i]++;
        }
    }

    /* Min "distance" options */
     
    /* Init min distance constraints */
    debug_printf("------------- Uniform density CPD -------------\n");
    debug_printf("FOVz/FOVy = %.2f, ", FOVRatio);
    debug_printf("Dims: %d %d %d", dims[Y_DIM], dims[Z_DIM], dims[T_DIM]);
    debug_printf("\nMin. distance criterion: ");
    switch( shapeOpt ){
        case CROSS:
        {
            debug_printf("cross");
            break;
        }
        case L2_BALL:
        {
            debug_printf("ellipsoid");
            break;
        }
        case L1_BALL:
        {
            debug_printf("l1 ball");
            break;
        }
        case CONES:
        {
            debug_printf("cones");
            break;
        }
        default:
        {
            my_assert(0, "Unrecognized min. distance criterion");
        }
    }
    debug_printf("\nMax total samples = %d", maxTotalSamples);
    debug_printf("\nMax # samples per phase: ");

    for( long i = 0 ; i < dims[T_DIM] ; i++ ){
        debug_printf("%d ", maxSamplesPerPhase[i]);
    }
    debug_printf("\n-----------------------------------------------\n");
    struct samplingConstraints *constraints = init_constraints(maxSamplesPerPhase, maxTotalSamples, feasiblePoints);
    constraints->dist = init_dst(shapeOpt, FOVRatio);
    constraints->dist->dky_min = 100*(mindist_scale+0.001)*sqrt(MIN( (0.5f*FOVRatio) / maxSamplesPerPhase[0], MIN(9*pow(FOVRatio,3), 9.0f ) ));
    constraints->dist->dt_min  = C*constraints->dist->dky_min;


    int iter = 0;
    while( data->numSamples < constraints->maxTotalSamples ){
        addSamplesAtMinimumDistance(data, constraints);
        iter++;
        debug_printf("It: %2d, # Samples: %7d, ky min. dist.: %.1f\n", iter, 
                     data->numSamples, constraints->dist->dky_min); 

        if( data->numSamples < constraints->maxTotalSamples && !relax_distance_constraints( constraints->dist ) ){
            debug_printf("WARNING: could not fill sampling pattern...\n");
            break;
        }
    }
    
    for( long i = 0 ; i < data->nptsKT ; i++ ){
        pattern[i] = data->masks[i];
    }
    
    free(data->masks);
}

#ifdef __cplusplus
}
#endif

