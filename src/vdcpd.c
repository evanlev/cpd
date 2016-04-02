
struct pattern_s *init_data_str( const int *dims, const int isPeriodicInK, const int isPeriodicInT)
{
    struct pattern_s *data = (struct pattern_s *) xmalloc(sizeof(struct pattern_s));
    memcpy(data->dims, dims, DIMS*sizeof(int));
    data->numSamples = 0;
    data->nptsKT = md_calc_size(DIMS, data->dims);
    data->nptsK  = md_calc_size(KDIMS, data->dims);
    md_calc_strides(DIMS, data->strides, dims, 1);
    data->masks = xmalloc(data->nptsKT*sizeof(double)); 
    
    data->numSamplest = (long *) xmalloc(dims[T_DIM]*sizeof(long));
    memset(data->numSamplest, 0, dims[T_DIM]*sizeof(long));
    memset(data->masks, 0, data->nptsKT*sizeof(double));
    
    data->isPeriodicInK = isPeriodicInK;
    data->isPeriodicInT = isPeriodicInT;

    return data;
}

struct samplingConstraints*
init_constraints(const long maxSamplesPerPhase, const long maxTotalSamples, 
                const double *feasiblePoints)
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

struct distanceCriteria*
init_dst(const long *dims, const int shapeOpt, const int distRelaxationOpt, const double alpha)
{
    struct distanceCriteria *dist = (struct distanceCriteria *) xmalloc(sizeof(struct distanceCriteria));
    const long N = md_calc_size(DIMS, dims);
    
    /* initialize variables for min distance */
    dist->dky = (double *) xmalloc(N*sizeof(double));
    dist->dt  = (double *) xmalloc(N*sizeof(double));
    dist->alpha = alpha;
    memset(dist->dky, 0, N*sizeof(double));
    memset(dist->dt,  0, N*sizeof(double));
    dist->shapeOpt = shapeOpt;

    /* initialize options for relaxing min distance */
    dist->distRelaxationOpt = distRelaxationOpt;
    dist->dky_scale = SCALE_MINDISTANCE_K;
    dist->dt_scale = SCALE_MINDISTANCE_T;
    dist->dky_thresh = THRESH_MINDISTANCE_K;
    dist->dt_thresh = THRESH_MINDISTANCE_T;
    return dist;
}

/* draw from a 3D PDF conditional */
static void
drawFrom3DPDF( double *pdf, int *sample, const int *dims)
{
    double r = (double)rand()/(double)RAND_MAX; 
    long N = md_calc_size(DIMS, dims);
    double cumsum = 0;
    long i;
    for( i = 0 ; i < N ; i++ ){
        cumsum += pdf[i];
        if( cumsum > r ){
            ind2sub(DIMS, dims, sample, i);
            break;
        }
    }
}

/* ----------------- PDFs -------------------- */
/* generate PDF of the form exp(-lambda*kr) */
static void
genExponentialPDF(double *pdf, const int *dims, const double alphasq, const double krT){
    const long nptsK = md_calc_size(KDIMS, dims);
    const double muy = ((double)dims[Y_DIM])/2;
    const double muz = ((double)dims[Z_DIM])/2;
    int y, z;
    for( z = 0 ; z < dims[Z_DIM] ; z++ )
    for( y = 0 ; y < dims[Y_DIM] ; y++ ){
        double arg = pow((double)y-muy,2) + alphasq*pow((double)z - muz,2);
        double kr = sqrt(arg); 
        pdf[y + dims[Y_DIM]*z] = exp(-1*kr/krT);
    }
    normalize(nptsK, pdf, pdf);
}

/* generate conditional PDF of the form 1/(kr)^order */
static void
genPolynomialPDF(double *pdf, const int *dims, const double alphasq, 
                 const double order){
    const long nptsK = md_calc_size(KDIMS, dims);
    const double muy = ((double)dims[Y_DIM])/2;
    const double muz = ((double)dims[Z_DIM])/2;
    int y, z;
    for( z = 0 ; z < dims[Z_DIM] ; z++ )
    for( y = 0 ; y < dims[Y_DIM] ; y++ ){
        double kr = sqrt( pow((double)y-muy,2) + alphasq*pow((double)z-muz,2));
        pdf[y + dims[Y_DIM]*z] = (double)1/MAX(pow(kr, order), 4);
    }
    normalize(nptsK, pdf, pdf);
}

/* generate gaussian conditional PDF with sigma standard dev in kr */
static void
genGaussianPDF( double *pdf,  const int *dims, 
                const double alphasq, const double *sigmasq)
{
    const long nptsK = md_calc_size(KDIMS, dims);
    const double muy = ((double)dims[Y_DIM])/2;
    const double muz = ((double)dims[Z_DIM])/2;
    int y, z;
    for( z = 0 ; z < dims[Z_DIM] ; z++ )
    for( y = 0 ; y < dims[Y_DIM] ; y++ ){
        double deltaysq = pow(((double)y-muy),2);
        double deltazsq = pow(((double)z-muz),2);
        pdf[y + dims[Y_DIM]*z] = exp(-0.5f*( deltaysq/sigmasq[Y_DIM] + deltazsq/sigmasq[Z_DIM] ));
    }
    normalize(nptsK, pdf, pdf);
}

/* generate uniform conditional PDF */
static void
genUniformPDF( double *pdf, const int *dims){
    const long nptsK = md_calc_size(KDIMS, dims);
    unsigned long i;
    for( i = 0 ; i < nptsK ; i++ )
        pdf[i] = 1.0 / (double)nptsK;
}

/******** 3D PDFs *********/
static void
make3DpdfFrom2Dpdf( const int *dims, double *pdf3d, const double *pdf2d ){
    const long nptsK = md_calc_size(KDIMS, dims);
    int t;
    for( t = 0 ; t < dims[T_DIM] ; t++ ){
        memcpy(&pdf3d[nptsK*t], pdf2d, nptsK*sizeof(double));
    }
    normalize(nptsK*dims[T_DIM], pdf3d, pdf3d);
}

/* generate 3D PDF, polynomial in k, constant in t */
static void
genPolynomial3DPDF( double *pdf3d, const int *dims,
                    const double alphasq, const double order ){
    const long nptsK = md_calc_size(KDIMS, dims);
    double *pdf2d = (double *) xmalloc(nptsK*sizeof(double));
    genPolynomialPDF( pdf2d, dims, alphasq, order );
    make3DpdfFrom2Dpdf( dims, pdf3d, pdf2d);
    free(pdf2d);
}

/* generate 3D PDF, polynomial in k, constant in t */
static void
genUniform3DPDF( double *pdf3d, const int *dims)
{
    const long nptsK = md_calc_size(KDIMS, dims);
    double *pdf2d = (double *) xmalloc(nptsK*sizeof(double));
    genUniformPDF( pdf2d, dims);
    make3DpdfFrom2Dpdf( dims, pdf3d, pdf2d);
    free(pdf2d);
}

/* generate 3D PDF, exponential in k */
static void
genExponential3DPDF( double *pdf3d, const int *dims,
                     const double alphasq, const double kT){
    const long nptsK = md_calc_size(KDIMS, dims);
    double *pdf2d = (double *) xmalloc(nptsK*sizeof(double));
    genExponentialPDF( pdf2d, dims, alphasq, kT);
    make3DpdfFrom2Dpdf( dims, pdf3d, pdf2d);
    free(pdf2d);
}
/* generate 3D PDF, Gaussian in k */
static void
genGaussian3DPDF( double *pdf3d, const int *dims,
                  const double alpha, const double *sigma)
{
    const long nptsK = md_calc_size(KDIMS, dims);
    double *pdf2d = (double *) xmalloc(nptsK*sizeof(double));
    genGaussianPDF( pdf2d, dims, alpha, sigma);
    make3DpdfFrom2Dpdf( dims, pdf3d, pdf2d);
    free(pdf2d);
}


int compareDoubles (const void * a, const void * b)
{
    if ( *(double*)a <  *(double*)b ) return -1;
    if ( *(double*)a == *(double*)b ) return 0;
    if ( *(double*)a >  *(double*)b ) return 1;
}

/* make a map of variable minimum distance in k-space 
 *
 * Pre: pdf is 3D pdf, sum_k,t pdf(k,t) = 1, pdf(k,t) >= 0
 *
 * */
static void
calcMinTDistance( double *dt, const double *pdf, 
                                     const long maxTotalSamples, const int *dims, 
                                     const double dtScaling, const double *feasiblePoints)
{
    double totalR = sumd(md_calc_size(KDIMS,dims), feasiblePoints) / 
                        (double)maxTotalSamples;
    long N = md_calc_size(DIMS, dims);
    long nptsK = md_calc_size(KDIMS, dims);
    long i;
    for( i = 0 ; i < N ; i++ ){
        if( pdf[i] == 0 ){
            /* avoid divide by 0 */
            dt[i] = 0;
        }else{
            /* TODO FIXME check this formula */
             /*
             * dt is chosen for consistency with UD CPD, where p(k) = 1/(ny*nz*nt) and
             * T(k) = # Feasible points / # points in each phase
             *      = 1/(ny*nz*nt * p(k)) * (# Feas. pts. / # pts. in each phase) 
             */
            dt[i] = dtScaling*sumd(nptsK, feasiblePoints)/(double)(maxTotalSamples/dims[T_DIM]*N*pdf[i]);
            /*
             */
        }
        /* Limit minimum distance */
    }
    hardThreshold( N, dt, THRESH_MINDISTANCE_T);
    
    /*
    printf("|F| = %f, max total samples = %d, samples per phase = %ld, dims[T_DIM] = %d", sumd(nptsK, feasiblePoints), maxTotalSamples, maxTotalSamples/(long)dims[T_DIM], dims[T_DIM]);
    printf("\npdf = %f, N = %d, N*pdf = %f\n", pdf[0], N, N*pdf[0]);
    */

    /* Print some statistics */
    double *dt_sort = xmalloc(N*sizeof(double));
    memcpy(dt_sort, dt, N*sizeof(double));
    qsort(dt_sort, N, sizeof(double), compareDoubles);
    mexPrintf("\nMin. temporal dist.: Mean: %3.2e, Median: %3.2e, Min: %3.2e, Max: %3.2e", 
                sumd(N,dt_sort)/(double)N, dt_sort[N/2], dt_sort[0], dt_sort[N-1]);
    free(dt_sort);
    /*
    */
}


/* make a map of variable minimum distance in k-space 
 *
 * Pre: pdf is 3D pdf, sum_k,t pdf(k,t) = 1, pdf(k,t) >= 0
 *
 * */
static void
calcMinKDistance( double *dky, const double *pdf, const float R, 
                          const int *dims, const double betak, 
                          const double *feasiblePoints)
{
    long N = md_calc_size(DIMS, dims);
    long i;
    for( i = 0 ; i < N ; i++ ){
        if( pdf[i] == 0 ){
            /* avoid divide by 0 */
            dky[i] = 0;
        }else{
            /* check htis formula */
            dky[i] = sqrt(R)*betak/sqrt(pdf[i]*N);
            /* */
        }
    }
    hardThreshold( N, dky, THRESH_MINDISTANCE_K);
    /* Print some statistics */
    double dky_sort[N];
    memcpy(dky_sort, dky, N*sizeof(double));
    qsort(dky_sort, N, sizeof(double), compareDoubles);
    mexPrintf("\nMin. k-space dist.:  Mean: %3.2e, Median: %3.2e, Min: %3.2e, Max: %3.2e", 
                sumd(N,dky_sort)/(double)N, dky_sort[N/2], dky_sort[0], dky_sort[N-1]);
    /*
    */
}


/* Multiply an ellipsoid by 0 */
static void
zeroOutKTNeighborsf( double *masks, const int c[],
                  const struct distanceCriteria *dist, const int *dims, 
                  const long *strides,
                  const int isPeriodicInK,
                  const int isPeriodicInT){
   
    int deltay, deltaz, deltat;
    long cInd = sub2ind(DIMS, strides, c);
    masks[cInd] = 0;
    
    double dkysq = pow(dist->dky[cInd],2);
    double dtsq  = pow(dist->dt[cInd],2);
    /* zero out neighbors in k-t space */
    /* loop over a box covering the ellipsoid (max deltat, deltak) and decide to zero */
    int ikt[DIMS];
    for( ikt[Y_DIM] = 0 ; ikt[Y_DIM] < dist->dky[cInd] ; ikt[Y_DIM]++  ){
        for( ikt[Z_DIM] = 0 ; ikt[Z_DIM] < dist->dky[cInd]/dist->alpha ; ikt[Z_DIM]++ ){
            for( ikt[T_DIM] = 0 ; ikt[T_DIM] < dims[T_DIM]  ; ikt[T_DIM]++ ){
                /* scale indices y,z,t */
                deltaz = (double)ikt[Z_DIM] * dist->alpha;
                deltay = (double)ikt[Y_DIM];
                deltat = (double)ikt[T_DIM];
                int isZero = 0;
                switch( dist->shapeOpt ){
                    case CROSS:
                    {
                        /* Cross shape */
                        double normDistl2k = (pow(deltaz,2) + pow(deltay,2))/dkysq;
                        if( normDistl2k == 0.0 || ( ikt[T_DIM] == 0 && normDistl2k < 1) ){
                            isZero = 1;
                        }
                        break;
                    }
                    case L1_BALL:
                    {
                        /* l1 ball */
                        double normDistl1 = (deltaz + deltay)/dist->dky[cInd] + deltat/dist->dt[cInd];
                        if( normDistl1 <= 1 ){
                            isZero = 1;
                        }
                        break;
                    }
                    case ELLIPSOID:
                    {
                        /* l2 ball */
                        double normDistl2k = (pow(deltaz,2) + pow(deltay,2))/dkysq;
                        double normDistl2  = normDistl2k +pow(deltat,2)/dtsq;
                        if( normDistl2 <= 1 ){
                            isZero = 1;
                        }
                        break;
                    }
                    case CONES:
                    {
                        /* Cones shape - ellipse in ky,kz, line in t with line through ky,kz = 0 */
                        double normDistl2k = (pow(deltaz,2) + pow(deltay,2))/dkysq;
                        if( normDistl2k <= 1 - fabs(ikt[T_DIM])/dist->dt[cInd] ){
                            isZero = 1;
                        }
                        break;
                    }
                    case PLANE_AND_CONES:
                    {
                        /* Union of plane at t= 0 and cone at y,z,t = 0 for radial */
                        double normDistl2k = (pow(deltaz,2) + pow(deltay,2))/dkysq;
                        if( ikt[T_DIM] == 0 || normDistl2k <= 1 - fabs(ikt[T_DIM])/dist->dt[cInd] ){
                            isZero = 1;
                        }
                        break;
                    }

                }
                /* Shape has symmetry such that all octants are the same */
                if( isZero ){
                    int cNew[3];
                    int tsgn, ysgn, zsgn;
                    for( tsgn = -1 ; tsgn <= 1 ; tsgn += 2 )
                    for( ysgn = -1 ; ysgn <= 1 ; ysgn += 2 )
                    for( zsgn = -1 ; zsgn <= 1 ; zsgn += 2 ){
                        cNew[Y_DIM] = c[Y_DIM] + ysgn*ikt[Y_DIM];
                        cNew[Z_DIM] = c[Z_DIM] + zsgn*ikt[Z_DIM];
                        cNew[T_DIM] = c[T_DIM] + tsgn*ikt[T_DIM];
                        
                        if( isPeriodicInK ){
                            cNew[Y_DIM] = mod( cNew[Y_DIM], dims[Y_DIM]);
                            cNew[Z_DIM] = mod( cNew[Z_DIM], dims[Z_DIM]);
                        }
                        if( isPeriodicInT ){
                            cNew[T_DIM] = mod( cNew[T_DIM], dims[T_DIM]);
                        }
                        if( in_bounds(DIMS, cNew, dims) ){
                            masks[sub2ind(DIMS, strides, cNew)] = 0;
                        }
                    }
                }
            } /* end for t = 0 ; t < */
        } /* end for z = 0 ; */
    } /* end for y = 0 ; y < ... */
}

/* return cpdf = p(ky, kz, t | feasiblePoints, dt, dk, nSamples(t) ) 
 * 
 * pdf = PDF p(ky, kz, t) 
 * M   = sampling mask
 * feasilbePoints  = 
 *
 * If it is overconstrained, returns 0 and cpdf is all zero, otherwise returns 1
 *
 * */
static int
calcCPDF(double *cpdf, const double *pdf, const struct pattern_s *data,
         struct samplingConstraints *constraints){
         
    const long nptsK = data->nptsK;
    const long N = data->nptsKT;
   
    /* max total samples */
    if( data->numSamples >= constraints->maxTotalSamples ){
        memset(cpdf, 0, N*sizeof(double));
        return 0;
    }
    
    /* cpdf = pdf */
    memcpy( cpdf, pdf, N*sizeof(double) );
    
    /* Condition on feasible points */
    int t;
    for( t = 0 ; t < data->dims[T_DIM] ; t++ ){
        mul(nptsK, &cpdf[t*nptsK], constraints->feasiblePoints, &pdf[t*nptsK]);
    }
    
    /* max samples per phase */
    for( t = 0 ; t < data->dims[T_DIM] ; t++ ){
        if( data->numSamplest[t] >= constraints->maxSamplesPerPhase ){
            memset(&cpdf[t*nptsK], 0, sizeof(double)*nptsK);
            /*
            writeDataRealDoubles("cpdfmx", "cpdfmx.dims", cpdf, DIMS, data->dims);
            */
        }
    }
    
    /* min distance constraints */
    long i;
    for( i = 0 ; i < N ; i++ ){
        if( data->masks[i] != 0.0 ){
            int s[DIMS];
            ind2sub(DIMS, data->dims, s, i);
            /* M(ky, kz, t) constraint: can't sample this point or its neighbors 
             * zero out cross or ellipsoid centered here 
             */
            zeroOutKTNeighborsf( cpdf, s, constraints->dist, data->dims, data->strides, data->isPeriodicInK, data->isPeriodicInT); 
        }
    }
    /* normalize cpdf */
    return (normalize(N, cpdf, cpdf) > 0);
}

/* distRelaxationOpt = 1 for k-space and t-distance
 * return 1/0 = success/failure, failure if distance constraints are fully relaxed
 */
static int
relaxDistanceConstraints( const long N, struct distanceCriteria *dist )
{
    smpy( N, dist->dky_scale, dist->dky, dist->dky);
    hardThreshold( N, dist->dky, dist->dky_thresh);
    switch( dist->distRelaxationOpt ){
        case DIST_KT:
        {
            smpy( N, dist->dt_scale, dist->dt, dist->dt);
            hardThreshold( N, dist->dt,  dist->dt_thresh);
            if( arrIsZero(N, dist->dky) && arrIsZero(N, dist->dt) ){
                return 0;
            }else{
                return 1;
            }
            break;
        }
        case DIST_K:
        {
            if( arrIsZero(N, dist->dky) ){
                return 0;
            }else{
                return 1;
            }
            break;
        }
        default:
        {
            mexErrMsgTxt("Unrecognized option.");
            break;
        }
    }
}

/* main algorithm- 3D = 2D + time */
static void
genVDCPD( struct pattern_s *data, const double *pdf, 
          struct samplingConstraints *constraints,
          const int periodizeSamples){
        
    /* Print the relative change in min distance */
    double dkyInit = constraints->dist->dky[0]; 
    double dtInit = constraints->dist->dt[0]; 
    
    /* Make a feasbile cpdf from the PDF */
    double *cpdf = (double *) xmalloc(data->nptsKT*sizeof(double));
    if( 0 == calcCPDF(cpdf, pdf, data, constraints) ){
        mexErrMsgTxt("Failed to make conditional PDF. Must be infeasible, so check the parameters.");
    }

    /* write initial data */
    mexPrintf("\nBegin randomized sampling...");
    int stoppingCriteria = 0;
    long numSamplesT = 0;
    int sample[DIMS];
    while( !stoppingCriteria ){
        
        /* Draw (sy, sz, st) ~ cpdf */    
        drawFrom3DPDF( cpdf, sample, data->dims);
        
        /* Add (sy, sz, st) */
        long sampleInd = sub2ind(DIMS, data->strides, sample);
        data->masks[sampleInd] = 1;
        data->numSamples++;
        data->numSamplest[sample[T_DIM]]++;
       
       if( periodizeSamples ){
            /* Add samples {(sy, sz, st + n T(sy,sz))}_n, 
             * T is chosen to be correct in UD CPD, where p(k) = 1/(ny*nz*nt) and
             * T(k) = R
             *      = 1/(ny*nz*nt * p(k)) * R 
             *      = 1/(ny*nz*nt * p(k)) * (ny*nz*nt / maxpts),
             *      = 1/(p(k)*maxpts)
             * where R is the reduction factor */
            double deltaTInc = sumd(data->nptsK, constraints->feasiblePoints)
                                /(double)(constraints->maxTotalSamples/data->dims[T_DIM]*data->nptsKT*
                                 pdf[sampleInd]);
            double deltat;
            for( deltat = deltaTInc ; deltat < data->dims[T_DIM] ; deltat += deltaTInc)
            {
                int tsgn;
                for( tsgn = -1 ; tsgn <= 1 ; tsgn += 2 ){ 
                    int sPlusT[3] = {sample[Y_DIM], sample[Z_DIM], 
                                        sample[T_DIM] + (int) tsgn*roundf(deltat)};
                    long sPlusTInd = sub2ind(DIMS, data->strides, sPlusT);
                    if( in_bounds(DIMS, sPlusT, data->dims) && data->numSamples < constraints->maxTotalSamples && !data->masks[sPlusTInd] && data->numSamplest[sPlusT[T_DIM]] < constraints->maxSamplesPerPhase ){
                        data->masks[sPlusTInd] = 1;
                        data->numSamples++;
                        data->numSamplest[sPlusT[T_DIM]]++;
                    }
                }
            } /* end deltat loop */
        }        

        /* recompute the conditional PDF */
        while( !stoppingCriteria &&
               0 == calcCPDF(cpdf, pdf, data, constraints) ){
            /* relax distance constraints */ 
            if( 0 == relaxDistanceConstraints(data->nptsKT, constraints->dist ) ){
                /* infeasible */
                mexPrintf(" (stopping)");
                stoppingCriteria = 1;
            }
            /*
            printf("scaling dky, sum = %f...\n", sumd(data->nptsKT, constraints->dist->dky));
            */
        }

        /* Determine if it needs to break */
        stoppingCriteria = stoppingCriteria || (numSamplesT == data->numSamples) || (data->numSamples >= constraints->maxTotalSamples);
        if( data->numSamples % (constraints->maxTotalSamples/10) == 0 ){
            
            mexPrintf("\nProgress: %3d%%, ky min. dist.: %.2f%%, kt min. dist.: %.2f%%", 
                        (100*data->numSamples + constraints->maxTotalSamples - 1)/constraints->maxTotalSamples,
                        100*constraints->dist->dky[0]/dkyInit,
                        100*constraints->dist->dt[0]/dtInit );
            mexEvalString("fprintf('')");
            /*
            mexEvalString("disp(sprintf('\r%d ', (100*numSamples)/maxTotalSamples))");
            printf("numSamples = %d\n", data->numSamples);
            */
        }
        numSamplesT = data->numSamples;
        
    } /* end While */
    
    /* clean up */
    free(cpdf);
    mexPrintf("\n");
}


