#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "misc.h"
#include "udcpd.h"
#include "cfl/mmio.h"
#include "cfl/io.h"
#include <complex.h>

int main(int argc, char *argv[])
{
    if(argc != 7 ){ /* Check the number of arguments */
        my_assert(0, "Wrong number of input arguments.");
    }

    /* INPUTS */
    const int nt = atoi(argv[1]);
    const double FOVRatio = atof(argv[2]);
    long dims[DIMS];
    /*
    _Complex float *feasiblePoints = load_cfl(argv[3], KDIMS, dims);
    */

    dims[T_DIM] = nt;
    const long shapeOpt = atoi(argv[4]);
    verbose = atoi(argv[5]);
    const double C = atof(argv[6]);
    const char *pattern_cfl_name = argv[7];

    const int isPeriodicInK = 0; /* use for periodic boundary conditions in k-dimension */
    struct pattern_s *data = init_data_str(dims, isPeriodicInK);
    
    double *feasiblePointsD = xmalloc(md_calc_size(KDIMS, dims)*sizeof(double));

    // Run CPD
    genUDCPD(data, feasiblePointsD, FOVRatio, C, shapeOpt);

    
    /*
    // Output as a cfl
    _Complex float *pattern_cfl = create_cfl(pattern_cfl_name, 3, dims);
    for( long i = 0 ; i < data->nptsKT ; i++ ){
        pattern_cfl[i] = (_Complex float) data->masks[i];
    }
    unmap_cfl(DIMS, dims, pattern_cfl);
    */

    free(data->masks);
    free(feasiblePointsD);

}



