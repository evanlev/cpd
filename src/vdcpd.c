#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>


#include "misc.h"
#include "udcpd.h"
#include "vdcpd.h"

#ifdef __cplusplus
extern "C" {
#endif


void genVDCPD(const long dims[], int *pattern, 
              const int *feasiblePointsD, const float FOVRatio, 
              const float C, const int shapeOpt, const float mindist_scaling, 
              const float vd_exp, const int maxR){
    
    /* Minimum size of a region */
    const int min_region_dim = 2; 
    const int corner_cut = 1;
    
    my_assert(FOVRatio > 0, "FOVRatio <= 0");
    my_assert(maxR > 0, "maxR <= 0");
    my_assert(vd_exp > 0, "vd_exp <= 0");
    int skip_reg = (int) ceilf((double) min_region_dim / ((double) MIN(dims[0], dims[1]) / (double) maxR));
    debug_printf("---------- Variable Density CPD -------- \n");
    debug_printf("Rmax:\t%d\n", maxR);
    debug_printf("FOVRatio:\t%f\n", FOVRatio);
    debug_printf("Density falloff:\t1/%f\n", vd_exp);
    debug_printf("C:\t%f\n", C);
    debug_printf("Shape option:\t%d\n", shapeOpt);
    debug_printf("Region skip:\t%d\n", skip_reg);
    debug_printf("---------------------------------------- \n");
    
    long f_size = dims[0]*dims[1];
    
    /* Clear out the pattern */
    memset(pattern, 0, sizeof(int)*f_size*dims[T_DIM]);
     
    double alpha_sq = pow( (double) dims[Y_DIM] / (double) dims[Z_DIM], 2);
    double * kr = (double *) xmalloc(dims[Y_DIM]*dims[Z_DIM]*sizeof(double));
    int *region     = (int *) xmalloc(dims[Y_DIM]*dims[Z_DIM]*sizeof(int));
    debug_printf("Region assignments...\n");

    for( int iky = 0; iky < dims[Y_DIM] ; iky++ )
    for( int ikz = 0; ikz < dims[Z_DIM] ; ikz++ ){
        long ik = iky + ikz * dims[Y_DIM];
        kr[ik] = sqrt( pow( (double) iky - (double) dims[Y_DIM] / 2.0,2) + alpha_sq*pow((double) ikz - (double) dims[Z_DIM] / 2.0,2));
        if( !corner_cut || kr[ik] <= (double) dims[Y_DIM]/2 ){
            region[ik] = MIN(1 + (int) (skip_reg * floorf( (double) maxR / (double) skip_reg * pow(kr[ik]/(dims[Y_DIM]/2),vd_exp))), maxR);
        }else{
            // Do not sample
            region[ik] = -1;
        }
    }
    int *feasiblePointsReg = xmalloc(f_size*sizeof(int));
    debug_printf("Region-wise CPD\n");
    for( int reg = 1 ; reg <= maxR ; reg++ ){
        debug_printf("%d/%d\n", reg, maxR);
        int *masks_reg = xmalloc(f_size*reg*sizeof(int));
        memcpy(feasiblePointsReg, feasiblePointsD, f_size*sizeof(int));
        /*
        */
        for( long ik = 0 ; ik < f_size ; ik++ ){
            feasiblePointsReg[ik] = (region[ik] == reg);
        }
        long dims_tmp[DIMS];
        for( unsigned int i = 0 ;i < DIMS ; i++ )
            dims_tmp[i] = dims[i];
        
        dims_tmp[T_DIM] = reg;
        
        genUDCPD(dims_tmp, masks_reg, feasiblePointsReg, FOVRatio, C, shapeOpt, mindist_scaling);
        for( long t = 0 ; t < dims[T_DIM] ; t++ )
        for( long ik = 0 ; ik < dims[Y_DIM]*dims[Z_DIM] ; ik++ ){
            if( masks_reg[ik + dims[Y_DIM]*dims[Z_DIM]*(t % reg)] )
                pattern[ik + dims[Y_DIM]*dims[Z_DIM]*t] = 1;
        }
        /*
        */
        free(masks_reg);
    }
    debug_printf("Cleanup...\n");
    free(feasiblePointsReg);
    free(kr); 
    free(region);
}

#ifdef __cplusplus
}
#endif
    
