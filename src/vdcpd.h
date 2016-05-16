#ifndef VDCPD_H
#define VDCPD_H 1

#ifdef __cplusplus
extern "C" {
#endif



extern void genVDCPD(const long dims[], int *pattern_cfl, 
              const int *feasiblePointsD, const float FOVRatio, 
              const float C, const int shapeOpt, const float initial_mindist,
              const float vd_exp, const int maxR);


#ifdef __cplusplus
}
#endif

#endif // VDCPD_H

