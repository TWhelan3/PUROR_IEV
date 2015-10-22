#ifndef PUROR_H
#define PUROR_H


#include <cmath>
#include "Mask.h"
#include <vector>
#include "gadgetron_PUROR_export.h"
#include <math.h>
#include "Gadget.h"
#include "GadgetStreamInterface.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/dataset.h"
#include "ismrmrd/xml.h"
#include "mri_core_data.h"
#include "mri_core_def.h"
#include "omp.h"
#include <chrono>

#include <boost/filesystem.hpp>

const double PI=3.1416;
const double PI2=2*PI;//gets used a lot
float mean( std::vector<float> const& vec );
//template<typename T> T mean(T* sample, int length);
float mean(float *sample, int length);
double mean(double *sample, int length);
void unwrap(float* phase, int length);
void dePULM_1D(std::vector<float> &phase_1D, std::vector<int>& signal, std::vector<int>& index_ls);
void dePULM_1D_brain(std::vector<float> &phase_1D, std::vector<int>& index_ls);
std::vector<int>  dePULM_1D_ls(std::vector<float> &phase_1D,float phi_good, std::vector<int>& g_seg);
std::vector<int>  dePULM_1D_ls_brain(std::vector<float> &phase_1D,float phi_good);
void dePULM_2D_merging(std::vector<float>& line_phase,std::vector<int>& index_ls, std::vector<float>& re_phase,int res, float*, float*, float*);
float* filter(float* data, int xres, int yres);
float** filter2(std::complex<float>* temp, int xres, int yres);
void DericheSmoothing(float* pData, size_t N, float* mem, float sigma, size_t offset);//from hoNDImage_util.cpp
float stdev(float*, int);

#endif //PUROR_H
