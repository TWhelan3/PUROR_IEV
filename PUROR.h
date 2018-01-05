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

float* filter(float* data, int xres, int yres);
float** filter2(std::complex<float>* temp, int xres, int yres);
void DericheSmoothing(float* pData, size_t N, float* mem, float sigma, size_t offset);//from hoNDImage_util.cpp
float stdev(float*, int);
void unwrap(std::vector<float>& phase);
void unwrap(float* phase, int length);
void unwrap_rows(std::vector<float>&, int, int);
void unwrap_columns(std::vector<float>&, int, int);
void shift_to_mean(std::vector<float>&, MaskData&, int, int, int);
void shift_to_mean_brain(std::vector<float>&, MaskData&, int, int, int);

void calc_quality_x(std::vector<float>&, std::vector<int>&, std::vector<float>&,int&,int&, int, int);

void calc_quality_y(std::vector<float>&, std::vector<int>&, std::vector<float>&,int&,int&, int, int);
void calc_quality_y_noref(std::vector<float>&, int,int, std::vector<float>&);
void center_x(std::vector<float>&, std::vector<float>&, std::vector<int>&, int,int,int,int);
void center_y(std::vector<float>&, std::vector<float>&, std::vector<int>&, int,int,int,int);

void final_compare(std::vector<float>&, std::vector<float>&, int, MaskData&,int,int);

void final_compare_brain(std::vector<float>&, std::vector<float>&, int, int, int);

void diff_x(std::vector<float>&,bool,MaskData&, int,int);
void diff_y(std::vector<float>&,bool, MaskData&, int, int);

void dePULM_1D(std::vector<float> &phase_1D, std::vector<int>& index_ls);
void dePULM_1D_brain(std::vector<float> &phase_1D, std::vector<int>& index_ls);
std::vector<int>  dePULM_1D_ls(std::vector<float> &,float,std::vector<int>&);
std::vector<int>  dePULM_1D_ls_brain(std::vector<float> &,float);
void dePULM_2D_merging(std::vector<float>&,std::vector<int>& , std::vector<float>&,int,std::vector<float>&, std::vector<float>& );
void dePULM_2D_itoh(std::vector<int>&, std::vector<float>&, int, int);
std::vector<float> dePULM_2D_mean_y(std::vector<float>&,MaskData&,int,int);

template<typename T> T mean(T*, int);
template<typename T> T mean(T *sample, int length)
{
	T mean;
	mean=0;
	for(int ii=0; ii<length; ii++)
		mean+=sample[ii];
	mean/=length;
	if(length==0){
		mean=0;//Not correct for all applications
	}
	return mean;
}
template<typename T> T mean(std::vector<T>);
template<typename T> T mean(std::vector<T> sample)
{
	return mean(sample.data(), sample.size());
}

#endif //PUROR_H
