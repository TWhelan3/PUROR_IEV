//IEVChannelSumGadget.h
//Written by Tim Whelan 2015
#ifndef IEVChannelSumGadget_H_
#define IEVChannelSumGadget_H_

#include "PUROR.h"

namespace Gadgetron{

  class EXPORTGADGETSPUROR IEVChannelSumGadget:
public Gadget1<ISMRMRD::ImageHeader>
    {

	enum class OUTPUT {LFS, PHASE};
	enum class IEV {NO, YES};

    public:
      GADGET_DECLARE(IEVChannelSumGadget);

    protected:
       	virtual int process_config(ACE_Message_Block* mb);
	virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);
	GADGET_PROPERTY(iev, int, "Weight channels by IEV?", 1); 
	GADGET_PROPERTY(output, int, "What to output?",1); //0=LFS 1=Phase
	
	void medianFilter(float* src_array, int xres, int yres);//should be able to do this without passing yres/xres, but it causes a segfault
	int numEchos;
	int num_slices;
	int yres;
	int xres;
	int num_ch;
	std::vector<float> echoTimes;

	float* freq_ptr;
	float* unfiltered_phase_ptr;
	ISMRMRD::ImageHeader* hdr_ptr;
	ISMRMRD::MetaContainer* attributes;
	

	//add array for pixels and headers
	
      
     };
}

#endif /* IEVChannelSumGadget_H_ */


