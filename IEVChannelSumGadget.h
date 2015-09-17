//IEVChannelSumGadget.h
//Written by Tim Whelan 2015
#ifndef IEVChannelSumGadget_H_
#define IEVChannelSumGadget_H_

#include "PUROR.h"
enum OUTPUT {LFS, PHASE};
namespace Gadgetron{

  class EXPORTGADGETSPUROR IEVChannelSumGadget:
public Gadget1<ISMRMRD::ImageHeader>
    {

    public:
      GADGET_DECLARE(IEVChannelSumGadget);

    protected:
       	virtual int process_config(ACE_Message_Block* mb);
	virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);
	GADGET_PROPERTY(iev, int, "Weight channels by IEV?", 0); 
	GADGET_PROPERTY(norm, int, "Normalize channel(pixel) weights?", 0); 
	GADGET_PROPERTY(output, int, "What to output?",0 ); //0=LFS 1=Phase
	
	void medianFilter(float* src_array,int xres,int yres);
	int numEchos;
	int num_slices;
	std::vector<float> echoTimes;

      
     };
}

#endif /* IEVChannelSumGadget_H_ */


