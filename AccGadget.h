//AccGadget.h
//Written by Tim Whelan 2015
#ifndef AccGadget_H_
#define AccGadget_H_

#include "PUROR.h"
namespace Gadgetron{

  class EXPORTGADGETSPUROR AccGadget:
public Gadget1<ISMRMRD::ImageHeader>
    {

    public:
      GADGET_DECLARE(AccGadget);

    protected:
       	virtual int process_config(ACE_Message_Block* mb);
	virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);
	GADGET_PROPERTY(iev, int, "Weight channels by IEV?", 0); 
	GADGET_PROPERTY(norm, int, "Normalize channel(pixel) weights?", 0); 
	
	void medianFilter(float* src_array,int xres,int yres);
	int numEchos;
	int num_slices;
	std::vector<float> echoTimes;

      
     };
}

#endif /* AccGadget_H_ */


