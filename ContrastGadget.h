//ContrastGadget.h
//Written by Tim Whelan 2015

#ifndef ContrastGadget_H_
#define ContrastGadget_H_

#include "PUROR.h"
enum PassOn {NO, YES};
namespace Gadgetron{

	class EXPORTGADGETSPUROR ContrastGadget:
	public Gadget1<ISMRMRD::ImageHeader>
	{
		public:
		GADGET_DECLARE(ContrastGadget);

		protected:
	 	virtual int process_config(ACE_Message_Block* mb);
		virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);
		GADGET_PROPERTY(pass_on, int, "Should individual echos be passed on?", 0);//0=No, 1=Yes

		int numEchos;
		int numSlices;
		int numPixels;

		std::vector<hoNDArray< float >*> sum_arrays;
		std::vector<std::string> types;
		//hoNDArray< float > *phase_sum, *mag_sum;
		//float *phase_sum_ptr, *mag_sum_ptr;
	};
}

#endif /* ContrastGadget_H_ */
