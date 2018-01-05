#ifndef UnwrapGadget_H_
#define UnwrapGadget_H_

#include "PUROR.h"
#include <vector>
namespace Gadgetron{

	class EXPORTGADGETSPUROR UnwrapGadget:
	public Gadget1<ISMRMRD::ImageHeader>
	{
		public:
		GADGET_DECLARE(UnwrapGadget);

		protected:
		virtual int process_config(ACE_Message_Block*);
		virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);
		GADGET_PROPERTY(filename, std::string, "Name of ISMRMRD Dataset file", "output");
		GADGET_PROPERTY(savephase, int, "Save Before Passing On", 0);
		ISMRMRD::Dataset* dsToWrite;


		int xres;
		int yres;
		int num_ch;
		int num_echos;	

		float q_th; //might want to set this is config

	};
}

#endif /* UnwrapGadget_H_ */
