#ifndef DecImageIndex_H_
#define DecImageIndex_H_

#include "PUROR.h"
namespace Gadgetron{

	class EXPORTGADGETSPUROR DecImageIndex:
	public Gadget1<ISMRMRD::ImageHeader>
    	{
		public:
		GADGET_DECLARE(DecImageIndex);

		virtual int process_config(ACE_Message_Block*);
		virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);
	};
}
#endif
