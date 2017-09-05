//DecImageIndex.cpp
//Written by Tim Whelan 2015
//Input ImageHeader-> ???
//Output ImageHeader-> ??? 
#include "DecImageIndex.h"

namespace Gadgetron{

int DecImageIndex::process_config(ACE_Message_Block* mb)
{
	//this gadget will be very fast, putting a high water mark would 
	//do nothing unless the next one is very slow and did not address backlog
	//in future might want to check read xml in this function
	
	return GADGET_OK;
}

int DecImageIndex::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	if(m1->getObjectPtr()->image_index==0)
		GERROR("Image index is already 0!\n");

	m1->getObjectPtr()->image_index=m1->getObjectPtr()->image_index-1;



	if(-1==this->next()->putq(m1))
	{
		m1->release();
		GERROR("Unable to pass on message\n");
		return GADGET_FAIL;
	}
	return GADGET_OK;
}
GADGET_FACTORY_DECLARE(DecImageIndex)
}
