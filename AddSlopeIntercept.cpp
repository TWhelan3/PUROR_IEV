//AddSlopeIntercept.cpp
//Written by Tim Whelan 2015
//Input ImageHeader-> ???
//Output ImageHeader-> ??? 
#include "AddSlopeIntercept.h"
#include "DicomFinishGadget.h"
namespace Gadgetron{

int AddSlopeIntercept::process_config(ACE_Message_Block* mb)
{
	//this gadget will be very fast, putting a high water mark would 
	//do nothing unless the next one is very slow and did not address backlog
	//in future might want to check read xml in this function
	
	return GADGET_OK;
}

int AddSlopeIntercept::process(GadgetContainerMessage<DcmFileFormat> * m1)
{

	DcmFileFormat * dcm = m1->getObjectPtr();

	GadgetContainerMessage<std::string> * f;
	GadgetContainerMessage<ISMRMRD::MetaContainer>* meta;

	if(dcm)
	{	
		f=AsContainerMessage<std::string>(m1->cont());
	}
	else
	{	GERROR("No filename set for DICOM file\n");
		return GADGET_FAIL;
	}
	
	if(f)
	{	
		meta= AsContainerMessage<ISMRMRD::MetaContainer>(f->cont());
	}
	else
	{	
		GERROR("No meta data found for DICOM\n");
		return GADGET_FAIL;
	}

	unsigned int BUFSIZE = 1024;
        char *buf = new char[BUFSIZE];
	OFCondition status;
	DcmTagKey key;
	DcmDataset *dataset = dcm->getDataset();

	float rescaleIntercept=	meta->getObjectPtr()->as_double("Intercept");
	float rescaleSlope=	meta->getObjectPtr()->as_double("Scale");

	rescaleIntercept = -1.0*rescaleIntercept*rescaleSlope;
	
	rescaleSlope= 1.0/rescaleSlope;
	
	key.set(0x0028,0x1052);
	ACE_OS::snprintf(buf, BUFSIZE, "%f", rescaleIntercept);//
	WRITE_DCM_STRING(key, buf);

	key.set(0x0028,0x1053);
	ACE_OS::snprintf(buf, BUFSIZE, "%f", rescaleSlope);//meta->getObjectPtr()->as_double("Intercept"));
	WRITE_DCM_STRING(key, buf);

	delete[] buf;
	
	//add try catch 

	if(-1==this->next()->putq(m1))
	{
		m1->release();
		GERROR("Unable to pass on message\n");
		return GADGET_FAIL;
	}
	return GADGET_OK;
}
GADGET_FACTORY_DECLARE(AddSlopeIntercept)
}
