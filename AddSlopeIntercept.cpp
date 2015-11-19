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

	static bool studyUIDmade=false;
	std::time_t rawtime;
         std::time(&rawtime);
         std::tm *timeinfo = std::localtime(&rawtime);
	

	rescaleIntercept = -1.0*rescaleIntercept*rescaleSlope;
	
	rescaleSlope= 1.0/rescaleSlope;
	
	key.set(0x0028,0x1052);
	ACE_OS::snprintf(buf, BUFSIZE, "%f", rescaleIntercept);//
	WRITE_DCM_STRING(key, buf);

	key.set(0x0028,0x1053);
	ACE_OS::snprintf(buf, BUFSIZE, "%f", rescaleSlope);//meta->getObjectPtr()->as_double("Intercept"));
	WRITE_DCM_STRING(key, buf);

	key.set(0x0008,0x1030); //Study Description
	if(!dataset->tagExistsWithValue(key))
	{
		ACE_OS::snprintf(buf, BUFSIZE, "%s", "Gadgetron^IEV");
		WRITE_DCM_STRING(key, buf);
	}


	key.set(0x0008,0x103E); //Series Description
	if(!dataset->tagExistsWithValue(key))
	{
		ACE_OS::snprintf(buf, BUFSIZE, "%s", "IEV phase");
		WRITE_DCM_STRING(key, buf);
	}

	key.set(0x0020,0x0010);//Study IUD
	if(!dataset->tagExistsWithValue(key))
	{
		if(studyUIDmade)
			WRITE_DCM_STRING(key, generatedStudyUID);
		else
		{
			dcmGenerateUniqueIdentifier(generatedStudyUID, SITE_STUDY_UID_ROOT);
			studyUIDmade=true;
		}
		
		// be sure to use the same one for all series you generate
	}
	
	std::strftime(buf, 100, "%Y%m%d", timeinfo);

	key.set(0x0008,0x0020);//Study Date
	if(!dataset->tagExistsWithValue(key))
	{
		WRITE_DCM_STRING(key, buf);
	}

	key.set(0x0008,0x0021);//Series Date
	if(!dataset->tagExistsWithValue(key))
	{
		WRITE_DCM_STRING(key, buf);
	}

	key.set(0x0008,0x0012);//Instance Creation Date		
	if(!dataset->tagExistsWithValue(key))
	{
		WRITE_DCM_STRING(key, buf);
	}	
	std::strftime(buf, 100, "%H%M%S.000", timeinfo);

	key.set(0x0008,0x0030);//Study Time
	if(!dataset->tagExistsWithValue(key))
	{
		WRITE_DCM_STRING(key, buf);
	}
	

	key.set(0x0008,0x0031);//Series Time
	if(!dataset->tagExistsWithValue(key))
	{
		WRITE_DCM_STRING(key, buf);
	}	


	key.set(0x0008,0x0013);//Instance Creation Time
	if(!dataset->tagExistsWithValue(key))
	{
		WRITE_DCM_STRING(key, buf);
	}	

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
