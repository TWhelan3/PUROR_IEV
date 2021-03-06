//AddMetaData.cpp
//Written by Tim Whelan 2015
//Input ImageHeader-> ???
//Output ImageHeader-> ???
#include "AddMetaData.h"
#include "DicomFinishGadget.h"
namespace Gadgetron{

int AddMetaData::process_config(ACE_Message_Block* mb)
{
	//taken from DicomFixMetaData
	ISMRMRD::IsmrmrdHeader hdr;
	deserialize(mb->rd_ptr(), hdr);

	//Fix incorrectly stored parameters
	//pixel spacing
	ISMRMRD::EncodingSpace r_space = hdr.encoding[0].reconSpace;
	pixel_spacing_X = r_space.fieldOfView_mm.x / r_space.matrixSize.x;
	pixel_spacing_Y = r_space.fieldOfView_mm.y / r_space.matrixSize.y;

	slice_spacing = r_space.fieldOfView_mm.z / r_space.matrixSize.z;

	return GADGET_OK;
}

int AddMetaData::process(GadgetContainerMessage<DcmFileFormat> * m1)
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

	if(!meta)
	{
		GERROR("No meta data found for DICOM\n");
		return GADGET_FAIL;
	}

	unsigned int BUFSIZE = 1024;
	char *buf = new char[BUFSIZE];
	const char* checkbuf;
	OFCondition status;
	DcmTagKey key;
	DcmDataset *dataset = dcm->getDataset();

	float rescaleIntercept;//=	meta->getObjectPtr()->as_double(GADGETRON_IMAGE_SCALE_OFFSET);
	float rescaleSlope;//=	meta->getObjectPtr()->as_double(GADGETRON_IMAGE_SCALE_RATIO);

	std::time_t rawtime;
	std::time(&rawtime);
	std::tm *timeinfo = std::localtime(&rawtime);

	if(meta->getObjectPtr()->exists("GADGETRON_IMAGE_SCALE_OFFSET") && meta->getObjectPtr()->exists("GADGETRON_IMAGE_SCALE_RATIO"))
	{
		rescaleIntercept=meta->getObjectPtr()->as_double(GADGETRON_IMAGE_SCALE_OFFSET);

		// Window Center
		key.set(0x0028, 0x1050);
		dataset->remove(key);

		// Window Width
		key.set(0x0028, 0x1051);
		dataset->remove(key);

		rescaleIntercept = -1.0*rescaleIntercept*rescaleSlope;

		rescaleSlope= 1.0/rescaleSlope;

		key.set(0x0028,0x1052);
		ACE_OS::snprintf(buf, BUFSIZE, "%f", rescaleIntercept);//
		write_dcm_string(dataset,key, buf);

		key.set(0x0028,0x1053);
		ACE_OS::snprintf(buf, BUFSIZE, "%f", rescaleSlope);//meta->getObjectPtr()->as_double("Intercept"));
		write_dcm_string(dataset,key, buf);


		key.set(0x0028, 0x0030); //Not sure this needs to be in this conditional
		ACE_OS::snprintf(buf, BUFSIZE, "%.6f\\%.6f", pixel_spacing_Y, pixel_spacing_X);
		write_dcm_string(dataset,key, buf);
	}

	key.set(0x0008,0x1030); //Study Description

	dataset->findAndGetString(key, checkbuf, false);
	if(checkbuf==NULL || !strcmp(checkbuf, "XXXXXXXX"))
	{
		ACE_OS::snprintf(buf, BUFSIZE, "%s", "Gadgetron^IEV");
		write_dcm_string(dataset,key, buf);
	}

	key.set(0x0008,0x103E); //Series Description

	std::string type;
	if(meta->getObjectPtr()->exists("GADGETRON_DATA_ROLE"))
		type=meta->getObjectPtr()->as_str(GADGETRON_DATA_ROLE);
	else
		type="MRI_Images";

	ACE_OS::snprintf(buf, BUFSIZE, "%s", type.c_str());
	write_dcm_string(dataset,key, buf);
	key.set(0x0020,0x0010);//Study ID
	dataset->findAndGetString(key, checkbuf, false);
	if(checkbuf==NULL || !strcmp(checkbuf, "XXXXXXXX"))
	{
		write_dcm_string(dataset,key, "1");
		// be sure to use the same one for all series you generate
	}

	//Study UID should be created in IEVChannelSumGadget.
	key.set(0x0020,0x000D);//Study UID
	dataset->findAndGetString(key, checkbuf, false);
	if(checkbuf==NULL || !strcmp(checkbuf, "XXXXXXXX")) //X's are put there by DicomFinish
	{
		if(meta->getObjectPtr()->exists("StudyInstanceUID"))
			write_dcm_string(dataset,key, meta->getObjectPtr()->as_str("StudyInstanceUID"));

		// If there's a better UID in meta, use it
		// be sure to use the same one for all series you generate
	}

	std::strftime(buf, 100, "%Y%m%d", timeinfo);

	key.set(0x0008,0x0020);//Study Date
	dataset->findAndGetString(key, checkbuf, false);
	if(checkbuf==NULL || !strcmp(checkbuf, "19000101"))
	{
		write_dcm_string(dataset,key, buf);
	}

	key.set(0x0008,0x0030);//Study Time
	dataset->findAndGetString(key, checkbuf, false);
	if(checkbuf==NULL || !strcmp(checkbuf, "121212"))
	{
		write_dcm_string(dataset,key, buf);
	}

	key.set(0x0008,0x0021);//Series Date
	if(!dataset->tagExistsWithValue(key))
	{
		write_dcm_string(dataset,key, buf);
	}

	key.set(0x0008,0x0012);//Instance Creation Date
	if(!dataset->tagExistsWithValue(key))
	{
		write_dcm_string(dataset,key, buf);
	}
	std::strftime(buf, 100, "%H%M%S", timeinfo);

	key.set(0x0008,0x0031);//Series Time
	if(!dataset->tagExistsWithValue(key))
	{
		write_dcm_string(dataset,key, buf);
	}

	key.set(0x0008,0x0013);//Instance Creation Time
	if(!dataset->tagExistsWithValue(key))
	{
		write_dcm_string(dataset,key, buf);
	}

	key.set(0x0018,0x0081);//Echo Time

	if(meta->getObjectPtr()->exists("TE"))
	{
		write_dcm_string(dataset,key, meta->getObjectPtr()->as_str("TE"));
	}

	delete[] buf;

	if(-1==this->next()->putq(m1))
	{
		m1->release();
		GERROR("Unable to pass on message\n");
		return GADGET_FAIL;
	}
	return GADGET_OK;
}
GADGET_FACTORY_DECLARE(AddMetaData)
}
