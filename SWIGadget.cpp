//SWIGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Complex Float 3D Array (Image Data) -> [MetaContainer]-> //
//Output ImageHeader->Float 3D Array (Phase Data) -> Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer] ->//
#include "SWIGadget.h"
using namespace Gadgetron;

int SWIGadget::process_config(ACE_Message_Block* mb)
{
	boost::filesystem::path fname(magnitudefilename.value());

	if (!boost::filesystem::exists(fname))
	{
		GERROR("Missing magnitude source.\n");
		return GADGET_FAIL;
	}
	try {
		pDataset = new ISMRMRD::Dataset(fname.c_str(), groupname.value().c_str(), false);
	}
	catch (...) {
		GDEBUG("Unable to open: %s\n", fname.c_str());
	}

	ISMRMRD::IsmrmrdHeader hdr;
	ISMRMRD::deserialize(mb->rd_ptr(),hdr);

	numEchos=hdr.encoding[0].encodingLimits.contrast().maximum +1;
	numChan =hdr.acquisitionSystemInformation.get().receiverChannels();

	echoTimes=hdr.sequenceParameters.get().TE.get();//should be doing checks, structures are optional

	for(int e=0;e<numEchos;e++)
		echoTimes[e]/=1000; //still unclear if TEs are stored in seconds or milliseconds, but latter seems more common

	return GADGET_OK;
}

int SWIGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	//skeleton code. Missing a lot of things, wanted to get the file i/o and control flow figured out 
      
	GadgetContainerMessage< hoNDArray< float > >* m2 = AsContainerMessage< hoNDArray< float > > (m1->cont());

	if(!m2){
		GERROR("Phase expected and not found.\n");

		return GADGET_FAIL;
	}

	GadgetContainerMessage<ISMRMRD::MetaContainer> *meta = AsContainerMessage<ISMRMRD::MetaContainer>(m2->cont());

	if(!meta || !meta->getObjectPtr()->exists(GADGETRON_DATA_ROLE))
	{
		GERROR("Cannot identify frequency map. Passing on.");

		if (this->next()->putq(m1) < 0) { //pass to next gadget
			GERROR("Failed to pass on\n");
			m1->release();
			return GADGET_FAIL;
		}
	}

	std::string data_role = meta->getObjectPtr()->as_str(GADGETRON_DATA_ROLE);

	if(std::strcmp(data_role.c_str(), GADGETRON_IMAGE_FREQMAP)) //If this is not a frequency map pass it on
	{	//GADGETRON_IMAGE_FREQMAP is a macro for "FREQMAP". 

		if (this->next()->putq(m1) < 0) { //pass to next gadget
			GERROR("Failed to pass on magnitude\n");
			m1->release();
			return GADGET_FAIL;
		}
		return GADGET_OK;
	}

	yres = m2->getObjectPtr()->get_size(0);
	xres = m2->getObjectPtr()->get_size(1);
	int thisEcho = m1->getObjectPtr()->contrast;
	float thisTE = echoTimes[thisEcho];

	int c;
	ISMRMRD::Image<complex_float_t> image;

	GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
	GadgetContainerMessage<hoNDArray< float > > *cm2 = new GadgetContainerMessage<hoNDArray< float > >();
	hoNDArray<float> phase_mask;

	*cm1->getObjectPtr() = *m1->getObjectPtr();
	cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;

	cm1->cont(cm2);

	static int index = 0;
	try {
		pDataset->readImage(varname.value(), index++, image); //Assuming here files are ordered same way
	}
	catch (std::exception &ex) {
		GERROR("Error reading image %d from %s: %s\n",index, varname.value().c_str(), ex.what());
		return -1;
	}

	boost::shared_ptr< std::vector<size_t> > dims = m2->getObjectPtr()->get_dimensions();
	try{cm2->getObjectPtr()->create(yres, xres);}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to allocate array in SWI Gadget.\n");
		return GADGET_FAIL;
	}

	try{phase_mask.create(yres, xres);}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to allocate array in SWI Gadget.\n");
		return GADGET_FAIL;
	}

	std::complex<float>* src = image.getDataPtr();
	float*  phase_ptr = m2->getObjectPtr()->get_data_ptr();
	float*  dst = cm2->getObjectPtr()->get_data_ptr();

	for (int i = 0; i < xres*yres; i++)
	{
		phase_mask[i]=0;
		if(phase_ptr[i]*thisTE>0)
			phase_mask[i]=1;
		else
			phase_mask[i]=(phase_ptr[i]*thisTE+PI)/PI; //frequency already has *2PI component

		if(phase_mask[i]<0)
			phase_mask[i]=0;
	}

	for (int i = 0; i < xres*yres; i++) 
		dst[i]=abs(src[i]);

	for(c = 1; c<numChan; c++)
	{
		int ch_offset= xres*yres*c;

		for (int i = 0; i < xres*yres; i++) 
		{
			dst[i] += abs(src[i+ch_offset]);
		}
	}

	for (int i = 0; i < xres*yres; i++) 
		dst[i]*=std::pow(phase_mask[i],4);

	cm1->getObjectPtr()->channels=1;
	cm1->getObjectPtr()->image_type=ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE; //No SWI type
	cm1->getObjectPtr()->image_series_index++; //Make it one more than the LFS

	GadgetContainerMessage<ISMRMRD::MetaContainer>* SWImeta = new GadgetContainerMessage<ISMRMRD::MetaContainer>(*(meta->getObjectPtr()));

	SWImeta->getObjectPtr()->set(GADGETRON_DATA_ROLE, "SWI");

	cm2->cont(SWImeta);

	if (this->next()->putq(cm1) < 0) { //pass to next gadget
		GERROR("Failed to pass on magnitude\n");
		cm1->release();
		return GADGET_FAIL;
	}
	if (this->next()->putq(m1) < 0) { //pass to next gadget
		GERROR("Failed to pass on magnitude\n");
		m1->release();
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(SWIGadget)

