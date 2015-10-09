//SWIGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Complex Float 3D Array (Image Data) -> [MetaContainer]-> //
//Output ImageHeader->Float 3D Array (Phase Data) -> Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer] ->//
#include "SWIGadget.h"
using namespace Gadgetron;


int SWIGadget::process_config(ACE_Message_Block* mb)
{	
//this->msg_queue()->high_water_mark(128);//This helps with memory. It's not a hard limit though. 

	boost::filesystem::path fname(filename.value());
           
      if (boost::filesystem::exists(fname)) {
         try {
            pDataset = new ISMRMRD::Dataset(fname.c_str(), groupname.value().c_str(), true);
         }
         catch (...) {
            GDEBUG("Unable to open: %s\n", fname.c_str());
         }
      }

	ISMRMRD::IsmrmrdHeader hdr;
        ISMRMRD::deserialize(mb->rd_ptr(),hdr);

	numEchos=hdr.encoding[0].encodingLimits.contrast().maximum +1;
	numChan =hdr.acquisitionSystemInformation.get().receiverChannels();

return GADGET_OK;
}



int SWIGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{

    	//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
      
	GadgetContainerMessage< hoNDArray< float > >* m2 = AsContainerMessage< hoNDArray< float > > (m1->cont());
	
	if(!m2){
		GERROR("Complex image array expected and not found.\n");
		
		return GADGET_FAIL;
	}

	yres = m2->getObjectPtr()->get_size(0);
	xres = m2->getObjectPtr()->get_size(1);
	

	int c;
	ISMRMRD::Image<complex_float_t> image;

	GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
	GadgetContainerMessage<hoNDArray< float > > *cm2 = new GadgetContainerMessage<hoNDArray< float > >();


	*cm1->getObjectPtr() = *m1->getObjectPtr();
	cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;
	
	cm1->cont(cm2);

	int index =(m1->getObjectPtr()->image_index)*(numEchos)+1;
	       try {
            pDataset->readImage(varname.value(), index, image);
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
	
	std::complex<float>* src = image.getDataPtr();
	float*  pm = m2->getObjectPtr()->get_data_ptr();
	float*  dst = cm2->getObjectPtr()->get_data_ptr();

	for (int i = 0; i < xres*yres; i++) 
		dst[i]=abs(src[i]);
	//#pragma omp parallel for private(c)
	for(c = 1; c<numChan; c++)
	{
		
		int ch_offset= xres*yres*c;

		for (int i = 0; i < xres*yres; i++) 
		{
		dst[i] += abs(src[i+ch_offset]);
		}
	}
	for(int t=0; t<timetimes.value(); t++)
	{
	for (int i = 0; i < xres*yres; i++) 
		dst[i]*=pm[i];
	}
	cm1->getObjectPtr()->channels=1;
	cm1->getObjectPtr()->image_type=ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;
	
	 if (this->next()->putq(cm1) < 0) { //pass to next gadget
		GERROR("Failed to pass on magnitude\n");
		cm1->release();
	   return GADGET_FAIL;
	}
	m2->cont(NULL);//necessary so keep meta etc
	m1->release(); 
	return GADGET_OK;
	
}

GADGET_FACTORY_DECLARE(SWIGadget)

