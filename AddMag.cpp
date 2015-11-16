//AddMag.cpp
//Written by Tim Whelan 2015
//Throughput ImageHeader->Float 3D Array (Image Data) -> [MetaContainer]-> // 

//Output ImageHeader->Float 3D Array (Magnitude Data) -> [MetaContainer???] ->//

#include "AddMag.h"
using namespace Gadgetron;


int AddMag::process_config(ACE_Message_Block* mb)
{	
	this->msg_queue()->high_water_mark(128);//This helps with memory. It's not a hard limit though. 

	Gadget* g=this->get_controller()->find_gadget("LoadIsmrmrdDatasetImages");
	
	
	if(!g)
		return GADGET_FAIL;
	
	this->set_parameter("filename", g->get_string_value("filename")->c_str());
	this->set_parameter("groupname", g->get_string_value("groupname")->c_str());
	this->set_parameter("varname", g->get_string_value("varname")->c_str());
	
	boost::filesystem::path fname(filename.value());
	boost::filesystem::path dname(workingDirectory.value());
	
	fname=dname/fname;

	if (boost::filesystem::exists(fname)) {
         try {
            pDataset = new ISMRMRD::Dataset(fname.c_str(), groupname.value().c_str(), false);
         }
         catch (...) {
            GDEBUG("Unable to open: %s\n", fname.c_str());
         }
        }           
	
	ISMRMRD::IsmrmrdHeader hdr;
        ISMRMRD::deserialize(mb->rd_ptr(),hdr);
	numEchos= hdr.encoding[0].encodingLimits.contrast().maximum +1;
	numChan = hdr.acquisitionSystemInformation.get().receiverChannels();
	numSlices  =  hdr.encoding[0].reconSpace.matrixSize.z;

return GADGET_OK;
}


int AddMag::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	GadgetContainerMessage< hoNDArray< float > >* m2 = AsContainerMessage< hoNDArray< float > > (m1->cont());
	


	if(!m2){
		GERROR("Float image array expected and not found.\n");
		
		return GADGET_FAIL;
	}



	GadgetContainerMessage<ISMRMRD::MetaContainer> *meta = nullptr;
	ACE_Message_Block* m3 =m2->cont();
	while(m3!=NULL)
	{
		if(AsContainerMessage<ISMRMRD::MetaContainer>(m3)!=0)
		{
			meta=AsContainerMessage<ISMRMRD::MetaContainer>(m3);
			break;
		}
		m3=m3->cont();
	}
	if(meta)
	{
		meta->getObjectPtr()->set(GADGETRON_DATA_ROLE, "PHASE");
	}

	yres = m2->getObjectPtr()->get_size(0);
	xres = m2->getObjectPtr()->get_size(1);
	

	int c;
	ISMRMRD::Image<complex_float_t> image;

	GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
	GadgetContainerMessage<hoNDArray< float > > *cm2 = new GadgetContainerMessage<hoNDArray< float > >();


	*cm1->getObjectPtr() = *m1->getObjectPtr();
	cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;
	cm1->getObjectPtr()->image_type=ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;
	cm1->getObjectPtr()->contrast+=130; 
	cm1->getObjectPtr()->image_series_index+=130;
	
	cm1->cont(cm2);

	int index =(m1->getObjectPtr()->slice)*(numEchos)+m1->getObjectPtr()->contrast;
	      try {
            pDataset->readImage(varname.value(), index, image);
         }
         catch (std::exception &ex) {
		
            GERROR("Error reading image %d from %s: %s\n", index, varname.value().c_str(), ex.what());
            return -1;
         }

	boost::shared_ptr< std::vector<size_t> > dims = m2->getObjectPtr()->get_dimensions();
	try{cm2->getObjectPtr()->create(yres, xres);}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to allocate array in AddMag Gadget.\n");
		return GADGET_FAIL;
	}
	
	std::complex<float>* src = image.getDataPtr();
	float*  pm = m2->getObjectPtr()->get_data_ptr();
	float*  dst = cm2->getObjectPtr()->get_data_ptr();

	for (int i = 0; i < xres*yres; i++) 
		dst[i]=norm(src[i]);
	//#pragma omp parallel for private(c)
	for(c = 1; c<numChan; c++)
	{
		
		int ch_offset= xres*yres*c;

		for (int i = 0; i < xres*yres; i++) 
		{
		dst[i] += norm(src[i+ch_offset]);
		}
	}
	for (int i = 0; i < xres*yres; i++) 
		dst[i]=sqrtf(dst[i]);

	cm1->getObjectPtr()->channels=1;
	
	if(meta)
	{

		GadgetContainerMessage<ISMRMRD::MetaContainer> *newmeta = new GadgetContainerMessage<ISMRMRD::MetaContainer>();
		*(newmeta->getObjectPtr()) = *(meta->getObjectPtr());//actually copy instead of clone

		newmeta->getObjectPtr()->set(GADGETRON_DATA_ROLE, "MAG");
		cm2->cont(newmeta);
	}


	
	if (this->next()->putq(cm1) < 0) { //pass to next gadget
		GERROR("Failed to pass on magnitude\n");
		cm1->release();
	   return GADGET_FAIL;
	}
	if(this->next()->put(m1)<0){
	 //pass to next gadget
		GERROR("Failed to pass on phase\n");
		m1->release();
	   return GADGET_FAIL;
	}
	return GADGET_OK;
	
}

GADGET_FACTORY_DECLARE(AddMag)

