//HPFGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Float 3D Array (Phase Data) -> [MetaContainer]//
//Output ImageHeader->Float 3D Array (Filtered Phase Data) -> Float 3D Array (Unfiltered Phase Data) -> [MetaContainer] //

#include "HPFGadget.h"

namespace Gadgetron{

int HPFGadget::process_config(ACE_Message_Block* mb)
{
	ISMRMRD::IsmrmrdHeader hdr;
        ISMRMRD::deserialize(mb->rd_ptr(),hdr);
	FOVx=hdr.encoding[0].reconSpace.fieldOfView_mm.x/1000;
	FOVy=hdr.encoding[0].reconSpace.fieldOfView_mm.y/1000;
	//This is an optional part of the header and not sure if it will always be appropriate. # channels in acquisition may not be # channels in image gadget receives if they are already combined
	//for now leacing cres (actual image channels) separate
	num_ch=hdr.acquisitionSystemInformation.get().receiverChannels();
	this->msg_queue()->high_water_mark(128);//This helps with memory. It's not a hard limit though. 
	//this is awkward
	yres=hdr.encoding[0].reconSpace.matrixSize.x; //match my (and MATLABs) unfortunate convention 
	xres=hdr.encoding[0].reconSpace.matrixSize.y;

	//do FFT filter setup if that's what they want
		
	F.resize(xres);
	x2.resize(yres);
	y2.resize(xres);


	for(int i=0; i<yres; i++)
		x2[i]=pow((i-yres/2)*FOVx/yres,2);

	for(int i=0; i<xres; i++)
	{
		y2[i]=pow((i-xres/2)*FOVy/xres,2);
		F[i].resize(yres);
	}
		
	float sigma2=sigma.value()*sigma.value();
	for(int i=0; i<xres; i++)
	{
		for(int j=0; j<yres; j++)
		{
			F[i][j]=exp((y2[i]+x2[j])*-1/sigma2);//F[i][j]=exp((y2[i]+x2[j])*-1/.00009);
		}
	}

	

	return GADGET_OK;
}
int HPFGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	GadgetContainerMessage<hoNDArray<float > > *unfiltered_phase_msg =AsContainerMessage<hoNDArray<float >> (m1->cont());

	if(!unfiltered_phase_msg){
		GERROR("Phase array (float/single) expected and not found.\n");
		
		return GADGET_FAIL;
	}

	//GadgetContainerMessage<ISMRMRD::MetaContainer> *meta = AsContainerMessage<ISMRMRD::MetaContainer>(unfiltered_phase_msg->cont()->cont());
	GadgetContainerMessage<ISMRMRD::MetaContainer> *meta = AsContainerMessage<ISMRMRD::MetaContainer>(unfiltered_phase_msg->cont());
	ISMRMRD::Image<float> image;
	ISMRMRD::ImageHeader *header=m1->getObjectPtr();

	int cres = header->channels;
	float *pixel = unfiltered_phase_msg->getObjectPtr()->get_data_ptr();
	float* im_ptr;
	
	GadgetContainerMessage<hoNDArray< float > > *filtered_phase_msg;
	boost::shared_ptr< std::vector<size_t> > dims; 
	float *filtered_image;
	
	filtered_phase_msg = new GadgetContainerMessage<hoNDArray< float > >();
	dims = unfiltered_phase_msg->getObjectPtr()->get_dimensions();
		

	try{filtered_phase_msg->getObjectPtr()->create(dims.get());}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to create filtered image in Highpass Filter Gadget");
		return GADGET_FAIL;
	}


	filtered_image = filtered_phase_msg->getObjectPtr()->get_data_ptr();
	
	int ch;
	#pragma omp parallel for private(ch)	
	for(ch=0; ch<cres; ch++)
	{
		hoNDArray<std::complex<float> > *slicedata = new hoNDArray<std::complex<float> >();
		try{slicedata->create(header->matrix_size[0], header->matrix_size[1]);}
		catch (std::runtime_error &err){
		GEXCEPTION(err,"Creation Fail");
		//return GADGET_FAIL;
		}
		//This messes up the edges, how to get around it?
		std::complex<float> *slicedata_ptr=slicedata->get_data_ptr();

		for (int p = 0; p < xres*yres; p++) 
		{
		pixel[xres*yres*ch+p]=-pixel[xres*yres*ch+p];			//apply diagmagnetic convention
		slicedata_ptr[p] = std::complex<float>(pixel[xres*yres*ch+p],0);
		}
		//The actual FFT won't happen in parallel, but the rest of the data processing can
		//Mutexing happens within FFT library
		hoNDFFT<float>::instance()->ifftshift2D(*slicedata);
		hoNDFFT<float>::instance()->fft(slicedata,0);
		hoNDFFT<float>::instance()->fft(slicedata,1);
		hoNDFFT<float>::instance()->fftshift2D(*slicedata);			

		for(int i=0; i<xres; i++)
		{
			for(int j=0; j<yres; j++)
			{
				slicedata_ptr[i*yres+j]*=F[i][j];
			}
		}
		hoNDFFT<float>::instance()->ifftshift2D(*slicedata);
		hoNDFFT<float>::instance()->ifft(slicedata,1);
		hoNDFFT<float>::instance()->ifft(slicedata,0);
		hoNDFFT<float>::instance()->fftshift2D(*slicedata);
		
		for (int p =0; p < xres*yres; p++) 
		{
			filtered_image[xres*yres*ch+p] = pixel[xres*yres*ch+p] - slicedata_ptr[p].real(); //take out low freq signal, 
		}
		
			
		delete slicedata;

	}	
	//assumes a lot about the chain

	if(meta)
	{
		filtered_phase_msg->cont(meta);
	}
	unfiltered_phase_msg->cont(filtered_phase_msg);
	
	if (this->next()->putq(m1) == -1) {
		m1->release();
		GERROR("Unable to pass on filtered image.\n");
	    	return GADGET_FAIL;
	}



	return GADGET_OK;
}
GADGET_FACTORY_DECLARE(HPFGadget)
}












