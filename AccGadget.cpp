//Unwrap3DGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Float 3D Array (Phase Data) -> // Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer]
//Output ImageHeader->Float 3D Array (Phase Data) -> // Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer]

#include "AccGadget.h"


namespace Gadgetron{

int AccGadget::process_config(ACE_Message_Block* mb)
{
	ISMRMRD::IsmrmrdHeader hdr;
        ISMRMRD::deserialize(mb->rd_ptr(),hdr);
	numEchos=hdr.encoding[0].encodingLimits.contrast().maximum +1; //number of echos is one more than highest numbers (0-based)
	//this->msg_queue()->high_water_mark(128);
	echoTimes=hdr.sequenceParameters.get().TE.get();//should be doing checks, structures are optional
	num_slices=hdr.encoding[0].reconSpace.matrixSize.z; //number of slices (will this always work?)
	return GADGET_OK;
}
int AccGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	GadgetContainerMessage<hoNDArray< float > > *unwrapped_ptr =     AsContainerMessage<hoNDArray<float>>(m1->cont());
	
	//I think there should be a subtract first echo step involved that I am missing?
	static int c=0;	
	int e;
	int yres = unwrapped_ptr->getObjectPtr()->get_size(0);
	int xres = unwrapped_ptr->getObjectPtr()->get_size(1);
	int cres = unwrapped_ptr->getObjectPtr()->get_size(3);
	int image_series_index = m1->getObjectPtr()->image_series_index;
	int echo = m1->getObjectPtr()->contrast;
	float inv_echo_time;
	static hoNDArray< float >  *collapsed;
	static float* cptr;
	if(!unwrapped_ptr)
		GINFO("Wrong type\n");
	if(echo==0)
	{
		collapsed = new hoNDArray< float >();
		try{collapsed->create(xres,yres,cres,numEchos);}
			catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to create collapsed space");
		return GADGET_FAIL;
			}
		cptr=collapsed->get_data_ptr();
	}

	
	float* uptr= unwrapped_ptr->getObjectPtr()->get_data_ptr();
	
	m1->getObjectPtr()->channels=numEchos;
	inv_echo_time=1/echoTimes[echo];//to avoid millions of divisions per slice
	for (int i = 0; i < xres*yres*cres; i++) 
	{
		cptr[echo*xres*yres*cres+i] = uptr[i]*inv_echo_time;
				 
	}
	//need to save header (or possibly just the index, (contrast etc can come from that)) before deleting
	unwrapped_ptr->release();
	
	if(echo==(numEchos-1))
	{	
		
		float* weights= new float[xres*yres*cres];
		float** channel_weights= new float* [cres];
		float* to_normalize = new float[xres*yres];
		int ch,i;
		//GDEBUG("%d Saving Queue Length  %d My Queue Length\n", this->next()->msg_queue()->message_count(), this->msg_queue()->message_count());
		if(iev.value()==1)
		{
			#pragma omp parallel //expanded parallel --- to allow sample to be allocated once
			{
				float* sample = new float[numEchos];
				int start = omp_get_thread_num()/omp_get_num_threads()*xres*yres*cres;
				int end = (omp_get_thread_num()+1)/omp_get_num_threads()*xres*yres*cres;
				for(int i =start; i <end; i++)
				{
					/////////
					for(int j = 0; j < numEchos; j++)
						sample[j]=cptr[i+j*xres*yres*cres];     //assuming all(6-10 at at time) pages can be held in memory, this isn't terrible
					
					weights[i]=stdev(sample, numEchos);		//find standard deviation between echoes
					/////				
				}
				delete[] sample;
			}
			
			#pragma omp parallel for private(ch,i)
			for(ch = 0; ch < cres; ch++)
			{	
				float* temp_array;
		
			
				channel_weights[ch]=&weights[ch*xres*yres];
				medianFilter(channel_weights[ch],xres,yres);
				for (int i = 0; i < xres*yres; i++)
				{
				  channel_weights[ch][i]=1/(channel_weights[ch][i]*cres+FLT_MIN);		//weight as inverse, 
				}
				
			}
			if(norm.value()==1)
			{
				for (int i = 0; i < xres*yres; i++)
					to_normalize[i]	= 0;	

				for(int ch=0; ch< cres; ch++)		
					for (int i = 0; i < xres*yres; i++)
					{
					to_normalize[i]+=channel_weights[ch][i];
					}
		
				for(int ch=0; ch< cres; ch++)		
					for (int i = 0; i < xres*yres; i++)
					{
					channel_weights[ch][i]/=to_normalize[i];			//normalize weights
					}
			}
		}
		else
		{
			#pragma omp parallel for private(ch,i)
			for(int ch=0; ch< cres; ch++)	
			{
				channel_weights[ch]=&weights[ch*xres*yres];	
				for (int i = 0; i < xres*yres; i++)
				{
				channel_weights[ch][i]=1;			
				}
			}
		}
		for(e=0; e<numEchos; e++)
		{

			GadgetContainerMessage<ISMRMRD::ImageHeader>* h1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
		
			*h1->getObjectPtr()= *m1->getObjectPtr();


			GadgetContainerMessage<hoNDArray< float > > *outimage = new GadgetContainerMessage<hoNDArray< float > >();
		
			try{outimage->getObjectPtr()->create(xres,yres);}	

			catch (std::runtime_error &err){
			GEXCEPTION(err,"Unable to create output image");
			//return GADGET_FAIL;  omp cannot have two exits. Could flag this and return error after
			}
			
			float* output_ptr=outimage->getObjectPtr()->get_data_ptr();
			h1->getObjectPtr()->channels=1;
			h1->getObjectPtr()->contrast=e;
			h1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;
			h1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_PHASE;
			h1->getObjectPtr()->image_index = (h1->getObjectPtr()->image_index-1) % num_slices+1;//+e*num_slices;//not right, but distinguishes for now
		
			h1->getObjectPtr()->image_series_index = e;//h1->getObjectPtr()->image_series_index % num_slices+e*num_slices;//not right, but distinguishes for now
			h1->cont(outimage);


			for (int i = 0; i < xres*yres; i++)
				output_ptr[i]=cptr[e*xres*yres*cres+i]*channel_weights[0][i]; //instead of setting to 0 and adding first channel
			for(int ch=1; ch< cres; ch++)		
				for (int i = 0; i < xres*yres; i++)
				{
				output_ptr[i]+=cptr[e*xres*yres*cres+xres*yres*ch+i]*channel_weights[ch][i];
				}

			


			if (this->next()->putq(h1) == -1) {
			m1->release();
			GDEBUG("Unable to put collapsed images on next gadgets queue\n");
			//return GADGET_FAIL; omp cannot have two exits. Could flag this and return error after
			}
			

		
		}
	

		delete[] to_normalize;
		delete collapsed;
		delete[] weights;
		
		delete[] channel_weights;
		

	}
	


		
	return GADGET_OK;
}

void AccGadget::medianFilter(float* src_array,int xres,int yres)
{
	float array[25];
	float* filteredArray = new float[xres*yres];
	int column, row;
	int offset_index;
	int i,j,k;
	int x_offset, y_offset;
	float temp;
	int ext_xres=xres+4;
	int ext_yres=yres+4;
	k=0;

	float* dummyArray= new float[(ext_xres)*(ext_yres)];
	for(i=0;i<xres;i++)
		for(j=0;j<yres;j++)
			dummyArray[(i+2)*(ext_yres)+j+2]=src_array[i*yres+j];// 'real' part of array
	
	for(i=0;i<2*(ext_yres);i++)//pad left side with 0s
		dummyArray[i]=0;

	for(i=(ext_yres)*(ext_xres-2); i<(ext_yres)*(ext_xres);i++)//pad right side with 0s
		dummyArray[i]=0;

	for(i=2*(ext_yres);i<(ext_yres)*(ext_xres-2);i+=(ext_yres))//pad top and bottom with 0s
	{
		dummyArray[i]=0;
		dummyArray[i+1]=0;
		dummyArray[i+yres+2]=0;
		dummyArray[i+yres+3]=0;
	}
	
	for(column=0; column<xres; column++)//within each column
	{ 
		
		offset_index=(column+2)*(ext_yres)+2;
		
	 	for(row=0; row<yres; row++)// go down row by row
		{
		k=0;
		for(x_offset=-2;  x_offset<3;  x_offset++)
			for(y_offset=-2; y_offset<3; y_offset++) //get the values in the 5x5 square around the point
				array[k++]=dummyArray[offset_index+row+x_offset*(ext_yres)+y_offset]; 

		for(i=1; i<25; i++)
		{
			j=i;
			while(j>0 && array[j-1]>array[j])
			{
				temp=array[j];
				array[j]=array[j-1];
				array[j-1]=temp;
				j--;
			}
		}
		filteredArray[column*yres+row]=array[12];
		}
	}
	delete[] dummyArray;
	
	
	filteredArray[1]=(filteredArray[2]+filteredArray[yres+1])/2;
	filteredArray[yres]=(filteredArray[2*yres]+filteredArray[yres+1])/2;
	filteredArray[0]=(filteredArray[1]+filteredArray[yres])/2;
	
	filteredArray[yres-2]=(filteredArray[yres-3]+filteredArray[2*yres-2])/2;
	filteredArray[2*yres-1]=(filteredArray[2*yres-2]+filteredArray[3*yres-1])/2;
	filteredArray[yres-1]=(filteredArray[yres-2]+filteredArray[2*yres-1])/2;
	
	

	filteredArray[(xres-2)*yres]=(filteredArray[(xres-3)*yres]+filteredArray[(xres-2)*yres+1])/2;
	filteredArray[(xres-1)*yres+1]=(filteredArray[(xres-2)*yres+1]+filteredArray[(xres-1)*yres+2])/2;
	filteredArray[(xres-1)*yres]=(filteredArray[(xres-2)*yres]+filteredArray[(xres-1)*yres+1])/2;

	
	filteredArray[xres*yres-2]=(filteredArray[xres*yres-3]+filteredArray[(xres-1)*yres-2])/2;
	filteredArray[(xres-1)*yres-1]=(filteredArray[(xres-2)*yres-1]+filteredArray[(xres-1)*yres-2])/2;
	filteredArray[xres*yres-1]=(filteredArray[xres*yres-2]+filteredArray[(xres-1)*yres-1])/2;

	memcpy(src_array, filteredArray, xres*yres*sizeof(float));
	delete[] filteredArray;
}




GADGET_FACTORY_DECLARE(AccGadget)
}

