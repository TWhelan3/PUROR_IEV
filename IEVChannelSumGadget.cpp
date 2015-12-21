//Unwrap3DGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Float 3D Array (Phase Data) -> // Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer]
//Output ImageHeader->Float 3D Array (Filtered Phase Data) -> Float 3D Array (Unfiltered Phase Data) -> [MetaContainer] 

#include "IEVChannelSumGadget.h"


namespace Gadgetron{

int IEVChannelSumGadget::process_config(ACE_Message_Block* mb)
{
	ISMRMRD::IsmrmrdHeader hdr;
        ISMRMRD::deserialize(mb->rd_ptr(),hdr);
	numEchos=hdr.encoding[0].encodingLimits.contrast().maximum +1; //number of echos is one more than highest numbers (0-based)
	this->msg_queue()->high_water_mark(128);//This helps with memory. It's not a hard limit though. 
	echoTimes=hdr.sequenceParameters.get().TE.get();//should be doing checks, structures are optional


	for(int e=0;e<numEchos;e++)
		echoTimes[e]/=1000; //still unclear if TEs are stored in seconds or milliseconds, but latter seems more common

	num_slices=hdr.encoding[0].reconSpace.matrixSize.z; //number of slices (will this always work?)

	yres=hdr.encoding[0].reconSpace.matrixSize.x; //match my (and MATLABs) unfortunate convention 
	xres=hdr.encoding[0].reconSpace.matrixSize.y;
	num_ch=hdr.acquisitionSystemInformation.get().receiverChannels();

	//Since we cannot connect to server to get series numbers, log locally
	std::string reader;//buffer for lines read in
	std::string found_id;
	
	if(hdr.studyInformation.is_present())
		studyInstanceUID = hdr.studyInformation.get().studyInstanceUID.get();
	else
	{
		dcmGenerateUniqueIdentifier(generatedStudyUID, SITE_STUDY_UID_ROOT);
		studyInstanceUID=std::string(generatedStudyUID);
	}
	//try to open file
	int pos;
	
	std::fstream log("../series_number_log");
	series_id_offset=99;
	//see if number is there
	if(log.is_open())
	{
		do{

			std::getline(log, reader);
			pos=reader.find_first_of(":");	
			found_id=reader.substr(0,pos);

			if(found_id.compare(studyInstanceUID)==0){//earlier series were part of this study
				series_id_offset=atoi(reader.substr(pos+1, std::string::npos).c_str());
			}
		}while(!log.eof());
		log.close();
	}
		
	//if file doesn't exist, it will be created	
	log.open("../series_number_log", std::fstream::app | std::fstream::out);
	log<<studyInstanceUID<<":"<<series_id_offset+output_phase.value()+output_LFS.value()<<std::endl;
	log.close();


	freq_ptr=new float[xres*yres*num_ch*numEchos];
		
	//unfiltered_phase_ptr=new float[xres*yres*num_ch*numEchos];

	hdr_ptr=new ISMRMRD::ImageHeader[numEchos];
	
	attributes=new ISMRMRD::MetaContainer[numEchos];

	return GADGET_OK;
}
int IEVChannelSumGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	GadgetContainerMessage<hoNDArray< float > > *unfiltered_unwrapped_msg_ptr =     AsContainerMessage<hoNDArray<float>>(m1->cont());

	GadgetContainerMessage<hoNDArray< float > > *filtered_unwrapped_msg_ptr =   AsContainerMessage<hoNDArray<float>>(unfiltered_unwrapped_msg_ptr->cont());
	GadgetContainerMessage<ISMRMRD::MetaContainer> *meta;
	
	static int c=0;	
	int e;
	int echo = m1->getObjectPtr()->contrast;
	float inv_echo_time;

	if(!filtered_unwrapped_msg_ptr || !unfiltered_unwrapped_msg_ptr)
	{
		GERROR("Wrong types received in IEVChannelSumGadget. Filtered and unfiltered phase expected.\n");
		return GADGET_FAIL;
	}
	
	 meta = AsContainerMessage<ISMRMRD::MetaContainer>(filtered_unwrapped_msg_ptr->cont());
	
	float* filtered_phase_ptr= filtered_unwrapped_msg_ptr->getObjectPtr()->get_data_ptr();
	
	m1->getObjectPtr()->channels=1; //yes?
	inv_echo_time=1/echoTimes[echo];//to avoid millions of divisions per slice

	//memcpy(unfiltered_phase_ptr+yres*xres*num_ch*echo, unfiltered_unwrapped_msg_ptr->getObjectPtr()->get_data_ptr(), xres*yres*num_ch*sizeof(float));
	
	for (int i = 0; i < xres*yres*num_ch; i++) 
		freq_ptr[echo*xres*yres*num_ch+i] = filtered_phase_ptr[i]*inv_echo_time;

	hdr_ptr[echo]=*(m1->getObjectPtr());
	
	if(meta)
	{
		meta->getObjectPtr()->append("StudyInstanceUID", studyInstanceUID.c_str());//may be in xml header, may not be, in that case put it in xml so it can get to dicom
		char TE[10];
		sprintf(TE, "%f", echoTimes[echo]*1000);
		meta->getObjectPtr()->append("TE", TE);

		attributes[echo]=*(meta->getObjectPtr());
	}
	unfiltered_unwrapped_msg_ptr->release();//all data has been copied
	if(echo==(numEchos-1))
	{	
		
		float* weights= new float[xres*yres*num_ch];
		float** channel_weights= new float* [num_ch];
		float* to_normalize = new float[xres*yres];
		int ch;
		
		if(iev.value()==int(IEV::YES))//just to allow this to work without IEV/make it very easy to compare
		{
			#pragma omp parallel //expanded parallel --- to allow sample to be allocated once
			{
				float* sample = new float[numEchos];
				int start = omp_get_thread_num()/omp_get_num_threads()*xres*yres*num_ch;
				int end = (omp_get_thread_num()+1)/omp_get_num_threads()*xres*yres*num_ch;
				for(int i =start; i <end; i++)
				{
					/////////
					for(int j = 0; j < numEchos; j++)
						sample[j]=freq_ptr[i+j*xres*yres*num_ch];     //assuming all(6-10 at at time) pages can be held in memory, this isn't terrible
					
					weights[i]=stdev(sample, numEchos);		//find standard deviation between echoes
					/////				
				}
				delete[] sample;
			}
			
			#pragma omp parallel for private(ch)
			for(ch = 0; ch < num_ch; ch++)
			{	
				float* temp_array;
		
			
				channel_weights[ch]=&weights[ch*xres*yres];
				medianFilter(channel_weights[ch], xres, yres);
				for (int i = 0; i < xres*yres; i++)
				{
				  channel_weights[ch][i]=1/(channel_weights[ch][i]+FLT_MIN);		//weight as inverse, 
				}
			
			}
			
			for (int i = 0; i < xres*yres; i++)
				to_normalize[i]	= 0;	

			for(int ch=0; ch< num_ch; ch++)		
				for (int i = 0; i < xres*yres; i++)
				{
				to_normalize[i]+=channel_weights[ch][i];
				}
	
			for(int ch=0; ch< num_ch; ch++)		
				for (int i = 0; i < xres*yres; i++)
				{
				channel_weights[ch][i]/=to_normalize[i];			//normalize weights
				}
			
		}
		else
		{
			#pragma omp parallel for private(ch)
			for(int ch=0; ch< num_ch; ch++)	
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

			//GadgetContainerMessage<ISMRMRD::ImageHeader>* h1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>(hdr_ptr[e]);
		
			
			//GadgetContainerMessage<hoNDArray< float > > *outimage = new GadgetContainerMessage<hoNDArray< float > >();
		
			/*try{outimage->getObjectPtr()->create(xres,yres);}	

			catch (std::runtime_error &err){
			GEXCEPTION(err,"Unable to create output image\n");
			return GADGET_FAIL;  
			}*/
			
			//float* output_ptr=outimage->getObjectPtr()->get_data_ptr();
			hdr_ptr[e].channels=1;
			hdr_ptr[e].contrast=e;
			hdr_ptr[e].data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;
			hdr_ptr[e].image_type = ISMRMRD::ISMRMRD_IMTYPE_PHASE;//There is no frequency image type
			hdr_ptr[e].slice= (hdr_ptr[e].image_index) % num_slices;//was hdr_ptr[e].image_index___-1_____) % num_slices before decrementor was added upstream
						
			if(output_phase.value())
			{
				//
				GadgetContainerMessage<ISMRMRD::ImageHeader>* phase_hdr = new GadgetContainerMessage<ISMRMRD::ImageHeader>(hdr_ptr[e]);
				//*(phase_hdr->getObjectPtr()) =*(hdr_ptr[e]->getObjectPtr());
				GadgetContainerMessage<hoNDArray< float > > *comb_phase_msg = new GadgetContainerMessage<hoNDArray< float > >();
				phase_hdr->getObjectPtr()->image_series_index=series_id_offset+1;
				try{comb_phase_msg->getObjectPtr()->create(xres,yres);}	

				catch (std::runtime_error &err){
				GEXCEPTION(err,"Unable to create output image\n");
				return GADGET_FAIL;  
				}
				float* output_ptr=comb_phase_msg->getObjectPtr()->get_data_ptr();
				phase_hdr->cont(comb_phase_msg);
				//	
				for (int i = 0; i < xres*yres; i++)
				{
				output_ptr[i]=filtered_phase_ptr[i]*channel_weights[0][i]; //instead of setting to 0 and adding first channel
				}
				for(int ch=1; ch< num_ch; ch++)	
					for (int i = 0; i < xres*yres; i++)
					{
						output_ptr[i]+=filtered_phase_ptr[xres*yres*ch+i]*channel_weights[ch][i];; //instead of setting to 0 and adding first channel
					}						
				//
				if(meta)
				{
				GadgetContainerMessage<ISMRMRD::MetaContainer>* meta = new GadgetContainerMessage<ISMRMRD::MetaContainer>(attributes[e]); 
						
				comb_phase_msg->cont(meta);	
				
				}
				//
				if (this->next()->putq(phase_hdr) == -1) {
				m1->release();
					GERROR("Unable to put collapsed images on next gadget's queue\n");
				return GADGET_FAIL; 
				}
				
			}
			if(output_LFS.value())
			{
				//
				GadgetContainerMessage<ISMRMRD::ImageHeader>* freq_hdr=new GadgetContainerMessage<ISMRMRD::ImageHeader>(hdr_ptr[e]);
				//*(freq_hdr->getObjectPtr()) =*(hdr_ptr[e]->getObjectPtr());
				GadgetContainerMessage<hoNDArray< float > > *comb_freq_msg = new GadgetContainerMessage<hoNDArray< float > >();
				freq_hdr->getObjectPtr()->image_series_index=series_id_offset+output_phase.value()+1;
				try{comb_freq_msg->getObjectPtr()->create(xres,yres);}	

				catch (std::runtime_error &err){
				GEXCEPTION(err,"Unable to create output image\n");
				return GADGET_FAIL;  
				}
				float* output_ptr=comb_freq_msg->getObjectPtr()->get_data_ptr();
				freq_hdr->cont(comb_freq_msg);
				freq_hdr->getObjectPtr()->image_type = 6;
				//
				for (int i = 0; i < xres*yres; i++)
					output_ptr[i]=freq_ptr[e*xres*yres*num_ch+i]*channel_weights[0][i]; //instead of setting to 0 and adding first channel
				for(int ch=1; ch< num_ch; ch++)		
					for (int i = 0; i < xres*yres; i++)
					{
						output_ptr[i]+=freq_ptr[e*xres*yres*num_ch+xres*yres*ch+i]*channel_weights[ch][i];
					}
				//
				if(meta)
				{
				GadgetContainerMessage<ISMRMRD::MetaContainer>* meta = new GadgetContainerMessage<ISMRMRD::MetaContainer>(attributes[e]); 
				meta->getObjectPtr()->set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_FREQMAP);
				//*(meta->getObjectPtr())=*(attributes[e]->getObjectPtr());
				comb_freq_msg->cont(meta);	
				}
				//
				if (this->next()->putq(freq_hdr) == -1) {
				//m1->release();
					GERROR("Unable to put collapsed images on next gadget's queue\n");
				return GADGET_FAIL; 
				}
			
			}
			
			
		}
		delete[] to_normalize;
		delete[] weights;
		delete[] channel_weights;
		//if(output.value()==int(OUTPUT::PHASE))
		//	delete[] unfiltered_phase_ptr;
		

	}
	


		
	return GADGET_OK;
}

void IEVChannelSumGadget::medianFilter(float* src_array, int xres, int yres)
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




GADGET_FACTORY_DECLARE(IEVChannelSumGadget)
}


