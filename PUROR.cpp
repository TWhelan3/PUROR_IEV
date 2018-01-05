#include "PUROR.h"

float indexed_mean(std::vector<float> full_sample, std::vector<int>& indices)
{
	std::vector<float> sample(indices.size());

	for(int ii=0; ii<indices.size(); ii++)
	{
		sample[ii]=full_sample[indices[ii]];
	}

	return mean(sample);
}
float stdev(float* sample, int length)
{
	
	int ii;
	float var=0;
	float datamean; 
	float diff;		
	datamean=mean(sample, length);

	for(ii=0; ii<length; ii++)
	{	
		diff=sample[ii]-datamean;
		var+=diff*diff;
	}
	var/=(length-1);
	
	return sqrt(var);



}
void unwrap(std::vector<float>& phase)
{
	unwrap(phase.data(), phase.size());
}
void unwrap(float* phase, int length)
{
	float dps;
	std::vector<float> diff(length);
	std::vector<float> dp_corr(length);

	dp_corr[0]=0;

	for(int ii=1; ii<length; ii++)		
		diff[ii]=phase[ii]-phase[ii-1]; //incremental phase step

	for(int ii=1; ii<length; ii++)
	{
		if(fabs(diff[ii])<PI)
		{
			dp_corr[ii]=0;
		}
		else
		{
			dps=(float)remainder((diff[ii]+PI),PI2)-PI; // should be equivalent step in[-pi, pi], but negative values for first arg give negative results

			if(dps==-PI && diff[ii]>0)//check pi vs -pi, compiler will complain this is unsafe, but does get hit
				dps=PI;

			if(dps<(-1*PI))//shouldn't happen, but does because of the way mod/remainder is done
				dps+=PI2;

			dp_corr[ii]=dps-diff[ii];//correction amount
		}

		dp_corr[ii]+=dp_corr[ii-1];//integrate correction
		phase[ii]+=dp_corr[ii];//apply correction
	}

	return;
}

void unwrap_rows(std::vector<float>& toUnwrap, int xres, int yres)
{
	std::vector<float> phase_tmp(xres);

	for(int ii=0; ii<yres; ii++)//Unwrap each row, alternating directions
	{
		if(ii%2==0)
		{
			for(int jj=0; jj<xres; jj++)
				phase_tmp[xres-1-jj]=toUnwrap[jj*yres+ii];

			unwrap(phase_tmp);

			for(int jj=0; jj<xres; jj++)
				toUnwrap[jj*yres+ii]=phase_tmp[xres-1-jj];

		}
		else
		{
			for(int jj=0; jj<xres; jj++)
				phase_tmp[jj]=toUnwrap[jj*yres+ii];

			unwrap(phase_tmp);

			for(int jj=0; jj<xres; jj++)
				toUnwrap[jj*yres+ii]=phase_tmp[jj];
		}
	}
}
void unwrap_columns(std::vector<float>& toUnwrap, int xres, int yres)
{
	std::vector<float> phase_tmp(yres);

	for(int ii=0; ii<xres; ii++)//Unwrap each column alternating directions
	{
		if(ii%2==0)
		{
			for(int jj=0; jj<yres; jj++)
				phase_tmp[yres-jj-1]=toUnwrap[ii*yres+jj];

			unwrap(phase_tmp);

			for(int jj=0; jj<yres; jj++)
				toUnwrap[ii*yres+jj]=phase_tmp[yres-jj-1];
		}
		else
		{
			for(int jj=0; jj<yres; jj++)
				phase_tmp[jj]=toUnwrap[ii*yres+jj];

			unwrap(phase_tmp);

			for(int jj=0; jj<yres; jj++)
				toUnwrap[ii*yres+jj]=phase_tmp[jj];
		}
	}
}

void dePULM_1D(std::vector<float> &phase_1D, std::vector<int>& signal,std::vector<int>& index_ls)
{
	//This function reconnects segments identified by dePULM_1D_ls
	float test, diff_1D;
	float inter_diff_l;
	int ii, kk, length, tmp_1D;
	int k;
	std::vector<float> tmp_smooth(phase_1D);
	kk=0;
	std::vector<float> phase_1D_tmp(phase_1D);

	//for(ii=0; ii<phase_1D.size(); ii++)
		//phase_1D_tmp[ii]=phase_1D[ii];

	if(index_ls.size()>2)//if at least one segment identified
	{
		for(ii=1; ii<index_ls.size()-2; ii+=2)//ii is end of a segment, ii+1 is beginning of next
		{
			length=index_ls[ii+1]-index_ls[ii];
			if(length <=6)
			{
				test = phase_1D[index_ls[ii + 1]] + phase_1D[index_ls[ii + 1] + 1] + phase_1D[index_ls[ii + 1] + 2] - phase_1D[index_ls[ii]] - phase_1D[index_ls[ii] - 1] - phase_1D[index_ls[ii] - 2];
				test = test/3; 
				//test=difference between the average of the first three points of next segment and last three points of this segment

				if(fabs(test)>PI)//represents a pole
				{
					test=round(test/(PI2))*PI2;
					kk=index_ls[ii]+2;
					for(k=kk; k<phase_1D.size(); k++)//shift phase data to eliminate discontinuity
					{
						phase_1D[k]-=test;
					}
				}
			}
			if(kk!=0)//has kk been changed(was there a pole)
			{
				inter_diff_l=phase_1D[kk]-phase_1D[kk-1];
				inter_diff_l=inter_diff_l>0?floor(inter_diff_l/PI2):ceil(inter_diff_l/PI2);//'fix' (biased round to zero)
				if(fabs(inter_diff_l)>0)//check for wrap in beginning of segment, fix it
					phase_1D[kk]-=PI2*round(inter_diff_l);

			}
		}


		tmp_1D=index_ls[0];//eliminate wraps at beginning of segment

		if(tmp_1D!=0)
		{
			if(tmp_1D==1)
			{
				diff_1D=phase_1D[tmp_1D]-phase_1D[0];
				phase_1D[0]=phase_1D_tmp[0] + PI2*round(diff_1D/PI2);
			}
			else
			{
				for(ii=0; ii<tmp_1D; ii++)
					tmp_smooth[ii] = phase_1D_tmp[ii];

				diff_1D=phase_1D[tmp_1D]-mean(tmp_smooth.data(), tmp_1D);//awkward here, using part of a vector
				diff_1D= PI2*round(diff_1D/PI2);
				for(ii=0; ii<tmp_1D; ii++)
					phase_1D[ii] = tmp_smooth[ii] + diff_1D;
			}

		}

		tmp_1D = index_ls[index_ls.size()-1];//eliminate wraps at end of segment

		if (tmp_1D != phase_1D.size()-1)
		{

			if (tmp_1D == phase_1D.size() - 2)
			{
				diff_1D = phase_1D[tmp_1D] - phase_1D_tmp[tmp_1D + 1];
				phase_1D[tmp_1D + 1] =phase_1D_tmp[tmp_1D+1]+PI2*round(diff_1D/PI2);  ///////    
			}
			else
			{
				for(ii=tmp_1D+1; ii<phase_1D.size(); ii++)
					tmp_smooth[ii-tmp_1D-1] = phase_1D_tmp[ii];
				diff_1D = phase_1D[tmp_1D] - mean(tmp_smooth.data(), phase_1D.size()-tmp_1D-1);
				diff_1D= PI2*round(diff_1D/PI2);
				for(ii=tmp_1D+1; ii<phase_1D.size(); ii++)
					phase_1D[ii]=phase_1D_tmp[ii]+diff_1D;
			}
		}
	}

	return;
}

float* filter(float* temp, int xres, int yres)//mean 2d filter returning 1d pointer for images
{
	float a,b,c;
	int ii,jj;
	float *data;
	data= new float[yres*xres];
	
	for(ii=1; ii<yres-1; ii++)
	{
		b=temp[ii-1]+temp[ii]+temp[ii+1];
		c=temp[ii-1+yres]+temp[ii+yres]+temp[ii+1+yres];

		for(jj=1; jj<xres-1; jj++)
		{
			a=b; b=c;
			c=temp[ii-1+yres*(jj+1)]+temp[ii+yres*(jj+1)]+temp[ii+1+yres*(jj+1)];
			data[ii+yres*jj]=(a+b+c)*.1111;//don't need to divide
		}
	}

	for(ii=1; ii<yres-1; ii++)//Get borders
	{
		b=temp[ii-1]+temp[ii]+temp[ii+1];
		c=temp[ii-1+yres]+temp[ii+yres]+temp[ii+1+yres];


		data[ii]=(b+c)*.1111;//don't need to divide

		a=temp[ii-1+yres*(xres-2)]+temp[ii+yres*(xres-2)]+temp[ii+1+yres*(xres-2)];
		b=temp[ii-1+yres*(xres-1)]+temp[ii+yres*(xres-1)]+temp[ii+1+yres*(xres-1)];

		data[ii+yres*jj]=(a+b)*.1111;//don't need to divide

	}

	b=temp[0]+temp[1];

	c=temp[0+yres]+temp[1+yres];

	data[0]=(b+c)*.1111;//don't need to divide

	for(jj=1; jj<xres-1; jj++)
	{
		a=b; b=c;
		c=temp[yres*(jj+1)]+temp[yres*(jj+1)+1];

		data[yres*jj]=(a+b+c)*.1111;//don't need to divide

	}

	data[yres*(xres-1)]=(b+c)*.1111;//don't need to divide

	b=temp[yres-1]+temp[yres-2];

	c=temp[yres-1+yres]+temp[yres-2+yres];

	data[yres-1]=(b+c)*.1111;//don't need to divide 

	for(jj=1; jj<xres-1; jj++)
	{
		a=b; b=c;
		c=temp[yres-1+yres*(jj+1)]+temp[yres-2+yres*(jj+1)];

		data[ii+yres*jj]=(a+b+c)*.1111;//don't need to divide

	}

	data[yres-1+yres*(xres-1)]=(b+c)*.1111;//don't need to divide

	return data;
}
float** filter2(std::complex<float>* input, int dim1, int dim2)//mean 2d filter returning absolute from complex returning 2d pointer (for mask creation)
{
	std::complex<float> a,b,c;
	int ii,jj;
	float **data;
	std::complex<float> **temp;
	data= new float*[dim1];
	temp= new std::complex<float>* [dim1];
	for(ii=0; ii<dim1; ii++)
		data[ii]= new float[dim2];

	for(ii=0; ii<dim1; ii++)
		temp[ii]= input+dim2*ii;
	

	for(ii=1; ii<dim1-1; ii++)
	{
		b=temp[ii-1][0]+temp[ii][0]+temp[ii+1][0];
		c=temp[ii-1][1]+temp[ii][1]+temp[ii+1][1];

		for(jj=1; jj<dim2-1; jj++)
		{
			a=b; b=c;
			c=temp[ii-1][jj+1]+temp[ii][jj+1]+temp[ii+1][jj+1];
			data[ii][jj]=abs(a+b+c)*.1111;//don't need to divide
		}
	}

	for(ii=1; ii<dim1-1; ii++)//Get borders
	{
		b=temp[ii-1][0]+temp[ii][0]+temp[ii+1][0];
		c=temp[ii-1][1]+temp[ii][1]+temp[ii+1][1];


		data[ii][0]=abs(b+c)*.1111;//don't need to divide

		a=temp[ii-1][dim2-2]+temp[ii][dim2-2]+temp[ii+1][dim2-2];
		b=temp[ii-1][dim2-1]+temp[ii][dim2-1]+temp[ii+1][dim2-1];

		data[ii][jj]=abs(a+b)*.1111;//don't need to divide

	}

	b=temp[0][0]+temp[1][0];

	c=temp[0][1]+temp[1][1];

	data[0][0]=abs(b+c)*.1111;//don't need to divide

	for(jj=1; jj<dim2-1; jj++)
	{
		a=b; b=c;
		c=temp[0][jj+1]+temp[1][jj+1];

		data[0][jj]=abs(a+b+c)*.1111;//don't need to divide

	}

	data[0][dim2-1]=abs(b+c)*.1111;//don't need to divide

	b=temp[dim1-1][0]+temp[dim1-2][0];

	c=temp[dim1-1][1]+temp[dim1-2][1];

	data[dim1-1][0]=abs(b+c)*.1111;//don't need to divide //Oct8 reswitch

	for(jj=1; jj<dim2-1; jj++)
	{
		a=b; b=c;
		c=temp[dim1-1][jj+1]+temp[dim1-2][jj+1];

		data[ii][jj]=abs(a+b+c)*.1111;//don't need to divide

	}

	data[dim1-1][dim2-1]=abs(b+c)*.1111;//don't need to divide

	return data;
}

void DericheSmoothing(float* pData, size_t N, float* mem, float sigma, size_t offset)//from hoNDImage_util.cpp
{
   

    // following the note of http://en.wikipedia.org/wiki/Deriche_edge_detector

    float alpha = (1.4105/sigma); // this value 1.4105 is from equation 37 of ref [1]
    float e_alpha = ( std::exp( (double)(-alpha) ) );
    float e_alpha_sqr = e_alpha*e_alpha;
    float k = ( (1-e_alpha)*(1-e_alpha) ) / ( 1 + 2*alpha*e_alpha - e_alpha_sqr );

    float a1 = k;
    float a2 = k * e_alpha * (alpha-1);
    float a3 = k * e_alpha * (alpha+1);
    float a4 = -k * e_alpha_sqr;

    float b1 = 2 * e_alpha;
    float b2 = -e_alpha_sqr;

    // compute the left to right filtering and the right to left filtering
    // for the speed, just use the zero boundary condition
    // TODO: try out other boundary conditions
    float* forward = mem;
    float* reverse = mem + N;

    if ( offset == 0 )
    {
        forward[0] = a1 * pData[0];
        reverse[N-1] = 0;

        size_t ii;

        if ( N > 1 )
        {
            forward[1] = a1 * pData[1] + a2*pData[0] + b1 * forward[0];
            reverse[N-2] = a3 * pData[N-1] + b1 * reverse[N-1];

            for ( ii=2; ii<N; ii++ )
            {
                forward[ii] = (a1*pData[ii] + a2*pData[ii-1]) + (b1*forward[ii-1] + b2*forward[ii-2]);
                reverse[N-1-ii] = (a3*pData[N-ii] + a4*pData[N-ii+1]) + (b1*reverse[N-ii] + b2*reverse[N-ii+1]);
            }
        }

        // Gadgetron::math::add(N, forward, reverse, pData);

        for ( ii=0; ii<N; ii++ )
        {
            pData[ii] -= (forward[ii]+ reverse[ii]);//get difference and negate
        }
    }
    else
    {
        forward[0] = a1 * pData[0];
        reverse[N-1] = 0;

        if ( N > 1 )
        {
            forward[1] = a1 * pData[offset] + a2*pData[0] + b1 * forward[0];
            reverse[N-2] = a3 * pData[(N-1)*offset] + b1 * reverse[N-1];

            size_t ii;
            for ( ii=2; ii<N; ii++ )
            {
                forward[ii] = (a1*pData[ii*offset] + a2*pData[(ii-1)*offset]) + (b1*forward[ii-1] + b2*forward[ii-2]);
                reverse[N-1-ii] = (a3*pData[(N-ii)*offset] + a4*pData[(N-ii+1)*offset]) + (b1*reverse[N-ii] + b2*reverse[N-ii+1]);
            }

            for ( ii=0; ii<N; ii++ )
            {
                pData[ii*offset] -= (forward[ii] + reverse[ii]);//get difference and negate
            }
        }
    }
}

void shift_to_mean(std::vector<float>& phase, MaskData& md, int lineLength, int lineCount, int dirFlag)
{
	std::vector<float> phase_sample(lineLength);
	std::vector<float> mean_unwrap(lineCount);
	std::vector<float> mean_connect(lineCount);
	std::vector<int>   index_ls;
	std::vector<std::vector<int>> connect;
	std::vector<std::vector<int>> signal;
	std::vector<std::vector<int>> segments;
	float correction, phi_good;

	int midline=round(lineCount/2)-1;
	int readStride, ptrStride;

	bool fix_first;

	float* phase_offset;

	if(dirFlag == 1) //shift to mean (x)
	{
		readStride = lineCount;
		ptrStride = 1;
		connect = md.connectXH;
		signal = md.signalX;
		segments = md.segX;
	}
	else //shift to mean (y)
	{
		readStride = 1;
		ptrStride = lineLength;
		connect = md.connectYH;
		signal = md.signalY;
		segments = md.segY;
	}

	for(int index=0; index<lineCount; index++)
	{
		float myMean;
		std::vector<int> good_segments=segments[index];
		int numSignalPoints=signal[index].size();
		int numConnectedPoints=connect[index].size();
		std::vector<int> lineSignal = signal[index];

		phase_offset=phase.data()+index*ptrStride;

		if(numSignalPoints>lineLength/4)//if mask covers more than 1/4 of the column
		{
			for(int ii=0; ii<numSignalPoints; ii++)
				phase_sample[ii]=phase_offset[lineSignal[ii]*readStride];

			phi_good=PI/2;
			index_ls = dePULM_1D_ls(phase_sample,phi_good,good_segments);//Find segments in phase_sample with phase gradient < phi_good

			dePULM_1D(phase_sample,index_ls);//Reconnect these segments
		}//took out else because it's fine to operate on original phase, copies can be made outside this function

		//calculate connect mean

		if(numSignalPoints==0)
		{
			if(index==0)
			{
				myMean = 0;
				fix_first = true;
			}
			else
				myMean=mean_connect[index-1];
		}
		else
		{
			//Find mean of each line
			if(numConnectedPoints<=3)
				myMean=indexed_mean(phase_sample, signal[index]);
			else
				myMean=indexed_mean(phase_sample, connect[index]);

			//bring the average of the line to [-pi,pi]
			correction=PI2*round(myMean/PI2);

			for(int ii=0; ii<lineLength; ii++)
				phase_sample[ii]-=correction;
			if(numConnectedPoints<=3)
				myMean-=correction; //don't need to recalculate mean to know
			else
				myMean=indexed_mean(phase_sample, connect[index]);
		}

		for(int ii=0; ii<lineLength; ii++)
				phase_offset[ii*readStride]=phase_sample[ii];

		mean_connect[index]=myMean;
	}

	if(fix_first)
		mean_connect[0]=mean_connect[1];

	mean_unwrap=mean_connect;

	//unwrap the mean values before global shifting the data
	unwrap(mean_unwrap);
	correction=round(mean_unwrap[midline]/PI2)*PI2;
	for(int ii=0; ii<lineCount; ii++)
		mean_unwrap[ii]-=correction;

	//shift phase data
	for(int index = 0; index<lineCount; index++)
	{
		float diff_test = mean_unwrap[index] - mean_connect[index];
		std::vector<int>& lineSignal = signal[index];
		if (fabs(diff_test) > PI)
		{
			phase_offset = phase.data() + index*ptrStride;
			correction=PI2*round(diff_test/(PI2));
			for(int ii=0; ii<lineSignal.size(); ii++)
				phase_offset[lineSignal[ii]*readStride] += correction;
		}
	}
}
void shift_to_mean_brain(std::vector<float>& phase, MaskData& md, int lineLength, int lineCount, int dirFlag)//brain version assumes full signal
{
	std::vector<float> phase_sample(lineLength);
	std::vector<float> mean_unwrap(lineCount);
	std::vector<float> mean_connect(lineCount);
	std::vector<int>   index_ls;
	std::vector<std::vector<int>> connect;
	float diff_test, correction, phi_good;

	int midline=round(lineCount/2)-1;
	int readStride, ptrStride;

	float* phase_offset;

	if(dirFlag == 1) //shift to mean (x)
	{
		readStride = lineCount;	//each iteration of loop below reads a line orthogonal to memory set up
		ptrStride = 1;
		connect = md.connectXH;
	}
	else //shift to mean (y)
	{
		readStride = 1;
		ptrStride = lineLength;	//each iteration of loop below reads a line contiguous in memory
		connect = md.connectYH;
	}

	for(int index=0; index<lineCount; index++)
	{
		float myMean;
		phase_offset=phase.data()+index*ptrStride;

		for(int ii=0; ii<lineLength; ii++)
			phase_sample[ii]=phase_offset[ii*readStride];

		phi_good=PI/2;

		index_ls = dePULM_1D_ls_brain(phase_sample,phi_good);//Find segments in phase with phase gradient < phi_good

		dePULM_1D_brain(phase_sample,index_ls);//Reconnect these segments

		//Find mean of each line
		if(connect[index].size()<=3)
			myMean=mean(phase_sample);
		else
			myMean=indexed_mean(phase_sample, connect[index]);

		//bring the average of the line to [-pi,pi]
		correction=PI2*round(myMean/PI2);

		for(int ii=0; ii<lineLength; ii++)
			phase_sample[ii]-=correction;

		if(connect[index].size()<=3)
			myMean-=correction; //don't need to recalculate mean to know
		else
			myMean=indexed_mean(phase_sample, connect[index]);

		for(int ii=0; ii<lineLength; ii++)
				phase_offset[ii*readStride]=phase_sample[ii];

		mean_connect[index]=myMean;
	}

	mean_unwrap=mean_connect;

	//unwrap the mean values before global shifting the data
	unwrap(mean_unwrap);
	correction=round(mean_unwrap[midline]/PI2)*PI2;
	for(int ii=0; ii<lineCount; ii++)
		mean_unwrap[ii]-=correction;

	//shift phase data
	for(int index = 0; index<lineCount; index++)
	{
		diff_test = mean_unwrap[index] - mean_connect[index];

		if (fabs(diff_test) > PI)
		{
			correction=PI2*round(diff_test/(PI2));

			for(int ii=0; ii<lineLength; ii++)
				phase[index*ptrStride+ii*readStride] +=correction;
		}
	}
}

std::vector<float> dePULM_2D_mean_y(std::vector<float>& phase_original, MaskData& md, int xres, int yres)
{
	//shift_to_mean clobbers input phase, so this function makes a copy and uses that
	std::vector<float> phase_new(phase_original);

	shift_to_mean_brain(phase_new, md, xres, yres, 2);

	return phase_new;
}

void calc_quality_y(std::vector<float>& phase_x, std::vector<int>& mask, std::vector<float> &quality_y, int &xy_start_dw, int &xy_start_up, int xres, int yres)//eqv to matlabs quality_ch
{
	int len_dw, len_up;

	int flag_start, flag_end;
	int g_start, g_end;

	std::vector<int> pointsToUse;

	float q_th = 0.75;  //*Changed over time?


	for(int row_index=0; row_index<yres; row_index++)
	{
		len_dw=0; len_up=0;
		if(row_index!=0 && row_index!=(yres-1))
		{
			for(int ii=0;ii<xres; ii++)//find spots in t_mask which has a point above and below to check quality
			{
				//if(t_mask[(row_index+1)*xres+ii] && t_mask[(row_index-1)*xres+ii] && t_mask[row_index*xres+ii])
				if(mask[row_index+1+yres*ii] && mask[row_index-1+yres*ii] && mask[row_index+yres*ii])
					pointsToUse.push_back(ii);
			}
			for(int ii=0; ii<pointsToUse.size(); ii++)
			{
				//count number of jumps between rows
				if((fabs(phase_x[yres*pointsToUse[ii]+row_index]-phase_x[yres*pointsToUse[ii]+row_index+1]))>PI)
					len_dw++;
			}
			for(int ii=0; ii<pointsToUse.size(); ii++)
			{
				//count number of jumps between rows
				if((fabs(phase_x[yres*pointsToUse[ii]+row_index]-phase_x[yres*pointsToUse[ii]+row_index-1]))>PI)
					len_up++;
			}
			if(!pointsToUse.empty())//find quality of row
			{
				quality_y[row_index]=(float)(pointsToUse.size()-len_up)*(float)(pointsToUse.size()-len_dw)/(float)(pointsToUse.size()*pointsToUse.size());
			}
			else
			{
				quality_y[row_index]=-1;
			}
		}
		if(row_index==0)//special case
		{
			len_dw=0;
			for(int ii=0;ii<xres; ii++)
			{
				if(mask[row_index+1+yres*ii] && mask[row_index+yres*ii])
					pointsToUse.push_back(ii);
			}

			for(int ii=0; ii<pointsToUse.size(); ii++)
			{
				if(fabs(phase_x[yres*pointsToUse[ii]] -  phase_x[yres*pointsToUse[ii]+1])>PI);
				len_dw++;
			}

			if(!pointsToUse.empty())
				quality_y[row_index]=(float)(pointsToUse.size()-len_dw)/((float)pointsToUse.size());
			else
				quality_y[row_index]=-1;
		}
		if(row_index==(yres-1))//special case
		{
			len_up=0;
			for(int ii=0;ii<xres; ii++)
			{
				if(mask[row_index-1+yres*ii] && mask[row_index+yres*ii])
					pointsToUse.push_back(ii);

			}
			for(int ii=0; ii<pointsToUse.size(); ii++)
			{
				//count number of jumps between rows
				if((fabs(phase_x[yres*pointsToUse[ii]+row_index]-phase_x[yres*pointsToUse[ii]+row_index-1]))>PI)
					len_up++;
			}

			if(!pointsToUse.empty())
				quality_y[row_index]=(float)(pointsToUse.size()-len_up)/(pointsToUse.size());
			else
				quality_y[row_index]=-1;

		}

		pointsToUse.clear();
	}

	xy_start_dw=0;
	xy_start_up=0;
	flag_start=1;
	flag_end=2;

	for(int ii=0; ii<yres; ii++)//find highest quality strip
	{
		if(quality_y[ii]>q_th && flag_start==1)//quality above threshold, mark start
		{
			g_start=ii;
			flag_start=0;
			flag_end=1;
		}
		if(quality_y[ii]<q_th && flag_end ==1)//quality below threshold, if started, mark end
		{
			g_end=ii-1;
			flag_start=1;
			flag_end=0;
		}
		if(ii==(yres-1) && flag_end ==1)//last row, if started, save end
		{
			if (quality_y[ii]>=q_th)
				g_end=ii;
			else
				g_end=g_start;
			flag_end=0;
		}
		if(flag_end==0 && ii>0)//check if this strip is bigger than earlier strips
			if((g_end-g_start)>=(xy_start_up-xy_start_dw))
			{
				xy_start_dw=g_start;
				xy_start_up=g_end;
			}
	}
}
void calc_quality_y_noref(std::vector<float>& phase_tmp, int xres, int yres, std::vector<float>& quality_br)
{
	/*This function works the same as dePULM_2D_quality, but assumes no mask.*/
	float q_th;
	int flag_start, flag_end;
	int g_start, g_end;
	float len_dw, len_up;/*floats for multiplication etc*/

	quality_br.resize(yres+2);/*Holds start/end in last two spots*/

	for(int index_y=0; index_y<yres; index_y++)
	{
		len_dw=0; len_up=0;
		if(index_y!=0 && index_y!=(yres-1))
		{
			len_up=len_dw;/*not safe if this loop runs on multiple threads.*/
			len_dw=0;
			for(int ii=0; ii<xres; ii++)
			{
				/*count number of jumps*/

				if((fabs(phase_tmp[index_y+yres*ii]-phase_tmp[index_y+1+yres*ii]))>PI)
					len_dw++;
			}
			/*get quality of row*/
			quality_br[index_y]=((float)xres-len_up)*((float)xres-len_dw)/(float)(xres*xres);

		}
		if(index_y==0)/*special case*/
		{
			len_dw=0;
			for(int ii=0; ii<xres; ii++)
			{
				if((fabs(phase_tmp[index_y+yres*ii]-phase_tmp[index_y+1+yres*ii]))>PI)
					len_dw++;
			}

			quality_br[index_y]=((float)xres-len_dw)/(float)(xres);
		}
		if(index_y==(yres-1))/*special case*/
		{
			len_up=len_dw;
			quality_br[index_y]=((float)xres-len_up)/(float)(xres);
		}
	}

	quality_br[yres]=1;
	quality_br[yres+1]=1;
	flag_start=1;
	flag_end=2;
	q_th=0.9;

	for(int ii=0; ii<yres; ii++)/*find highest quality strip*/
	{
		if(quality_br[ii]>q_th && flag_start==1)/*quality above threshold, mark start*/
		{
			g_start=ii;
			flag_start=0;
			flag_end=1;
		}
		if(quality_br[ii]<q_th && flag_end ==1)/*quality below threshold, if started, mark end*/

		{
			g_end=ii-1;
			flag_start=1;
			flag_end=0;
		}
		if(ii==(yres-1) && flag_end ==1)/*last row, if started, save end*/
		{
			if (quality_br[ii]>=q_th)
				g_end=ii;
			else
				g_end=g_start;
			flag_end=0;
		}
		if(flag_end==0 && ii>0)/*check if this strip is bigger than earlier strips*/
			if((g_end-g_start)>=(quality_br[yres+1]-quality_br[yres]))
			{
				quality_br[yres]=g_start;
				quality_br[yres+1]=g_end;
			}
	}
}

void calc_quality_x(std::vector<float>& phase_y, std::vector<int>& mask, std::vector<float> &quality_x,int &xy_start_L, int &xy_start_R, int xres, int yres)//eqv to matlabs quality_ch_y
{
	int len_dw, len_up;

	int flag_start, flag_end;
	int g_start, g_end;

	std::vector<int> pointsToUse;

	float q_th = 0.75;  //*Changed over time?

	for(int col_index=0; col_index<xres; col_index++)
	{
		len_dw=0; len_up=0;
		if(col_index!=0 && col_index!=(xres-1))
		{
			for(int ii=0;ii<yres; ii++)//find spots in mask which has a point left and right to check quality
			{
				if(mask[(col_index+1)*yres+ii] && mask[(col_index-1)*yres+ii] && mask[col_index*yres+ii])
					pointsToUse.push_back(ii);
			}
			for(int ii=0; ii<pointsToUse.size(); ii++)
			{
				//count jumps between columns
				if(fabs(phase_y[pointsToUse[ii]+yres*col_index] -  phase_y[pointsToUse[ii]+yres*(col_index+1)])>PI)
					len_dw++;
				if(fabs(phase_y[pointsToUse[ii]+yres*col_index] - phase_y[pointsToUse[ii]+yres*(col_index-1)])>PI)
					len_up++;
			}

			if(!pointsToUse.empty())//get quality of column
				quality_x[col_index]=(float)(pointsToUse.size()-len_up)*(pointsToUse.size()-len_dw)/(pointsToUse.size()*pointsToUse.size());

			else

				quality_x[col_index]=-1;

		}
		if(col_index==0)//special case of above loop
		{
			len_dw=0;
			for(int ii=0;ii<yres; ii++)
			{
				if(mask[(col_index+1)*yres+ii] && mask[col_index*yres+ii])
					pointsToUse.push_back(ii);
			}
			for(int ii=0; ii<pointsToUse.size(); ii++)
			{
				if(fabs(phase_y[pointsToUse[ii]+yres*col_index] -  phase_y[pointsToUse[ii]+yres*(col_index+1)])>PI)
				len_dw++;
			}

			if(!pointsToUse.empty())
				quality_x[col_index]=(pointsToUse.size()-len_dw)/(float)(pointsToUse.size());
			else
				quality_x[col_index]=-1;

		}
		if(col_index==(xres-1))//special case of above loop
		{
			len_up=0;
			for(int ii=0;ii<yres; ii++)
			{
				if(mask[(col_index-1)*yres+ii] && mask[col_index*yres+ii])
					pointsToUse.push_back(ii);
			}
			for(int ii=0; ii<pointsToUse.size(); ii++)
			{
				if(fabs(phase_y[pointsToUse[ii]+yres*col_index] -  phase_y[pointsToUse[ii]+yres*(col_index-1)])>PI)
				len_up++;
			}

			if(!pointsToUse.empty())
				quality_x[col_index]=((float)pointsToUse.size()-len_up)/(pointsToUse.size());
			else
				quality_x[col_index]=-1;
		}
		pointsToUse.clear();

	}

	xy_start_L=0;
	xy_start_R=0;
	flag_start=1;
	flag_end=2;

	for(int ii=0; ii<xres; ii++)//find highest quality strip
	{
		if(quality_x[ii]>q_th && flag_start==1)//quality above threshold, mark start
		{
			g_start=ii;
			flag_start=0;
			flag_end=1;
		}
		if(quality_x[ii]<q_th && flag_end ==1)//quality below threshold, if started, mark end
		{
			g_end=ii-1;
			flag_start=1;
			flag_end=0;
		}
		if(ii==(xres-1) && flag_end ==1)//last column, if started, save end
		{
			if (quality_x[ii]>=q_th)
				g_end=ii;
			else
				g_end=g_start;
			flag_end=0;
		}
		if(flag_end==0 && ii>0)//check if this strip is bigger than earlier strips
			if((g_end-g_start)>=(xy_start_R-xy_start_L))
			{
				xy_start_L=g_start;
				xy_start_R=g_end;
			}
	}

}

void center_x(std::vector<float>& phase_x, std::vector<float>& phase_y, std::vector<int>& mask, int xy_start_dw, int xy_start_up, int xres, int yres)
{
	float diffmean;

	std::vector<int> index_l;
	std::vector<int> index_s;
	std::vector<float> diff;

	index_l.reserve(yres);
	index_s.reserve(yres);
	diff.reserve(yres);
	for(int col_index=0; col_index<xres; col_index++)
	{
		for(int ii=xy_start_dw;ii<=xy_start_up; ii++)
		{
			if(mask[col_index*yres+ii])
				index_s.push_back(ii);

		}
		for(int ii=0;ii<yres; ii++)
		{
			if(std::abs(phase_x[ii+yres*col_index])>0)
				index_l.push_back(ii);
		}
		if(index_s.empty())//no good strip, use whole column
		{
			for(int ii=0;ii<yres; ii++)
			{
				if(mask[col_index*yres+ii])
					index_s.push_back(ii);
			}

			for(int ii=0;ii<index_s.size(); ii++)//find column mean difference
				diff.push_back(phase_x[index_s[ii]+yres*col_index]-phase_y[index_s[ii]+yres*col_index]);

			diffmean=mean(diff);
			diffmean=(PI2)*round(diffmean/(PI2));
			for(int ii=0;ii<index_l.size(); ii++)//adjust column
				phase_y[index_l[ii]+yres*col_index]+=diffmean;
		}
		else//cross reference with good strip
		{
			//find column mean difference
			for(int ii=0;ii<index_s.size(); ii++)
				diff.push_back(phase_x[index_s[ii]+yres*col_index]-phase_y[index_s[ii]+yres*col_index]);

			diffmean=mean(diff);
			if(fabs(diffmean)>PI)//adjust column
			{
				diffmean=(PI2)*round(diffmean/(PI2));
				for(int ii=0;ii<index_l.size(); ii++)
					phase_y[index_l[ii]+yres*col_index]+=diffmean;
			}
		}

		index_l.clear();
		index_s.clear();
		diff.clear();
	}
}

void center_y(std::vector<float>& phase_y, std::vector<float>& phase_x, std::vector<int>& t_mask, int xy_start_L, int xy_start_R, int xres, int yres)
{
	float diffmean;

	std::vector<int> index_l;
	std::vector<int> index_s;
	std::vector<float> diff;

	index_l.reserve(xres);
	index_s.reserve(xres);
	diff.reserve(xres);
	for(int row_index=0; row_index<yres; row_index++)
	{
		for(int ii=xy_start_L;ii<=xy_start_R; ii++)
		{
			if(t_mask[row_index*xres+ii])
				index_s.push_back(ii);
		}
		for(int ii=0;ii<xres; ii++)
		{
			if(std::abs(phase_y[row_index+yres*ii])>0)
				index_l.push_back(ii);
		}

		if(index_s.empty())//no good strip use whole row
		{
			for(int ii=0;ii<xres; ii++)
			{
				if(t_mask[row_index*xres+ii])
					index_s.push_back(ii);
			}
			for(int ii=0;ii<index_s.size(); ii++)//find mean row difference
				diff.push_back(phase_y[row_index+yres*index_s[ii]]-phase_x[row_index+yres*index_s[ii]]);

			diffmean=mean(diff);
			diffmean=(PI2)*round(diffmean/(PI2));

			for(int ii=0;ii<index_l.size(); ii++)//adjust row
				phase_x[row_index+yres*index_l[ii]]+=diffmean;

		}
		else//cross reference with good strip
		{
			for(int ii=0;ii<index_s.size(); ii++)//find mean row difference
				diff.push_back(phase_y[row_index+yres*index_s[ii]]-phase_x[row_index+yres*index_s[ii]]);

			diffmean=mean(diff);

			if(fabs(diffmean)>PI)
			{	diffmean=(PI2)*round(diffmean/(PI2));
				//adjust row

				for(int ii=0;ii<index_l.size(); ii++)
					phase_x[row_index+yres*index_l[ii]]+=diffmean;
			}
		}
		index_l.clear();
		index_s.clear();
		diff.clear();
	}
}

void final_compare(std::vector<float>& phase_x, std::vector<float>& phase_y, int iterations, MaskData& md, int xres, int yres)
{
	float seg_phi=PI2;

	//Preallocate space for speed.
	std::vector<float> diff_tmp;
	std::vector<float> seg_mean_shift;

	std::vector<int> g_seg;
	std::vector<int> index_s;

	std::vector<float> ref_tmp;
	std::vector<float> trl_tmp;
	std::vector<int> index_ls_TR;
	if(xres>=yres)
	{
		diff_tmp.resize(xres);
		seg_mean_shift.resize(xres);
	}
	else
	{
		diff_tmp.resize(yres);
		seg_mean_shift.resize(yres);
	}
	for(int loopcount=0;loopcount<iterations; loopcount ++)
	{
		ref_tmp.resize(yres);
		trl_tmp.resize(yres);
		seg_phi/=2;
		for(int index_x=0; index_x<xres; index_x++)
		{
			if(!(md.segY[index_x].empty()) || (md.signalY[index_x].empty()))
			{
				g_seg=md.segY[index_x];
				index_s=md.signalY[index_x];

				for(int ii=0; ii<index_s.size(); ii++)
				{
					ref_tmp[ii]=phase_x[index_s[ii]+yres*index_x];
					trl_tmp[ii]=phase_y[index_s[ii]+yres*index_x];
				}

				index_ls_TR=dePULM_1D_ls(trl_tmp, seg_phi, g_seg);//find phase segments with gradient < seg_phi
				dePULM_2D_merging(trl_tmp, index_ls_TR, ref_tmp, index_s.size(), seg_mean_shift, diff_tmp);//Cross reference those segmentsmd.Ylengths[index_x][0] was yres

				for(int ii=0; ii<index_s.size(); ii++)
					phase_y[index_s[ii]+index_x*yres]=trl_tmp[ii];

			}

		}
	//////////////////////
		ref_tmp.resize(xres);
		trl_tmp.resize(xres);
		for(int index_y=0; index_y<yres; index_y++)
		{
			if(!(md.segX[index_y].empty()) || (md.signalX[index_y].empty()))
			{
				g_seg=md.segX[index_y];
				index_s=md.signalX[index_y];

				for(int ii=0; ii<index_s.size(); ii++)
				{
					ref_tmp[ii]=phase_y[index_y+index_s[ii]*yres];
					trl_tmp[ii]=phase_x[index_y+index_s[ii]*yres];
				}

				index_ls_TR=dePULM_1D_ls(trl_tmp,seg_phi,g_seg);//Identify phase segments with gradient < seg_phi

				dePULM_2D_merging(trl_tmp,index_ls_TR,ref_tmp,index_s.size(), seg_mean_shift, diff_tmp);//Cross reference those strips

				for(int ii=0; ii<index_s.size(); ii++)
					phase_x[index_y+index_s[ii]*yres]=trl_tmp[ii];
			}

		}
	}
}

void final_compare_brain(std::vector<float>& phase_x, std::vector<float>& phase_y, int iterations, int xres, int yres)
{
	float seg_phi=PI2;

	//Preallocate space for speed.
	std::vector<float> diff_tmp;
	std::vector<float> seg_mean_shift;

	std::vector<float> ref_tmp;
	std::vector<float> trl_tmp;
	std::vector<int> index_ls_TR;
	if(xres>=yres)
	{
		diff_tmp.resize(xres);
		seg_mean_shift.resize(xres);
	}
	else
	{
		diff_tmp.resize(yres);
		seg_mean_shift.resize(yres);
	}
	for(int loopcount=0;loopcount<iterations; loopcount++)
	{
		ref_tmp.resize(yres);
		trl_tmp.resize(yres);
		seg_phi/=2;
		for(int index_x=0; index_x<xres; index_x++)
		{
				memcpy(ref_tmp.data(), phase_x.data()+index_x*yres, sizeof(float)*yres);
				memcpy(trl_tmp.data(), phase_y.data()+index_x*yres, sizeof(float)*yres);

				index_ls_TR=dePULM_1D_ls_brain(trl_tmp, seg_phi);//find phase segments with gradient < seg_phi
				dePULM_2D_merging(trl_tmp, index_ls_TR, ref_tmp,yres, seg_mean_shift, diff_tmp);//Cross reference those segmentsmd.Ylengths[index_x][0] was yres

				memcpy(phase_y.data()+index_x*yres,trl_tmp.data(), sizeof(float)*yres);

		}
	//////////////////////
		ref_tmp.resize(xres);
		trl_tmp.resize(xres);
		for(int index_y=0; index_y<yres; index_y++)
		{
				for(int ii=0; ii<xres; ii++)
				{
					ref_tmp[ii]=phase_y[index_y+ii*yres];
					trl_tmp[ii]=phase_x[index_y+ii*yres];
				}

				index_ls_TR=dePULM_1D_ls_brain(trl_tmp,seg_phi);//Identify phase segments with gradient < seg_phi

				dePULM_2D_merging(trl_tmp, index_ls_TR,ref_tmp,xres, seg_mean_shift, diff_tmp);//Cross reference those strips

				for(int ii=0; ii<xres; ii++)
					phase_x[index_y+ii*yres]=trl_tmp[ii];
		}
	}
}

void diff_x(std::vector<float>& phase_x,bool fullsignal, MaskData& md, int xres, int yres)
{
	int index_y; //this index is declared because first/last need to be addressed differently
	int len_dw,len_up;
	float diff_up, diff_dw, quality_br_tmp;

	std::vector<std::vector<float> > trl_dw(yres, std::vector<float>(xres));

	std::vector<bool> bad_up(xres);
	std::vector<bool> bad_dw(xres);

	for(int ii=0; ii<yres; ii++)
	{
		for(int jj=0; jj<xres; jj++)
			trl_dw[ii][jj]=phase_x[ii+jj*yres];
	}

	index_y=0; len_dw=0;
	//Same as below with index_y = 0, so no check -1
	for(int ii=0; ii<xres; ii++)
	{
		if(fabs(phase_x[index_y+ii*yres] -  phase_x[index_y+1+ii*yres])>PI)
			len_dw++;
	}

	if(len_dw!=0)
	{
		quality_br_tmp=(float)((xres-len_dw)/(xres));

		if(quality_br_tmp<0.95)
		{
			for(int ii=0; ii<xres; ii++)
			{
				diff_dw=phase_x[index_y+1+ii*yres] - phase_x[index_y+ii*yres];
				phase_x[index_y+ii*yres]+=PI2*(round(diff_dw/PI2));
			}
		}
	}

	for(index_y=1; index_y<(yres-1);index_y++)
	{
		len_up=0;
		len_dw=0;

		for(int ii=0;ii<xres; ii++)
		{
			bad_up[ii]=0;
			bad_dw[ii]=0;

			if(fabs(trl_dw[index_y][ii] -  trl_dw[index_y+1][ii])>PI)
			{
				len_dw++;
				bad_dw[ii]=true;
			}

			if(fabs(trl_dw[index_y][ii] -  trl_dw[index_y-1][ii])>PI)
			{
				len_up++;
				bad_up[ii]=true;
			}

		}

		if(len_dw || len_up)//if there are any jumps
		{
			quality_br_tmp=(float)((xres-len_up)*(xres-len_dw)/(xres*xres));
			if(quality_br_tmp<.95)
			{
				for(int ii=0; ii<xres; ii++)
				{
					if(bad_up[ii] && bad_dw[ii])//if it's bad on both sides
					{
						diff_up=trl_dw[index_y-1][ii] -  trl_dw[index_y][ii];
						diff_dw=trl_dw[index_y+1][ii] -  trl_dw[index_y][ii];
						phase_x[index_y+ii*yres]+=PI2*(round(0.5*(diff_up+diff_dw)/PI2));//split the difference
					}
				}
			}
		}
	}
	//Same as below with index_y = 0, so no check yres
	index_y=yres-1; len_dw=0;
	len_up=0;

	for(int ii=0;ii<xres; ii++)
	{
		bad_up[ii]=0;

		if(fabs(trl_dw[index_y][ii] -  trl_dw[index_y-1][ii])>PI)
		{
			len_up++;
			bad_up[ii]=1;
		}
	}
	if(len_up!=0)
	{
		quality_br_tmp=(xres-len_up)/(xres);

		if(quality_br_tmp<0.95)
		{
			for(int ii=0; ii<xres; ii++)
			{
				diff_up=phase_x[index_y-1+ii*yres] - phase_x[index_y+ii*yres];
				phase_x[index_y+ii*yres]+=PI2*(round(diff_up/PI2));
			}
		}
	}
////////////////////
	/*float* filter_out_x=filter(phase_x, xres, yres);
	float tmp_mean;
	int count;
	for(int ii=0; ii<yres*xres;ii++)
	{
		filter_out_x[ii]=phase_x[ii]-filter_out_x[ii];
		if(fabs(filter_out_x[ii])>=PI && md.support_MASK[ii]==1)
			phase_x[ii]-=round(filter_out_x[ii]/PI2)*PI2;//unwrapped_phase_x should be equivalent to out_x by now

	}

	tmp_mean=0;
	count=0;
	for(int ii=0; ii<yres*xres; ii++)
	{
			if(md.support_MASK[ii]==1)
			{
				tmp_mean+=phase_x[ii];
				count++;
			}
	}
	tmp_mean/=count;
	for(int ii=0; ii<yres*xres; ii++)
	{
			if(fullsignal || md.MASK[ii]==1)//should shortcircuit with fullsignal (i.e. not try to check mask)
				phase_x[ii]-=PI2*round(tmp_mean/PI2);
	}

	delete[] filter_out_x; //dyn allocated in filter
	*/
}

///Below this line, were part of PUROR.cpp not UnwrapGadget.cpp
void dePULM_1D(std::vector<float> &phase_1D, std::vector<int>& index_ls)
{
	//This function reconnects segments identified by dePULM_1D_ls
	float test, diff_1D;
	float inter_diff_l;
	int kk, length, tmp_1D;
	int k;
	std::vector<float> tmp_smooth(phase_1D);
	kk=0;
	std::vector<float> phase_1D_tmp(phase_1D);

	//for(int ii=0; ii<phase_1D.size(); ii++)
		//phase_1D_tmp[ii]=phase_1D[ii];

	if(index_ls.size()>2)//if at least one segment identified
	{
		for(int ii=1; ii<index_ls.size()-2; ii+=2)//ii is end of a segment, ii+1 is beginning of next
		{
			length=index_ls[ii+1]-index_ls[ii];
			if(length <=6)
			{
				test = phase_1D[index_ls[ii + 1]] + phase_1D[index_ls[ii + 1] + 1] + phase_1D[index_ls[ii + 1] + 2] - phase_1D[index_ls[ii]] - phase_1D[index_ls[ii] - 1] - phase_1D[index_ls[ii] - 2];
				test = test/3;
				//test=difference between the average of the first three points of next segment and last three points of this segment

				if(fabs(test)>PI)//represents a pole
				{
					test=round(test/(PI2))*PI2;
					kk=index_ls[ii]+2;
					for(k=kk; k<phase_1D.size(); k++)//shift phase data to eliminate discontinuity
					{
						phase_1D[k]-=test;
					}
				}
			}
			if(kk!=0)//has kk been changed(was there a pole)
			{
				inter_diff_l=phase_1D[kk]-phase_1D[kk-1];
				inter_diff_l=inter_diff_l>0?floor(inter_diff_l/PI2):ceil(inter_diff_l/PI2);//'fix' (biased round to zero)
				if(fabs(inter_diff_l)>0)//check for wrap in beginning of segment, fix it
					phase_1D[kk]-=PI2*round(inter_diff_l);

			}
		}

		tmp_1D=index_ls[0];//eliminate wraps at beginning of segment

		if(tmp_1D!=0)
		{
			if(tmp_1D==1)
			{
				diff_1D=phase_1D[tmp_1D]-phase_1D[0];
				phase_1D[0]=phase_1D_tmp[0] + PI2*round(diff_1D/PI2);
			}
			else
			{
				for(int ii=0; ii<tmp_1D; ii++)
					tmp_smooth[ii] = phase_1D_tmp[ii];

				diff_1D=phase_1D[tmp_1D]-mean(tmp_smooth.data(), tmp_1D);//awkward here, using part of a vector
				diff_1D= PI2*round(diff_1D/PI2);
				for(int ii=0; ii<tmp_1D; ii++)
					phase_1D[ii] = tmp_smooth[ii] + diff_1D;
			}
		}

		tmp_1D = index_ls[index_ls.size()-1];//eliminate wraps at end of segment

		if (tmp_1D != phase_1D.size()-1)
		{
			if (tmp_1D == phase_1D.size() - 2)
			{
				diff_1D = phase_1D[tmp_1D] - phase_1D_tmp[tmp_1D + 1];
				phase_1D[tmp_1D + 1] =phase_1D_tmp[tmp_1D+1]+PI2*round(diff_1D/PI2);
			}
			else
			{
				for(int ii=tmp_1D+1; ii<phase_1D.size(); ii++)
					tmp_smooth[ii-tmp_1D-1] = phase_1D_tmp[ii];

				diff_1D = phase_1D[tmp_1D] - mean(tmp_smooth.data(), phase_1D.size()-tmp_1D-1);
				diff_1D= PI2*round(diff_1D/PI2);

				for(int ii=tmp_1D+1; ii<phase_1D.size(); ii++)
					phase_1D[ii]=phase_1D_tmp[ii]+diff_1D;
			}
		}
	}

	return;
}

void dePULM_1D_brain(std::vector<float> &phase_1D, std::vector<int>& index_ls)
{
	//This function reconnects segments identified by dePULM_1D_ls
	float test, diff_1D;
	float inter_diff_l;
	int tmp_1D;
	std::vector<float> phase_1D_tmp(phase_1D);
	std::vector<float> tmp_smooth;
	int kk = 0;
	if(index_ls.size()>2)//if at least one segment identified
	{
		for(int ii=1; ii<index_ls.size()-2; ii+=2)//ii is end of a segment, ii+1 is beginning of next
		{
			if((index_ls[ii+1]-index_ls[ii]) <=6) //if there is less than seven between
			{
				test = phase_1D[index_ls[ii + 1]] + phase_1D[index_ls[ii + 1] + 1] + phase_1D[index_ls[ii + 1] + 2] - phase_1D[index_ls[ii]] - phase_1D[index_ls[ii] - 1] - phase_1D[index_ls[ii] - 2];
				test = test/3;
				//test=difference between the average of the first three points of next segment and last three points of this segment

				if(fabs(test)>PI)//represents a pole
				{
					test=round(test/(PI2))*PI2;
					kk=index_ls[ii]+2;
					for(int k=kk; k<phase_1D.size(); k++)//shift phase data to eliminate discontinuity
					{
						phase_1D[k]-=test;
					}
				}
			}
			if(kk!=0)//has kk been changed(was there a pole)
			{
				inter_diff_l=phase_1D[kk]-phase_1D[kk-1];
				inter_diff_l=inter_diff_l>0?floor(inter_diff_l/PI2):ceil(inter_diff_l/PI2);//'fix' (biased round to zero)
				if(fabs(inter_diff_l)>0)//check for wrap in beginning of segment, fix it
					phase_1D[kk]-=PI2*round(inter_diff_l);

			}
		}

		tmp_1D=index_ls[0];//eliminate wraps at beginning of segment
		tmp_smooth.resize(tmp_1D);
		if(tmp_1D!=0)
		{
			if(tmp_1D==1)
			{
				diff_1D=phase_1D[tmp_1D]-phase_1D[0];
				phase_1D[0]=phase_1D_tmp[0] + PI2*round(diff_1D/PI2);
			}
			else
			{
				for(int ii=0; ii<tmp_1D; ii++)
					tmp_smooth[ii] = phase_1D_tmp[ii];

				diff_1D=phase_1D[tmp_1D]-mean(tmp_smooth.data(), tmp_1D);
				diff_1D= PI2*round(diff_1D/PI2);
				for(int ii=0; ii<tmp_1D; ii++)
					phase_1D[ii] = tmp_smooth[ii] + diff_1D;
			}

		}

		tmp_1D = index_ls[index_ls.size()-1];//eliminate wraps at end of segment

		if (tmp_1D != phase_1D.size()-1)
		{
			if (tmp_1D == phase_1D.size() - 2)
			{
				diff_1D = phase_1D[tmp_1D] - phase_1D_tmp[tmp_1D + 1];
				phase_1D[tmp_1D + 1] =phase_1D_tmp[tmp_1D+1]+PI2*round(diff_1D/PI2);
			}
			else
			{
				tmp_smooth.resize(phase_1D.size()-tmp_1D-1);
				for(int ii=tmp_1D+1; ii<phase_1D.size(); ii++)
					tmp_smooth[ii-tmp_1D-1] = phase_1D_tmp[ii];

				diff_1D = phase_1D[tmp_1D] - mean(tmp_smooth.data(), phase_1D.size()-tmp_1D-1);
				diff_1D= PI2*round(diff_1D/PI2);

				for(int ii=tmp_1D+1; ii<phase_1D.size(); ii++)
					phase_1D[ii]=phase_1D_tmp[ii]+diff_1D;
			}

		}

	}
	return;
}

std::vector<int> dePULM_1D_ls(std::vector<float> &phase_1D,float phi_good, std::vector<int>& g_seg)
{
	//This function identifies segments (within mask segments) in phase_1D with phase gradient < phi_good
	int connections= 0;//kk in MATLAB
	bool end_flag = false;
	float diff_nb;
	std::vector<int> index_ls;

	if (phase_1D.empty())
	{
		index_ls.push_back(0);
		index_ls.push_back(-1);//added -1, probably never gets called but that would have been a bug
	}
	if(g_seg.size()>=2)
	{
		for(int ii=0; ii<g_seg.size()-1; ii+=2)//using the segments defined in ini from mask
		{
			connections=0; end_flag=0;
			for(int jj=g_seg[ii]; jj<g_seg[ii+1]; jj++)//from the start of the segment to the end of the segment
			{
				diff_nb=phase_1D[jj+1] - phase_1D[jj];
				if(end_flag==0)//if a segments with 4 pixels has not been found
				{
					if(fabs(diff_nb)<phi_good)//if the next point is within phi_good (PI/2), add one to length
						connections++;
					//else//if it isn't, reset the length counter
					//{
					//	connections=0;
					//}
				}

				if (connections==3)//if a segment with 4 pixels is found
				{
					index_ls.push_back(jj-2);//record start, reset length counter, set end flag
					connections = 0;
					end_flag = true;
				}

				if( end_flag && (diff_nb>phi_good || -diff_nb>phi_good))//gradient too high, save end point
				{
					index_ls.push_back(jj);
					end_flag = false;

				}
				if (end_flag && jj == (g_seg[ii+1] - 1))//end of mask segment reached, save end point
				{
					index_ls.push_back(jj+1);
				}
			}
		}
	}
	return index_ls;
}
std::vector<int> dePULM_1D_ls_brain(std::vector<float> &phase_1D, float phi_good)
{
	//This function identifies segmentsin phase_1D with phase gradient < phi_good
	int connections; //kk in MATLAB
	bool end_flag;
	float diff_nb;
	std::vector<int> index_ls;

	connections=0; end_flag=false;
	for(int jj=0; jj<phase_1D.size()-1; jj++)//from the start of the segment to the end of the segment
	{
		diff_nb=phase_1D[jj+1] - phase_1D[jj];
		if(!end_flag)//if a segments with 4 pixels has not been found
		{
			if(fabs(diff_nb)<phi_good)//if the next point is within phi_good (PI/2), add one to length
				connections++;
		}
		if (connections==3)//if a segment with 4 pixels is found
		{
			index_ls.push_back(jj-2);//record start, reset connections counter, set end flag
			connections = 0;
			end_flag = true;
		}
		if( end_flag && (diff_nb>phi_good || -diff_nb>phi_good))//gradient too high, save end point
		{
			index_ls.push_back(jj);
			end_flag = 0;

		}
		if (end_flag && jj == phase_1D.size()-2)//end of mask segment reached, save end point
		{
			index_ls.push_back(jj+1);
		}
	}

	return index_ls;
}
void dePULM_2D_merging(std::vector<float>& line_phase,std::vector<int>& index_ls, std::vector<float>& re_phase, int res, std::vector<float>& seg_mean_shift, std::vector<float>& diff_tmp)
{
	int count; //kk in MATLAB
	double sum_diff,ave_diff;
	int gphi, index_start, index_end; //in_s and in_e in MATLAB

	gphi=0;
	if(index_ls.size()>2)
	{
		for(int ii=0; ii<res; ii++)//find differences between phase_tmp and phase_tmp_ref from final
			diff_tmp[ii]=re_phase[ii]-line_phase[ii];

		for(int ii=0; ii<(index_ls.size()-1); ii+=2)
		{
			count=0;
			sum_diff=0;
			for(int jj=index_ls[ii]; jj<index_ls[ii+1]+1; jj++)
			{
				sum_diff+=diff_tmp[jj];
				count++;
			}
			if(count>0)
			{
				ave_diff= sum_diff/count;
				if(fabs(ave_diff)>PI)
					gphi=1;
				seg_mean_shift[ii]=ave_diff;//find the average difference for this segment
				seg_mean_shift[ii+1]=ave_diff;
			}
		}

		if(gphi!=0)
		{
			for(int ii=0; ii<index_ls.size(); ii+=2)
			{
				ave_diff=seg_mean_shift[ii];
				ave_diff=PI2*round(ave_diff/PI2);
				for(int jj=index_ls[ii]; jj<=index_ls[ii+1]; jj++)//adjust segment
				{
					line_phase[jj]+=ave_diff;
				}
			}

			for(int ii=1; ii<(index_ls.size()-2); ii+=2)
			{
				index_start=index_ls[ii]+1;
				index_end=index_ls[ii+1]-1;
				if((index_end-index_start)>=0)
				{
					count=0;
					sum_diff=0;
					for(int jj=index_ls[ii]-1; jj<=index_ls[ii+1]+1; jj++)
					{
						sum_diff+=diff_tmp[jj];
						count++;
					}

					if(count==0)
						ave_diff=(seg_mean_shift[ii]+seg_mean_shift[ii+1])/2;
					else
					{
						ave_diff=sum_diff/count;

					}
					ave_diff=PI2*round(ave_diff/PI2);
					for(int jj=index_start; jj<=index_end; jj++)//adjust points between segments
					{
						line_phase[jj]+=ave_diff;
					}
				}
			}
		}
	}

	return;
}
void dePULM_2D_itoh(std::vector<int>& mask, std::vector<float>& phase_tmp, int xres, int yres)
{
	/*unwrap phase_tmp*/
	std::vector<float> phase_1D_tmp;

	std::vector<int> index_m;

	for(int index_x=0; index_x<xres; index_x++)
	{
		for(int ii=0; ii<yres; ii++)
		{
			if(mask[ii+yres*index_x])
				index_m.push_back(ii);
		}
		int count = index_m.size();
		phase_1D_tmp.resize(count);
		if(count>1)
		{
			if(index_x%2==0)
			{
				for(int ii=0; ii<count; ii++)
					phase_1D_tmp[count-1-ii]=phase_tmp[index_m[ii]+yres*index_x];

				unwrap(phase_1D_tmp);

				for(int ii=0; ii<count; ii++)
					phase_tmp[index_m[ii]+yres*index_x]=phase_1D_tmp[count-1-ii];

			}
			else
			{
				for(int ii=0; ii<count; ii++)
					phase_1D_tmp[ii]=phase_tmp[index_m[ii]+yres*index_x];;

				unwrap(phase_1D_tmp);

				for(int ii=0; ii<count; ii++)
					phase_tmp[index_m[ii]+yres*index_x]=phase_1D_tmp[ii];
			}
		}

	}

	return;
}

void diff_y(std::vector<float>& phase_y,bool fullsignal, MaskData& md, int xres, int yres)
{
	int index_x;
	int len_dw,len_up;
	float diff_up, diff_dw, quality_br_tmp;

	std::vector<std::vector<float> > trl_dw(yres, std::vector<float>(xres));

	std::vector<bool> bad_up(yres);
	std::vector<bool> bad_dw(yres);

	for(int ii=0; ii<yres; ii++)
	{
		for(int jj=0; jj<xres; jj++)
			trl_dw[ii][jj]=phase_y[ii+jj*yres];
	}

	index_x=0; len_dw=0;
	//Same as below with index_x = 0, no -1 check
	for(int ii=0; ii<yres; ii++)
	{
		if(fabs(phase_y[ii] -phase_y[yres+ii])>PI)
			len_dw++;
	}

	if(len_dw!=0)
	{
		quality_br_tmp=(float)((yres-len_dw)/(yres));

		if(quality_br_tmp<0.95)
		{
			for(int ii=0; ii<yres; ii++)
			{
				diff_dw=phase_y[yres+ii] - phase_y[ii];
				phase_y[ii]+=PI2*(round(diff_dw/PI2));
			}
		}
	}

	for(index_x=1; index_x<(xres-1);index_x++)
	{
		len_up=0;
		len_dw=0;

		for(int ii=0;ii<yres; ii++)
		{
			bad_up[ii]=0;
			bad_dw[ii]=0;

			if(fabs(trl_dw[ii][index_x] -trl_dw[ii][index_x+1])>PI)
			{
				len_dw++;
				bad_dw[ii]=true;
			}

			if(fabs(trl_dw[ii][index_x] -trl_dw[ii][index_x-1])>PI)
			{
				len_up++;
				bad_up[ii]=true;
			}
		}

		if(len_dw || len_up)//if there are any jumps
		{
			quality_br_tmp=(float)((yres-len_up)*(yres-len_dw)/(yres*yres));
			if(quality_br_tmp<.95)
			{
				for(int ii=0; ii<yres; ii++)
				{
					if(bad_up[ii] && bad_dw[ii])
					{
						diff_up=trl_dw[ii][index_x-1] -trl_dw[ii][index_x];
						diff_dw=trl_dw[ii][index_x+1] -trl_dw[ii][index_x];
						phase_y[index_x*yres+ii]+=PI2*(round(0.5*(diff_up+diff_dw)/PI2));//split the difference
					}
				}
			}
		}
	}

	index_x=xres-1; len_dw=0;
	len_up=0;
	//Same as above with no xres check
	for(int ii=0;ii<yres; ii++)
	{
		bad_up[ii]=0;

		if(fabs(trl_dw[ii][index_x] -trl_dw[ii][index_x])>PI)
		{
			len_up++;
			bad_up[ii]=1;
		}
	}
	if(len_up!=0)
	{
		quality_br_tmp=(xres-len_up)/(xres);

		if(quality_br_tmp<0.95)
		{
			for(int ii=0; ii<yres; ii++)
			{
				diff_up=phase_y[(index_x-1)*yres+ii] - phase_y[index_x*yres+ii];
				phase_y[index_x*yres+ii]+=PI2*(round(diff_up/PI2));
			}
		}
	}
////////////////////
	/*float* filter_out_x=filter(phase_y, xres, yres);
	float tmp_mean;
	int count;
	for(int ii=0; ii<yres*xres;ii++)
	{
			filter_out_x[ii]=phase_y[ii]-filter_out_x[ii];
			if(fabs(filter_out_x[ii])>=PI && md.support_MASK[ii]==1)
				phase_y[ii]-=round(filter_out_x[ii]/PI2)*PI2;//unwrapped_phase_y should be equivalent to out_x by now

	}

	tmp_mean=0;
	count=0;
	for(int ii=0; ii<yres*xres; ii++)
	{
			if(md.support_MASK[ii]==1)
			{
				tmp_mean+=phase_y[ii];
				count++;
			}

	}
	tmp_mean/=count;
	for(int ii=0; ii<yres*xres; ii++)
	{
			if(fullsignal || md.MASK[ii]==1)//should shortcircuit with fullsignal (i.e. not try to check mask)
				phase_y[ii]-=PI2*round(tmp_mean/PI2);

	}*/

	//delete[] filter_out_x; //dyn allocated in filter
}
