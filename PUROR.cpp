#include "PUROR.h"



/*template<typename T>
T mean(T* sample, int length)
{
	T mean;
	mean=0;
	if(length==0)
		return 0;

	for(int ii=0; ii<length; ii++)
		mean+=sample[ii];
	mean/=length;
	
	return mean;
}*/

float mean(float *sample, int length)
{
	int ii;
	float mean;
	mean=0;
	for(ii=0; ii<length; ii++)
		mean+=sample[ii];
	mean/=length;
	if(length==0){
		mean=0;//Not correct for all applications
	}
	return mean;
}
double mean(double *sample, int length)
{
	int ii;
	double mean;
	mean=0;
	for(ii=0; ii<length; ii++)
		mean+=sample[ii];
	mean/=length;
	if(length==0){
		mean=0;//Not correct for all applications
	}
	return mean;
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

void unwrap(float* phase, int length)
{
	//float* diff, *dps, *dp_corr;
	int ii;
	std::vector<float> diff(length);
	std::vector<float> dps(length);
	std::vector<float> dp_corr(length);
	//diff= new float[length];
	//dps= new float[length];
	//dp_corr= new float[length];

	dp_corr[0]=0;

	for(ii=1; ii<length; ii++)		
		diff[ii]=phase[ii]-phase[ii-1]; //incremental phase step

	for(ii=1; ii<length; ii++)
	{
		dps[ii]=(float)remainder((diff[ii]+PI),PI2)-PI; // should be equivalent step in[-pi, pi], but negative values for first arg give negative results

		if(dps[ii]==-PI && diff[ii]>0)//check pi vs -pi
			dps[ii]=PI;

		if(dps[ii]<(-1*PI))//shouldn't happen, but does because of the way mod is done
			dps[ii]+=PI2;

		dp_corr[ii]=dps[ii]-diff[ii];//correction amount
		if(fabs(diff[ii])<PI)
			dp_corr[ii]=0;
		dp_corr[ii]+=dp_corr[ii-1];//integrate correction
		phase[ii]+=dp_corr[ii];//apply correction
	}
	//delete[] dp_corr;
	//delete[] diff;
	//delete[] dps;

	return;
}


void dePULM_1D(std::vector<float> &phase_1D, std::vector<int>& signal,std::vector<int>& index_ls)
{
	//This function reconnects segments identified by dePULM_1D_ls
	float test, inder_diff_l, diff_1D;
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

void dePULM_1D_brain(std::vector<float> &phase_1D, std::vector<int>& index_ls)
{
	//This function reconnects segments identified by dePULM_1D_ls
	float test, inder_diff_l, diff_1D;
	float inter_diff_l;
	int ii, kk, length, tmp_1D;
	int k;
	kk=0;
	std::vector<float> phase_1D_tmp(phase_1D);
	std::vector<float> tmp_smooth;

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
				for(ii=0; ii<tmp_1D; ii++)
					tmp_smooth[ii] = phase_1D_tmp[ii];

				diff_1D=phase_1D[tmp_1D]-mean(tmp_smooth.data(), tmp_1D);
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
				tmp_smooth.resize(phase_1D.size()-tmp_1D-1);
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


std::vector<int>  dePULM_1D_ls(std::vector<float> &phase_1D,float phi_good, std::vector<int>& g_seg)
{
	//This function identifies segments (within mask segments) in phase_1D with phase gradient < phi_good
	int ii,jj,kk, end_flag,count;
	float diff_nb;
	std::vector<int> index_ls;


	if (phase_1D.empty())
	{
		index_ls.push_back(0);
		index_ls.push_back(-1);//added -1, probably never gets called but that would have been a bug
	}
	if(g_seg.size()>=2)
	{

		for(ii=0; ii<g_seg.size()-1; ii+=2)//using the segments defined in ini from mask
		{	
			kk=0; end_flag=0;
			for(jj=g_seg[ii]; jj<g_seg[ii+1]; jj++)//from the start of the segment to the end of the segment
			{
				diff_nb=phase_1D[jj+1] - phase_1D[jj];
				if(end_flag==0)//if a segments with 4 pixels has not been found
				{
					if(fabs(diff_nb)<phi_good)//if the next point is within phi_good (PI/2), add one to length
						kk++;
					//else//if it isn't, reset the length counter
					//{
					//	kk=0;
					//}

				}

				if (kk==3)//if a segment with 4 pixels is found
				{			
					index_ls.push_back(jj-2);//record start, reset length counter, set end flag
					kk = 0;
					end_flag = 1;
				}

				
				if( end_flag != 0 && (diff_nb>phi_good || -diff_nb>phi_good))//gradient too high, save end point
				{
					index_ls.push_back(jj);
					end_flag = 0;

				}
				if (end_flag != 0 && jj == (g_seg[ii+1] - 1))//end of mask segment reached, save end point
				{
					index_ls.push_back(jj+1);
				}
			}
		}

	}
	return index_ls;
}
std::vector<int>   dePULM_1D_ls_brain(std::vector<float> &phase_1D, float phi_good)
{
	//This function identifies segmentsin phase_1D with phase gradient < phi_good
	int ii,jj,kk, end_flag;
	float diff_nb;
	std::vector<int> index_ls;

	
	kk=0; end_flag=0;
	for(jj=0; jj<phase_1D.size()-1; jj++)//from the start of the segment to the end of the segment
	{
		diff_nb=phase_1D[jj+1] - phase_1D[jj];
		if(end_flag==0)//if a segments with 4 pixels has not been found
		{
			if(fabs(diff_nb)<phi_good)//if the next point is within phi_good (PI/2), add one to length
				kk++;
		}

		if (kk==3)//if a segment with 4 pixels is found
		{			
			index_ls.push_back(jj-2);//record start, reset length counter, set end flag
			kk = 0;
			end_flag = 1;
		}

		
		if( end_flag != 0 && (diff_nb>phi_good || -diff_nb>phi_good))//gradient too high, save end point
		{
			index_ls.push_back(jj);
			end_flag = 0;

		}
		if (end_flag != 0 && jj == phase_1D.size()-2)//end of mask segment reached, save end point
		{
			index_ls.push_back(jj+1);
		}
	}
	
	
	return index_ls;
}
void dePULM_2D_merging(std::vector<float>& line_phase,std::vector<int>& index_ls, std::vector<float>& re_phase, int res, float* line_tmp_phase, float*seg_mean_shift, float*diff_tmp)
{
	int kk, jj, ii;
	double sum_diff,ave_diff;
	int gphi, in_s, in_e;

	//for(ii=0; ii<res; ii++)
		//line_tmp_phase[ii]=line_phase[ii];

	gphi=0;
	if(index_ls.size()>2)
	{
		for(ii=0; ii<res; ii++)//find differences between phase_tmp and phase_tmp_ref from final
			diff_tmp[ii]=re_phase[ii]-line_phase[ii];

		for(ii=0; ii<(index_ls.size()-1); ii+=2)
		{
			kk=0;
			sum_diff=0;
			for(jj=index_ls[ii]; jj<index_ls[ii+1]+1; jj++)
			{
				sum_diff+=diff_tmp[jj];
				kk++;
			}
			if(kk>0)
			{
				ave_diff= sum_diff/kk;
				if(fabs(ave_diff)>PI)
					gphi=1;
				seg_mean_shift[ii]=ave_diff;//find the average difference for this segment
				seg_mean_shift[ii+1]=ave_diff;
			}

		}

		if(gphi!=0)
		{
			for(ii=0; ii<index_ls.size(); ii+=2)
			{
				ave_diff=seg_mean_shift[ii];
				ave_diff=PI2*round(ave_diff/PI2);
				for(jj=index_ls[ii]; jj<=index_ls[ii+1]; jj++)//adjust segment 
				{
					line_phase[jj]+=ave_diff;
					//line_tmp_phase[jj]=line_phase[jj]+ave_diff;
				}
			}

			for(ii=1; ii<(index_ls.size()-2); ii+=2)
			{

				in_s=index_ls[ii]+1;
				in_e=index_ls[ii+1]-1;
				if((in_e-in_s)>=0)
				{

					kk=0;
					sum_diff=0;
					for(jj=index_ls[ii]-1; jj<=index_ls[ii+1]+1; jj++)
					{

						sum_diff+=diff_tmp[jj];
						kk++;
					}

					if(kk==0)
						ave_diff=(seg_mean_shift[ii]+seg_mean_shift[ii+1])/2;
					else
					{
						ave_diff=sum_diff/kk;

					}
					ave_diff=PI2*round(ave_diff/PI2);
					for(jj=in_s; jj<=in_e; jj++)//adjust points between segments
					{
						line_phase[jj]+=ave_diff;
						//line_tmp_phase[jj]=line_phase[jj]+ave_diff;
					}
				}
			}
		}

	}


	//for(ii=0; ii<res; ii++)
		//line_phase[ii]=line_tmp_phase[ii];

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

