#include "Mask.h"
#include <iostream>
MaskData::MaskData(){}
MaskData::~MaskData()
{
		delete[] MASK;
		delete[] support_MASK;
		delete[] support_MASK_trans;
}
MaskData::MaskData(int yres, int xres)
{

	this->yres=yres;
	this->xres=xres;
	MASK = new int[yres*xres];
	support_MASK = new int[yres*xres];
	support_MASK_trans = new int[yres*xres];
	signalX.resize(yres);
	connectXH.resize(yres);
	segX.resize(yres);
	fullsignal=0;
	signalY.resize(xres);
	connectYH.resize(xres);
	segY.resize(xres);
}
void  MaskData::iniMask()
{
	int end_flag,ii,jj,kk;
	std::vector<int> sum_tmp(xres);
	if(!fullsignal)
	{
		for(ii=0; ii<yres;ii++)
		{
			for(jj=0; jj<xres;jj++)
			{
				if(MASK[ii+yres*jj])//transfer list of points to signal array
					signalX[ii].push_back(jj);
			}
			if(signalX[ii].size()>6)
			{
				kk=0; end_flag=0;
				for(jj=0; jj<signalX[ii].size()-1; jj++)//identify start and end points of line segments in rows
				{
					if(signalX[ii][jj+1]-signalX[ii][jj]==1)//if next point is also in mask:
					{
						if(end_flag==0)//if the identified segment is <=3 points long add one to the length
							kk++;	
						if(kk==3)//if the identified segment is 4 pixels long (three connections)
						{
							segX[ii].push_back(jj-2); //save the beginning point
							end_flag=1;	
							kk=0;	//reset the length counter
						}
						if(jj==(signalX[ii].size()-2) && end_flag==1)//end of mask reached for this row, segment ends here if it's long enough
						{
							segX[ii].push_back(jj+1);//save end point
							end_flag=0;
							kk=0;
						}

					}
					else //if next point is not in mask
					{
						if(end_flag==1)//if identified segment is 4+ pixels long
						{
							segX[ii].push_back(jj); //save the end point
							end_flag=0;
						}
					}
				}
			}
		}			
	}
	for(ii=0; ii<yres;ii++)
	{
		if(ii!=0 && ii!=yres-1)
		{
			for(jj=0; jj<xres; jj++)
				sum_tmp[jj]=support_MASK[(ii-1)+yres*jj] + support_MASK[ii+yres*jj] + support_MASK[(ii+1)+yres*jj];
			
			for(jj=1; jj<xres-1; jj++)	
			{
				if(sum_tmp[jj]==3 && ((support_MASK[ii+yres*(jj-1)] + support_MASK[ii+yres*(jj+1)])==2))//identify support_MASK points with points on four sides
					connectXH[ii].push_back(jj);//move index up and assign
			}
			if(connectXH[ii].size()<=xres/4-1)//if not enough are indentified, use centers of 1x3 segments 
			{
				//connectXH.clear();
				for(jj=1; jj<xres; jj++)
				{
					if(sum_tmp[jj]==3)
						connectXH[ii].push_back(jj);
				}
			}

		}
		else
		{
			if(ii==0)
				for(jj=0; jj<xres; jj++)
					sum_tmp[jj]=support_MASK[(ii+1)+yres*jj] + support_MASK[ii+yres*jj];
			else//ii==yres-1
				for(jj=0; jj<xres; jj++)
					sum_tmp[jj]=support_MASK[(ii-1)+yres*jj] + support_MASK[ii+yres*jj];
			for(jj=1; jj<xres-1; jj++)	
			{
				if(sum_tmp[jj]==2 && ((support_MASK[ii+yres*(jj-1)] + support_MASK[ii+yres*(jj+1)])==2))//identify mask points with points on three sides (on edge)
					connectXH[ii].push_back(jj);
			}
			if(connectXH[ii].size()<=xres/4-1)//if not enough are indentified, use centers of 1x2 segments 
			{
				//connectXH.clear();
				for(jj=1; jj<xres; jj++)
				{
					if(sum_tmp[jj]==2)
						connectXH[ii].push_back(jj);
				}
			}
		}	
	}
	return;
}
void MaskData::iniMask_y()
{
	int kk, end_flag,ii,jj;

	for(ii=0; ii<xres;ii++)
		for(jj=0; jj<yres;jj++)
			support_MASK_trans[jj*xres+ii]=support_MASK[jj+yres*ii];

	std::vector<int> sum_tmp(yres);
	if(!fullsignal)
	{

		for(ii=0; ii<xres;ii++)
		{
			for(jj=0; jj<yres;jj++)
			{
				if(MASK[jj+yres*ii])//transfer list of points to signal array
					signalY[ii].push_back(jj);
			}		

			if(signalY[ii].size()>6)
			{
				kk=0; end_flag=0;
				for(jj=0; jj<signalY[ii].size()-1; jj++)//identify start and end poiints of line segmetns in columns
				{
					if(signalY[ii][jj+1]-signalY[ii][jj]==1)//if next point also in mask
					{
						if(end_flag==0)//if the identified segment is <=3 points long add one to the length
							kk++;
						if(kk==3)//if the identified segment is 4 pixels long (three connections)
						{
							segY[ii].push_back(jj-2);//save beginning point
							end_flag=1;
							kk=0;//reset length counter
						}
						if(jj==(signalY[ii].size()-2) && end_flag==1)//end of mask reached for this column, segment ends here if it's long enough
						{
							segY[ii].push_back(jj+1);//save end point
							end_flag=0;
							kk=0;
						}
					}
					else //is next point is not in mask
					{
						if(end_flag==1)//if identified segment is 4+ pixels long

						{
							segY[ii].push_back(jj);//save end point
							end_flag=0;
						}
					}
				}
			}
		}
	}
	for(ii=0; ii<xres;ii++)
	{
		if(ii!=0 && ii!=xres-1)
		{
			for(jj=0; jj<yres; jj++)
				sum_tmp[jj]=support_MASK[jj+(ii-1)*yres]+support_MASK[jj+(ii*yres)]+support_MASK[jj+(ii+1)*yres];
			for(jj=1; jj<yres-1; jj++)	
			{
				if(sum_tmp[jj]==3 && ((support_MASK[(jj-1)+ii*yres]+support_MASK[(jj+1)+ii*yres])==2))//identify mask points with points on four sides
					connectYH[ii].push_back(jj);//move index up and assign
			}
			if(connectYH.size()<=yres/4-1)//if not enough found, use centers of 3x1 segments
			{
				//connectYH.clear();
				for(jj=1; jj<yres; jj++)
				{
					if(sum_tmp[jj]==3)
						connectYH[ii].push_back(jj);
				}
			}
		}
		else
		{
			if(ii==0)
				for(jj=0; jj<yres; jj++)
					sum_tmp[jj]=support_MASK[jj+(ii+1)*yres]+support_MASK[jj+ii*yres];
			else//ii==xres-1
				for(jj=0; jj<yres; jj++)
					sum_tmp[jj]=support_MASK[jj+(ii-1)*yres]+support_MASK[jj+ii*yres];

			for(jj=1; jj<yres-1; jj++)	
			{
				if(sum_tmp[jj]==2 && ((support_MASK[(jj-1)+ii*yres]+support_MASK[(jj+1)+ii*yres])==2))//identify mask points with points on three sides (on edge)
					connectYH[ii].push_back(jj);
			}
			if(connectYH.size()<=yres/4-1)//if not enough are indentified, use centers of 2x1 segments 
			{
				//connectYH.clear();
				for(jj=1; jj<yres; jj++)
				{
					if(sum_tmp[jj]==2)
						connectYH[ii].push_back(jj);
				}
			}
		}	

	}
return;
}

