#ifndef MASK_H
#define	MASK_H
#include "PUROR.h"
#include <vector>
class MaskData
{
public:
	~MaskData();
	MaskData();
	MaskData(int yres, int xres);

	void iniMask();
	void iniMask_y();
	bool fullsignal;//can speed up processing if every pixel is in the mask
	
		
int xres, yres;
int *MASK;
int *support_MASK;//MASK_h
int *support_MASK_trans;//transposed version of support_MASK

std::vector<std::vector<int>> signalY;
std::vector<std::vector<int>> connectYH;
std::vector<std::vector<int>> segY;

std::vector<std::vector<int>> signalX;
std::vector<std::vector<int>> connectXH;
std::vector<std::vector<int>> segX;
};
#endif //MASK_H
