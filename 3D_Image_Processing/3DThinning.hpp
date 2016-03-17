#ifndef __ml3DThinning_h
#define __ml3DThinning_h

#include <iostream>
#include <memory>
#include <vector>
#include <thread>
#include <math.h>
#include <algorithm>

#include "vtkSmartPointer.h"
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkMetaImageReader.h>
#include "vtkImageData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkUnsignedCharArray.h"
#include "vtkStringArray.h"
#include "vtkPointData.h"
#include "vtkDataSet.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkImageEuclideanDistance.h"
#include "vtkImageMarchingCubes.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkAppendPolyData.h"
#include "vtkImageGaussianSmooth.h"
#include "vtkLink.h"

#include "ml\util.hpp"
#include <boost/chrono.hpp>

#define VTK_NEW(CLASS_NAME) vtkSmartPointer<##CLASS_NAME>::New()

class ml3DThinning
{
public:
	ml3DThinning();
	~ml3DThinning();

	// Run to thinning, after setup filters
	void Update(); 
	void Delete();
	void AllFilterOn();
	void AllFilterOff();
	void SetInputFile(std::string _inputFile);
	void SetInputImageData(vtkImageData* _inputImageData);
	void SetVesselArray(unsigned int * _vesselArray, unsigned int _w, unsigned int _h, unsigned int _sn);
	void SetBinaryArray(unsigned char * _binaryArray, unsigned int _w, unsigned int _h, unsigned int _sn);

	// Median filter process
	void SetMedianFilterOn();
	void SetMedianFileterOff();
	// Default is OFF
	void SetMedianFileter(bool _isNeedMedian); 
	// Default is 3
	void SetMedianKernel(unsigned int _medianKernel); 

	// Linear smooth filter process
	void SetLinearSmoothFilterOn();
	void SetLinearSmoothFilterOff();
	// Default is OFF
	void SetLinearSmoothFilter(bool _isNeedLinearSmooth); 
	// Default is 26 (you can set 6, 18 or 26)
	void SetLinearSmoothKernel(unsigned int _linearSmoothKernel); 

	// Binarization process
	// Default is 1
	void SetBinarizationTheardshold(unsigned int _binarizationTheardshold); 

	// Nearest smooth filter process
	void SetNearestSmoothFilterOn();
	void SetNearestSmoothFilterOff();
	// Default is OFF
	void SetNearestSmoothFilter(bool _isNeedNearestSmooth); 
	// Default is 3
	void SetNearestSmoothKernel(unsigned int _nearestSmoothKernel); 

	// 3D closing 
	void SetClosingOn();
	void SetClosingOff();
	// Default is ON
	void SetClosing(bool _isNeedClosing); 
	// Default is 2
	void SetClosingKernel(unsigned int _closingKernel); 

	// 3D opening
	void SetOpeningOn();
	void SetOpeningOff();
	// Default is ON
	void SetOpening(bool _isNeedOpening); 
	// Default is 1
	void SetOpeningKernel(unsigned int _openingKernel); 

	// Second thinning process
	void SetSecondThinningOn();
	void SetSecondThinningOff();
	// Default is ON
	void SetSecondThinning(bool _isNeedSecondThinning); 
	// Default is 3
	void SetSecondThinningKernel(unsigned int _secondThinningKernel);

	// Dilation & Erosion
	void SetArraySize(unsigned int _w, unsigned int _h, unsigned int _sn);
	void ErosionOnArray(const int erosionSize, unsigned char * inputArray);
	void DilationOnArray(const int dilationSize, unsigned char * inputArray);

	// Second thinning process
	void SetDeleteErrorPaternOn();
	void SetDeleteErrorPaternOff();
	// Default is OFF
	void SetDeleteErrorPatern(bool _isNeedDeleteError); 
	// Default is 1
	void SetErrorPaternLenght(unsigned int _errorPaternLenght); 
	void SetErrorAlgorithm2(bool _needAlgorithm2);
	void SetRadiusThreshold(double _threshold);

	void DetectBifurcationFromArray(unsigned char *srcArray, unsigned char* dstArray, unsigned int _w, unsigned int _h, unsigned int _sn);
	void DetectErrorFromArray(unsigned char *binariArray, unsigned char* vesselArray, unsigned int _w, unsigned int _h, unsigned int _sn);
	void DeleteEdgesFromPoint(unsigned char *binariArray, unsigned char* vesselArray, unsigned int _w, unsigned int _h, unsigned int _sn, vector3d startP);
	void DeletePointsInArray(unsigned char *binariArray, unsigned char* vesselArray, unsigned int _w, unsigned int _h, unsigned int _sn, points3d deletePoints);

	points3d GetBifurcationPoints();
	vtkSmartPointer<vtkPolyData> GetThinnedDataWithSpacing();
	vtkSmartPointer<vtkPolyData> GetThinnedDataWithoutSpacing();

	void GetDimensions(unsigned int &_w, unsigned int &_h, unsigned int &_sn);
	void GetSpacing(double &_xScale, double &_yScale, double &_zScale);
	unsigned int GetMaxVoxelValue();
	unsigned char * GetThinBinaryArray(); // Return Array2
	unsigned char * GetThinVessselArray(); // Return Array3

	void DeleteOldArrays();
	unsigned int  *inputArray, *outputArray;
	unsigned char *vesselArray, *bifurcationArray;

private:
	points3d bifurcationPoints;
	vtkSmartPointer<vtkPolyData> thinnedDataWithSpacing, thinnedDataWithoutSpacing;
	unsigned int maxVoxelValue;

	void StartMedian();
	void MedianBetween(unsigned int fromIdx, unsigned int toIdx);

	void StartLinearSmooth();
	void LinearSmoothBetween(unsigned int fromIdx, unsigned int toIdx);

	void StartBinarization();
	void BinarizationBetween(unsigned int fromIdx, unsigned int toIdx);

	void StartNearestSmooth();
	void NearestSmoothBetween(unsigned int fromIdx, unsigned int toIdx);

	// Closing process : dilated then eroded
	void StartClosing(); 
	// Openning process : eroded then dilated
	void StartOpenning(); 
	// Delation process
	void DilationBetween(const int _n, const int fromIdx, const int toIdx, const int type); 
	// Erosion process
	void ErosionBetween(const int _n, const int fromIdx, const int toIdx, const int type); 
	// Pararell join process array to vessel array
	void Join2VesselArrayBetween(const int fromIdx, const int toIdx, const int type); 
	// Pararell thinning
	void StartThinning(); 
	void ThinningBetween(const int fromSlice, const int toSlice, const int type);
	int Deletable(int _x, int _y, int _z, unsigned char *Array2ForParalell, int *ArrayLabelling, const int type);
	int LiLabeling3D(const int _sizex, const int _sizey, const int _sizez, const int *pBin, int *pObj);
	int LiLabeling3D_2(const int _sizex, const int _sizey, const int _sizez, const int* pBin, int* pObj);
	// Dectect bifurcation and end points
	void DetectionBifurcationPoints(); 
	void ReDetectBifurcation();
	int neigborNumberFrom(vector3d p);

	void ErrorDeleteProcess2();
	void CheckTrullyError();

	void DirectionVectorFromPoint(const double startP[3], vtkSmartPointer<vtkPoints> pointList, vector3d &direction);
	void DirectionVectorFromPoint2(const double startP[3], vtkSmartPointer<vtkPoints> pointList, vtkSmartPointer<vtkPoints> detectedPointList, vector3d &direction);
	void DirectionVectorFromPoint3(const double startP[3], vtkSmartPointer<vtkPoints> pointList, vtkSmartPointer<vtkPoints> detectedPointList, vector3d &foot, vector3d &direction);

	double RadiusWithCenter(const double center[3]);
	double RadiusWithCenter2(const double center[3], vtkSmartPointer<vtkPoints> detectedPointList);
	double RadiusWithCenter3(const double center[3], vtkSmartPointer<vtkPoints> detectedPointList);

	// The error patern delete process
	// Error means that error network like a short branch
	void ErrorDeleteProcess(); 
	void DeleteErrorFromBifurcation();
	void DetectError(const int aPoint[3], bool &av);
	void DetectErrorFromBifurcation(const int aPoint[3], int &lenght);
	bool IsValidPoint(const int x, const int y, const int z);

	// Functions support DeleteEdgesFromPoint-Function
	void NeighborPointsFromPoint(const vector3d startPoint, points3d &foundedPoints, points3d &edgePoints);
	void FindNeighborPointsFromPoint(const vector3d fromPoint, points3d &edgePoints);
};

#endif