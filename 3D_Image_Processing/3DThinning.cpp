#include "3DThinning.hpp"


/*2014-9-15 : Params for test new error pattern processing*/
class Error
{
public:
	vtkSmartPointer<vtkPoints> pointList;
	vector<double> radiusList;
	double averageRadius;
	bool needDelete;

	Error(vtkSmartPointer<vtkPoints> _pointList) : pointList(_pointList){ needDelete = false; radiusList.clear(); averageRadius = 0; }
	~Error(){};
};
vector<Error> errorList;
vtkSmartPointer<vtkPoints> currentErrorPointList;
vtkSmartPointer<vtkModifiedBSPTree> bspTree;
vtkSmartPointer<vtkAppendPolyData> appendPoly;
/*2014-9-?? : End Params for test new error pattern processing*/




#pragma warning( disable: 4018 )

const int dir6[] = { 0, 0, -1, 0, -1, 0, -1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1 };
const int dir18[] = { 0, 0, -1, 0, -1, 0, -1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, -1, 1, 0, 1, -1, 0, -1, -1, 0,
0, 1, 1, 0, 1, -1, 0, -1, -1, 0, -1, 1, 1, 0, 1, 1, 0, -1, -1, 0, 1, -1, 0, -1 };
const int dir26[] = { 0, 0, -1, 0, -1, 0, -1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, -1, 1, 0, 1, -1, 0, -1, -1, 0,
0, 1, 1, 0, 1, -1, 0, -1, -1, 0, -1, 1, 1, 0, 1, 1, 0, -1, -1, 0, 1, -1, 0, -1, 1, 1, -1, 1, -1, -1,
-1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, 1 };

string volumePath;
//unsigned int  *inputArray, *outputArray;
//unsigned char *vesselArray, *bifurcationArray, *array21, *array22, *array23, *array24, *charOutputArray;
unsigned char *array21, *array22, *array23, *array24, *charOutputArray;

double spacing[3], bound[6], radiusThreshold; int ext[6];
unsigned int width, height, slice_num, size2D, size3D, medianKernel, linearSmoothKernel, binarizationTheardshold, nearestSmoothKernel, closingKernel, openingKernel, secondThinningKernel, errorPaternLenght;
bool isNeedMedian, isNeedLinearSmooth, isNeedNearestSmooth, isNeedClosing, isNeedOpening, isNeedSecondThinning, isNeedDeleteError, isNeedErrorAlgorithm2;

vtkImageData* convertedImg;
vtkSmartPointer<vtkPoints> oneEdgePoints4Detect, detectedPoints;

// Add on 2014-9-15
bool IsSameError(Error sourceError, Error targetError)
{
	if (sourceError.pointList->GetNumberOfPoints() != targetError.pointList->GetNumberOfPoints()){
		return false;
	}

	for (vtkIdType id = 0; id < targetError.pointList->GetNumberOfPoints(); id++){
		const double *targetPoint = targetError.pointList->GetPoint(id);
		if (!ExistedPointInSet(targetPoint, sourceError.pointList)){
			return false;
		}
	}

	return true;
}
bool IsAvaibleErrorInArray(Error targetError, vector<Error> _errorList)
{
	if (_errorList.size() <= 0){
		return false;
	}

	for (size_t index = 0; index < _errorList.size(); index++){
		Error sourceError = _errorList[index];
		if (!IsSameError(sourceError, targetError)){
			return false;
		}
	}

	return true;
}
// End

struct List{
	int i;
	int j;
	int k;
	int l;
};

ml3DThinning::ml3DThinning()
{
	width = height = slice_num = 0;
	size2D = size3D = 0;
	medianKernel = 3;
	linearSmoothKernel = 26;
	binarizationTheardshold = 1;
	nearestSmoothKernel = 3;
	closingKernel = 2;
	openingKernel = 1;
	secondThinningKernel = 3;
	errorPaternLenght = 7;
	isNeedMedian = isNeedLinearSmooth = false;
	isNeedNearestSmooth = isNeedDeleteError = isNeedClosing = isNeedOpening = isNeedSecondThinning = true;
	memset(spacing, 1, sizeof(double)* 3);

	// Add on 2014-9-24
	radiusThreshold = 0.4;
	isNeedErrorAlgorithm2 = false;
}

ml3DThinning::~ml3DThinning()
{
	printf("<ml3DThinning>Destructor was called\n");
}

void ml3DThinning::DeleteOldArrays()
{
	/*if (inputArray!=0)
	delete[] inputArray;
	if (outputArray!=0)
	delete[] outputArray;*/
	if (vesselArray != 0)
		delete[] vesselArray;
	if (bifurcationArray != 0)
		delete[] bifurcationArray;
}

void ml3DThinning::AllFilterOn()
{
	isNeedMedian = isNeedLinearSmooth = isNeedNearestSmooth = isNeedDeleteError = isNeedClosing = isNeedOpening = isNeedSecondThinning = true;
}

void ml3DThinning::AllFilterOff()
{
	isNeedMedian = isNeedLinearSmooth = isNeedNearestSmooth = isNeedDeleteError = isNeedClosing = isNeedOpening = isNeedSecondThinning = false;
}

// Add on 2014-9-24
void ml3DThinning::SetErrorAlgorithm2(bool _needAlgorithm2)
{
	isNeedErrorAlgorithm2 = _needAlgorithm2;
}

void ml3DThinning::SetRadiusThreshold(double _threshold)
{
	radiusThreshold = _threshold;
}
// End add

void ml3DThinning::Update()
{
	if (width == 0 && height == 0 && slice_num == 0)
		return;

	int TotalTime = 0;

	//1
	std::cout << "median processing... ";
	auto start = boost::chrono::steady_clock::now();
	if (isNeedMedian){
		ml3DThinning::StartMedian();
	}
	auto time = boost::chrono::steady_clock::now() - start; TotalTime = boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count();
	std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count() << " [msec]" << std::endl;

	//2
	std::cout << "linear processing... ";
	start = boost::chrono::steady_clock::now();
	if (isNeedLinearSmooth){
		ml3DThinning::StartLinearSmooth();
	}
	time = boost::chrono::steady_clock::now() - start; TotalTime += boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count();
	std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count() << " [msec]" << std::endl;

	//3
	std::cout << "binarization processing... ";
	start = boost::chrono::steady_clock::now();
	StartBinarization(); delete[] inputArray;
	time = boost::chrono::steady_clock::now() - start; TotalTime += boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count();
	std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count() << " [msec]" << std::endl;

	//4
	std::cout << "mean processing... ";
	start = boost::chrono::steady_clock::now();
	if (isNeedNearestSmooth){
		ml3DThinning::StartNearestSmooth();
	}
	time = boost::chrono::steady_clock::now() - start; TotalTime += boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count();
	std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count() << " [msec]" << std::endl;

	//5
	std::cout << "closing processing... ";
	start = boost::chrono::steady_clock::now();
	if (isNeedClosing){
		ml3DThinning::StartClosing();
	}
	time = boost::chrono::steady_clock::now() - start; TotalTime += boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count();
	std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count() << " [msec]" << std::endl;

	//6
	std::cout << "opening processing... ";
	start = boost::chrono::steady_clock::now();
	if (isNeedOpening)
	{
		ml3DThinning::StartOpenning();
	}
	time = boost::chrono::steady_clock::now() - start; TotalTime += boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count();
	std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count() << " [msec]" << std::endl;

	//7
	std::cout << "thinning processing... ";
	start = boost::chrono::steady_clock::now();
	ml3DThinning::StartThinning();
	time = boost::chrono::steady_clock::now() - start; TotalTime += boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count();
	std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count() << " [msec]" << std::endl;

	//8
	if (isNeedSecondThinning){
		std::cout << "thinning processing... ";
		start = boost::chrono::steady_clock::now();

		//ml3DThinning::DilationOnArray( secondThinningKernel, vesselArray );
		//ml3DThinning::ErosionOnArray( secondThinningKernel, vesselArray );
		closingKernel = secondThinningKernel;
		ml3DThinning::StartClosing();
		ml3DThinning::StartThinning();

		time = boost::chrono::steady_clock::now() - start; TotalTime += boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count();
		std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count() << " [msec]" << std::endl;
	}

	//9
	std::cout << "detection processing... ";
	start = boost::chrono::steady_clock::now();
	ml3DThinning::DetectionBifurcationPoints();
	time = boost::chrono::steady_clock::now() - start; TotalTime += boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count();
	std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count() << " [msec]" << std::endl;


	//10
	std::cout << "error pattern processing... ";
	start = boost::chrono::steady_clock::now();
	if (isNeedDeleteError){
		if (!isNeedErrorAlgorithm2){
			ml3DThinning::ErrorDeleteProcess();
		}
		else{
			errorList.clear();
			appendPoly = vtkSmartPointer<vtkAppendPolyData>::New();
			ml3DThinning::ErrorDeleteProcess2();
		}
	}
	time = boost::chrono::steady_clock::now() - start; TotalTime += boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count();
	std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count() << " [msec]" << std::endl;

	std::cout << "little error pattern processing... ";
	start = boost::chrono::steady_clock::now();
	ml3DThinning::ReDetectBifurcation();
	time = boost::chrono::steady_clock::now() - start; TotalTime += boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count();
	std::cout << boost::chrono::duration_cast<boost::chrono::milliseconds>(time).count() << " [msec]" << std::endl;

	std::cout << "TOTAL TIME = " << TotalTime << std::endl;
}

void ml3DThinning::SetInputFile(std::string _inputFile)
{
	volumePath = _inputFile;

	vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(_inputFile.c_str());
	reader->Update();

	ml3DThinning::SetInputImageData(reader->GetOutput());




	// No need to smooth using vtkImageGaussian, now must to smooth using itkAntiAlias
	//const double *scalar_range = reader->GetOutput()->GetScalarRange();
	//if (scalar_range[1] <= 1){
	//	ml3DThinning::SetInputImageData(reader->GetOutput());
	//}
	//else {
	//	vtkSmartPointer<vtkImageGaussianSmooth> smoothFilter = vtkSmartPointer<vtkImageGaussianSmooth>::New();
	//	smoothFilter->SetInputData(reader->GetOutput());
	//	smoothFilter->SetStandardDeviations(1, 1, 1);// (5, 5, 5);
	//	smoothFilter->SetRadiusFactors(0.5, 0.5, 0.5);
	//	smoothFilter->Update();

	//	ml3DThinning::SetInputImageData(smoothFilter->GetOutput());
	//}
}

void ml3DThinning::SetInputImageData(vtkImageData* _inputImageData)
{
	// Get volume dimensions, spacing, bound, extent ....
	int dims[3];
	_inputImageData->GetDimensions(dims);
	_inputImageData->GetExtent(ext);
	_inputImageData->GetBounds(bound);

	width = dims[0];
	height = dims[1];
	slice_num = dims[2];
	size2D = dims[0] * dims[1];
	size3D = dims[0] * dims[1] * dims[2];

	_inputImageData->GetSpacing(spacing);

	this->maxVoxelValue = 0;
	inputArray = new unsigned int[size3D];
	const double *scalar_range = _inputImageData->GetScalarRange();

	// Extract volume surface for new error processing algorithm
	if (isNeedDeleteError && isNeedErrorAlgorithm2){
		// Extract volume surface
		vtkSmartPointer<vtkImageMarchingCubes> marching_cubes = vtkSmartPointer<vtkImageMarchingCubes>::New();
		marching_cubes->ComputeScalarsOff();
		marching_cubes->ComputeNormalsOn();
		marching_cubes->SetValue(0, scalar_range[1] <= 1 ? 0.5 : binarizationTheardshold);
		//marching_cubes->SetValue(0, binarizationTheardshold);
		marching_cubes->SetInputData(_inputImageData);
		marching_cubes->Update();

		/*vtkSmartPointer<vtkXMLPolyDataWriter> wrTest = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		wrTest->SetInputConnection(marching_cubes->GetOutputPort());
		wrTest->SetFileName("surfaceTest.vtp");
		wrTest->Write();*/

		/*vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		smoothFilter->SetInputConnection(marching_cubes->GetOutputPort());
		smoothFilter->SetNumberOfIterations(500);
		smoothFilter->Update();*/

		vtkSmartPointer<vtkDataSetSurfaceFilter> extractSurfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
		extractSurfaceFilter->SetInputConnection(marching_cubes->GetOutputPort());
		extractSurfaceFilter->Update();

		bspTree = vtkSmartPointer<vtkModifiedBSPTree>::New();
		bspTree->SetDataSet(extractSurfaceFilter->GetOutput());
		bspTree->BuildLocator();
	}
	// End extract volume surface

	// The input volume was smooth by itkAntiAlias, so scalar type is UNSIGNED_CHAR
	for (int z = ext[4]; z < ext[5]; ++z){
		for (int y = ext[2]; y < ext[3]; ++y){
			for (int x = ext[0]; x < ext[1]; ++x){

				unsigned char* voxel = static_cast<unsigned char*>(_inputImageData->GetScalarPointer(x, y, z));
				if ((unsigned int)voxel[0] != 0){
					this->maxVoxelValue = this->maxVoxelValue >(unsigned int)voxel[0] ? this->maxVoxelValue : (unsigned int)voxel[0];
				}
				int _x = x - ext[0];
				int _y = y - ext[2];
				int _z = z - ext[4];

				inputArray[_x + _y*width + _z*size2D] = (unsigned int)voxel[0];
			}
		}
	}



	/*for (int z = ext[4]; z < ext[5]; ++z){
	for (int y = ext[2]; y < ext[3]; ++y){
	for (int x = ext[0]; x < ext[1]; ++x){

	if (_inputImageData->GetScalarType() == VTK_UNSIGNED_CHAR){

	unsigned char* voxel = static_cast<unsigned char*>(_inputImageData->GetScalarPointer(x, y, z));
	if ((unsigned int)voxel[0] != 0){
	this->maxVoxelValue = this->maxVoxelValue >(unsigned int)voxel[0] ? this->maxVoxelValue : (unsigned int)voxel[0];
	}
	int _x = x - ext[0];
	int _y = y - ext[2];
	int _z = z - ext[4];

	inputArray[_x + _y*width + _z*size2D] = (unsigned int)voxel[0];
	}
	else {
	if (scalar_range[1] <= 1){
	double* voxel = static_cast<double*>(_inputImageData->GetScalarPointer(x, y, z));
	if (*voxel > 0){
	int _x = x - ext[0];
	int _y = y - ext[2];
	int _z = z - ext[4];
	inputArray[_x + _y*width + _z*size2D] = 100;
	}
	}
	else {
	unsigned short* voxel = static_cast<unsigned short*>(_inputImageData->GetScalarPointer(x, y, z));
	if (*voxel > 0){
	unsigned int unsigned_char_value = 255 * (*voxel) / scalar_range[1];
	this->maxVoxelValue = this->maxVoxelValue > unsigned_char_value ? this->maxVoxelValue : unsigned_char_value;

	int _x = x - ext[0];
	int _y = y - ext[2];
	int _z = z - ext[4];
	inputArray[_x + _y*width + _z*size2D] = unsigned_char_value;
	}
	}
	}
	}
	}
	}

	cout << "Number of scalar components : " << _inputImageData->GetNumberOfScalarComponents() << endl;
	cout << "Max voxel value : " << this->maxVoxelValue << endl;*/

	/*vtkSmartPointer<vtkImageData> testImage = vtkSmartPointer<vtkImageData>::New();
	testImage->SetSpacing(spacing);
	testImage->SetExtent(ext);
	testImage->AllocateScalars(VTK_UNSIGNED_SHORT, 1);

	for (int z = ext[4]; z < ext[5]; ++z){
	for (int y = ext[2]; y < ext[3]; ++y){
	for (int x = ext[0]; x < ext[1]; ++x){
	unsigned short* voxel = static_cast<unsigned short*>(testImage->GetScalarPointer(x, y, z));
	int _x = x - ext[0];
	int _y = y - ext[2];
	int _z = z - ext[4];

	if (inputArray[_x + _y*width + _z*size2D] != 0){
	*voxel = 255;
	}
	}
	}
	}

	testImage->Modified();


	vtkSmartPointer<vtkXMLImageDataWriter> wr = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	wr->SetInputData(_inputImageData);
	wr->SetFileName("testOut.vti");
	wr->Write();*/
}

void ml3DThinning::SetVesselArray(unsigned int * _vesselArray, unsigned int _w, unsigned int _h, unsigned int _sn)
{
	width = _w;
	height = _h;
	slice_num = _sn;
	size2D = _w*_h;
	size3D = _w*_h*_sn;
	inputArray = _vesselArray;
}

void ml3DThinning::SetBinaryArray(unsigned char * _binaryArray, unsigned int _w, unsigned int _h, unsigned int _sn)
{
	width = _w;
	height = _h;
	slice_num = _sn;
	size2D = _w*_h;
	size3D = _w*_h*_sn;
	vesselArray = _binaryArray;
}

// Median filter process
void ml3DThinning::SetMedianFilterOn()
{
	isNeedMedian = true;
}

void ml3DThinning::SetMedianFileterOff()
{
	isNeedMedian = false;
}

void ml3DThinning::SetMedianFileter(bool _isNeedMedian)
{
	isNeedMedian = _isNeedMedian;
}

void ml3DThinning::SetMedianKernel(unsigned int _medianKernel)
{
	if (_medianKernel % 2 == 0){
		if (_medianKernel >= 4)
			medianKernel = _medianKernel - 1;
	}
	else {
		if (_medianKernel >= 3)
			medianKernel = _medianKernel;
	}
}

void ml3DThinning::StartMedian()
{
	if (medianKernel == 0)
		return;
	if (width == 0 && height == 0 && slice_num == 0)
		return;

	outputArray = new unsigned int[size3D];
	memset(outputArray, 0, sizeof(unsigned int)*size3D);

	std::vector<std::thread> parallelMedian;
	parallelMedian.resize(4);

	parallelMedian[0] = std::thread([this] {MedianBetween((medianKernel - 1) / 2, slice_num / 4 + 10); });
	parallelMedian[1] = std::thread([this] {MedianBetween(slice_num / 4 - 10, slice_num / 2 + 10); });
	parallelMedian[2] = std::thread([this] {MedianBetween(slice_num / 2 - 10, 3 * slice_num / 4 + 10); });
	parallelMedian[3] = std::thread([this] {MedianBetween(3 * slice_num / 4 - 10, slice_num - (medianKernel - 1) / 2); });

	for (auto it = parallelMedian.begin(); it != parallelMedian.end(); it++)
		it->join();

	parallelMedian.clear();

	memset(inputArray, 0, sizeof(unsigned int)*size3D);
	memcpy(inputArray, outputArray, sizeof(unsigned int)*size3D);
	delete[] outputArray;
}

void ml3DThinning::MedianBetween(unsigned int fromIdx, unsigned int toIdx)
{
	int kernelX3 = pow(medianKernel, 3);
	//int kernelX2 = pow( medianKernel, 2 ); //2D Median

	//int *xxxDATA = (int *)calloc( kernelX2, sizeof(int));
	int *xxxDATA = new int[kernelX3];
	memset(xxxDATA, 0, sizeof(int)*kernelX3);

	int distance = ((medianKernel - 1) / 2);

	for (unsigned int z = fromIdx; z < toIdx; z++){
		for (unsigned int y = distance; y < height - distance; y++){
			for (unsigned int x = distance; x < width - distance; x++){

				int index = 0;

				for (int p = -distance; p <= distance; p++){
					for (int q = -distance; q <= distance; q++){
						for (int r = -distance; r <= distance; r++){
							xxxDATA[index] = inputArray[(x + q) + (y + r)*width + z*size2D];
							index++;
						}
					}
				}
				std::sort(xxxDATA, xxxDATA + (kernelX3));
				outputArray[x + y*width + z*size2D] = xxxDATA[(kernelX3 - 1) / 2];

				/*for(int q = -distance; q <=distance; q++){
				for (int r = -distance; r <=distance; r++){
				xxxDATA[index] = inputArray[(x+q) +(y+r)*width + z*size2D];
				index++;
				}
				}
				std::sort(xxxDATA, xxxDATA+( kernelX2 ));
				outputArray[x + y*width + z*size2D] = xxxDATA[( kernelX2 - 1 )/2];*/
			}
		}
	}
	delete[] xxxDATA;
}

// Linear smooth filter process
void ml3DThinning::SetLinearSmoothFilterOn()
{
	isNeedLinearSmooth = true;
}

void ml3DThinning::SetLinearSmoothFilterOff()
{
	isNeedLinearSmooth = false;
}

void ml3DThinning::SetLinearSmoothFilter(bool _isNeedLinearSmooth)
{
	isNeedLinearSmooth = _isNeedLinearSmooth;
}

void ml3DThinning::SetLinearSmoothKernel(unsigned int _linearSmoothKernel)
{
	if (_linearSmoothKernel == 6 || _linearSmoothKernel == 18 || _linearSmoothKernel == 26){
		linearSmoothKernel = _linearSmoothKernel;
	}
	else{
		printf("Linear smooth kernel must be 6, 18 or 26\n");
		return;
	}
}

void ml3DThinning::StartLinearSmooth()
{
	if (linearSmoothKernel == 0)
		return;
	if (width == 0 && height == 0 && slice_num == 0)
		return;

	outputArray = new unsigned int[size3D];
	memset(outputArray, 0, sizeof(unsigned int)*size3D);

	std::vector<std::thread> parallelSmooth;
	parallelSmooth.resize(4);

	parallelSmooth[0] = std::thread([this] {LinearSmoothBetween(1, slice_num / 4 + 10); });
	parallelSmooth[1] = std::thread([this] {LinearSmoothBetween(slice_num / 4 - 10, slice_num / 2 + 10); });
	parallelSmooth[2] = std::thread([this] {LinearSmoothBetween(slice_num / 2 - 10, 3 * slice_num / 4 + 10); });
	parallelSmooth[3] = std::thread([this] {LinearSmoothBetween(3 * slice_num / 4 - 10, slice_num - 1); });

	for (auto it = parallelSmooth.begin(); it != parallelSmooth.end(); it++)
		it->join();

	parallelSmooth.clear();

	memset(inputArray, 0, sizeof(unsigned int)*size3D);
	memcpy(inputArray, outputArray, sizeof(unsigned int)*size3D);
	delete[] outputArray;
}

void ml3DThinning::LinearSmoothBetween(unsigned int fromIdx, unsigned int toIdx)
{
	int temp;
	for (unsigned int z = fromIdx; z < toIdx; z++){
		for (unsigned int y = 1; y < height - 1; y++){
			for (unsigned int x = 1; x < width - 1; x++){

				temp = 0;
				for (int p = -1; p <= 1; p++){
					for (int q = -1; q <= 1; q++){
						for (int r = -1; r <= 1; r++){

							switch (linearSmoothKernel)
							{
							case 6:  //平均フィルタ6近傍
							{
										 if ((p == 0 && q == 0 && r == 1) || (p == -1 && q == 0 && r == 0) || (p == 0 && q == 1 && r == 0) || (p == 1 && q == 0 && r == 0) || (p == 0 && q == -1 && r == 0) || (p == 0 && q == 0 && r == 0) || (p == 0 && q == 0 && r == -1))
											 temp = temp + inputArray[(x + p) + (y + q)*width + (z + r)*size2D] / 7;
										 temp = temp + 0 * inputArray[(x + p) + (y + q)*width + (z + r)*size2D];
							} break;

							case 18: //平均フィルタ18近傍
							{
										 if ((p == -1 && q == 1 && r == 1) || (p == 1 && q == 1 && r == 1) || (p == 1 && q == -1 && r == 1) || (p == -1 && q == -1 && r == 1) || (p == -1 && q == 1 && r == -1) || (p == 1 && q == 1 && r == -1) || (p == 1 && q == -1 && r == -1) || (p == -1 && q == -1 && r == -1))
											 temp = temp + 0 * inputArray[(x + p) + (y + q)*width + (z + r)*size2D];
										 temp = temp + inputArray[(x + p) + (y + q)*width + (z + r)*size2D] / 19;
							} break;

							case 26: //平均フィルタ26近傍
							{
										 temp = temp + inputArray[(x + p) + (y + q)*width + (z + r)*size2D] / 27;
							} break;

							default: break;
							}
						}
					}
				}
				outputArray[x + y*width + z*size2D] = temp;
			}
		}
	}
}

// Binarization process
void ml3DThinning::SetBinarizationTheardshold(unsigned int _binarizationTheardshold)
{
	binarizationTheardshold = _binarizationTheardshold;
}

void ml3DThinning::StartBinarization()
{
	//binarizedArray = (unsigned char *)(size3D, sizeof(unsigned char));
	vesselArray = new unsigned char[size3D];
	memset(vesselArray, 0, sizeof(unsigned char)*size3D);

	std::vector<std::thread> parallelBinarize;
	parallelBinarize.resize(4);

	parallelBinarize[0] = std::thread([this] {BinarizationBetween(0, slice_num / 4 + 10); });
	parallelBinarize[1] = std::thread([this] {BinarizationBetween(slice_num / 4 - 10, slice_num / 2 + 10); });
	parallelBinarize[2] = std::thread([this] {BinarizationBetween(slice_num / 2 - 10, 3 * slice_num / 4 + 10); });
	parallelBinarize[3] = std::thread([this] {BinarizationBetween(3 * slice_num / 4 - 10, slice_num); });

	for (auto it = parallelBinarize.begin(); it != parallelBinarize.end(); it++)
		it->join();
}

void ml3DThinning::BinarizationBetween(unsigned int fromIdx, unsigned int toIdx)
{
	for (unsigned int z = fromIdx; z<toIdx; ++z) {
		for (unsigned int y = 0; y<height; ++y) {
			for (unsigned int x = 0; x<width; ++x) {
				vesselArray[x + y*width + z*size2D] = inputArray[x + y*width + z*size2D] > binarizationTheardshold ? 1 : 0;
			}
		}
	}
}

// Nearest smooth filter process
void ml3DThinning::SetNearestSmoothFilterOn()
{
	isNeedNearestSmooth = true;
}

void ml3DThinning::SetNearestSmoothFilterOff()
{
	isNeedNearestSmooth = false;
}

void ml3DThinning::SetNearestSmoothFilter(bool _isNeedNearestSmooth)
{
	isNeedNearestSmooth = _isNeedNearestSmooth;
}

void ml3DThinning::SetNearestSmoothKernel(unsigned int _nearestSmoothKernel)
{
	if (_nearestSmoothKernel % 2 == 0){
		if (_nearestSmoothKernel >= 4)
			nearestSmoothKernel = _nearestSmoothKernel - 1;
	}
	else{
		if (_nearestSmoothKernel >= 3)
			nearestSmoothKernel = _nearestSmoothKernel;
	}
}

void ml3DThinning::StartNearestSmooth()
{
	if (nearestSmoothKernel == 0)
		return;
	if (width == 0 && height == 0 && slice_num == 0)
		return;

	charOutputArray = new unsigned char[size3D];
	memset(charOutputArray, 0, sizeof(unsigned char)*size3D);

	std::vector<std::thread> parallelSmooth;
	parallelSmooth.resize(4);

	parallelSmooth[0] = std::thread([this] {NearestSmoothBetween((nearestSmoothKernel - 1) / 2, slice_num / 4 + 10); });
	parallelSmooth[1] = std::thread([this] {NearestSmoothBetween(slice_num / 4 - 10, slice_num / 2 + 10); });
	parallelSmooth[2] = std::thread([this] {NearestSmoothBetween(slice_num / 2 - 10, 3 * slice_num / 4 + 10); });
	parallelSmooth[3] = std::thread([this] {NearestSmoothBetween(3 * slice_num / 4 - 10, slice_num - (nearestSmoothKernel - 1) / 2); });

	for (auto it = parallelSmooth.begin(); it != parallelSmooth.end(); it++)
		it->join();

	parallelSmooth.clear();

	memset(vesselArray, 0, sizeof(unsigned char)*size3D);
	memcpy(vesselArray, charOutputArray, sizeof(unsigned char)*size3D);
	delete[] charOutputArray;
}

void ml3DThinning::NearestSmoothBetween(unsigned int fromIdx, unsigned int toIdx)
{
	int kernelX3 = pow(nearestSmoothKernel, 3);
	int distance = (nearestSmoothKernel - 1) / 2;

	for (unsigned int z = fromIdx; z < toIdx; z++)
	for (unsigned int y = distance; y < height - distance; y++)
	for (unsigned int x = distance; x < width - distance; x++)
	{
		double heikin = 0;

		for (int p = -distance; p <= distance; p++){
			for (int q = -distance; q <= distance; q++){
				for (int r = -distance; r <= distance; r++){
					heikin = heikin + vesselArray[(x + p) + (y + q)*width + (z + r)*size2D];
				}
			}
		}

		if ((heikin / kernelX3) >= 0 && (heikin / kernelX3)<0.5)
			charOutputArray[x + y*width + z*size2D] = 0;
		else
			charOutputArray[x + y*width + z*size2D] = 1;
	}
}

// 3D closing 
void ml3DThinning::SetClosingOn()
{
	isNeedClosing = true;
}

void ml3DThinning::SetClosingOff()
{
	isNeedClosing = false;
}

void ml3DThinning::SetClosing(bool _isNeedClosing)
{
	isNeedClosing = _isNeedClosing;
}

void ml3DThinning::SetClosingKernel(unsigned int _closingKernel)
{
	closingKernel = _closingKernel;
}

// 3D opening
void ml3DThinning::SetOpeningOn()
{
	isNeedOpening = true;
}

void ml3DThinning::SetOpeningOff()
{
	isNeedOpening = false;
}

void ml3DThinning::SetOpening(bool _isNeedOpening)
{
	isNeedOpening = _isNeedOpening;
}

void ml3DThinning::SetOpeningKernel(unsigned int _openingKernel)
{
	openingKernel = _openingKernel;
}

// Second thinning process
void ml3DThinning::SetSecondThinningOn()
{
	isNeedSecondThinning = true;
}

void ml3DThinning::SetSecondThinningOff()
{
	isNeedSecondThinning = false;
}

void ml3DThinning::SetSecondThinning(bool _isNeedSecondThinning)
{
	isNeedSecondThinning = _isNeedSecondThinning;
}

void ml3DThinning::SetSecondThinningKernel(unsigned int _secondThinningKernel)
{
	secondThinningKernel = _secondThinningKernel;
}

void ml3DThinning::StartClosing()
{
	//
	//first : dilated
	//
	array21 = new unsigned char[size3D]; memset(array21, 0, sizeof(unsigned char)*size3D);
	array22 = new unsigned char[size3D]; memset(array22, 0, sizeof(unsigned char)*size3D);
	array23 = new unsigned char[size3D]; memset(array23, 0, sizeof(unsigned char)*size3D);
	array24 = new unsigned char[size3D]; memset(array24, 0, sizeof(unsigned char)*size3D);

	//pararell dilated
	std::vector<std::thread> parallelDilated;
	parallelDilated.resize(4);
	parallelDilated[0] = std::thread([this] {DilationBetween(closingKernel, 1, slice_num / 4 + 10, 1); });
	parallelDilated[1] = std::thread([this] {DilationBetween(closingKernel, slice_num / 4 - 10, slice_num / 2 + 10, 2); });
	parallelDilated[2] = std::thread([this] {DilationBetween(closingKernel, slice_num / 2 - 10, 3 * slice_num / 4 + 10, 3); });
	parallelDilated[3] = std::thread([this] {DilationBetween(closingKernel, 3 * slice_num / 4 - 10, slice_num - 1, 4); });
	for (auto it = parallelDilated.begin(); it != parallelDilated.end(); it++)
		it->join();
	parallelDilated.clear();

	//pararell join
	std::vector<std::thread> dilatedDataJoin;
	dilatedDataJoin.resize(4);
	dilatedDataJoin[0] = std::thread([this] {Join2VesselArrayBetween(1, slice_num / 4, 1); });
	dilatedDataJoin[1] = std::thread([this] {Join2VesselArrayBetween(slice_num / 4, slice_num / 2, 2); });
	dilatedDataJoin[2] = std::thread([this] {Join2VesselArrayBetween(slice_num / 2, 3 * slice_num / 4, 3); });
	dilatedDataJoin[3] = std::thread([this] {Join2VesselArrayBetween(3 * slice_num / 4, slice_num - 1, 4); });
	for (auto it = dilatedDataJoin.begin(); it != dilatedDataJoin.end(); it++)
		it->join();
	dilatedDataJoin.clear();

	//release memory
	delete[] array21; delete[] array22; delete[] array23; delete[] array24;

	//
	//second : eroded
	//
	array21 = new unsigned char[size3D]; memset(array21, 0, sizeof(unsigned char)*size3D);
	array22 = new unsigned char[size3D]; memset(array22, 0, sizeof(unsigned char)*size3D);
	array23 = new unsigned char[size3D]; memset(array23, 0, sizeof(unsigned char)*size3D);
	array24 = new unsigned char[size3D]; memset(array24, 0, sizeof(unsigned char)*size3D);

	//pararell eroded
	std::vector<std::thread> parallelEroded;
	parallelEroded.resize(4);
	parallelEroded[0] = std::thread([this] {ErosionBetween(closingKernel, 1, slice_num / 4 + 10, 1); });
	parallelEroded[1] = std::thread([this] {ErosionBetween(closingKernel, slice_num / 4 - 10, slice_num / 2 + 10, 2); });
	parallelEroded[2] = std::thread([this] {ErosionBetween(closingKernel, slice_num / 2 - 10, 3 * slice_num / 4 + 10, 3); });
	parallelEroded[3] = std::thread([this] {ErosionBetween(closingKernel, 3 * slice_num / 4 - 10, slice_num - 1, 4); });
	for (auto it = parallelEroded.begin(); it != parallelEroded.end(); it++)
		it->join();
	parallelEroded.clear();

	//pararell joint
	std::vector<std::thread> erodedDataJoin;
	erodedDataJoin.resize(4);
	erodedDataJoin[0] = std::thread([this] {Join2VesselArrayBetween(1, slice_num / 4, 1); });
	erodedDataJoin[1] = std::thread([this] {Join2VesselArrayBetween(slice_num / 4, slice_num / 2, 2); });
	erodedDataJoin[2] = std::thread([this] {Join2VesselArrayBetween(slice_num / 2, 3 * slice_num / 4, 3); });
	erodedDataJoin[3] = std::thread([this] {Join2VesselArrayBetween(3 * slice_num / 4, slice_num - 1, 4); });
	for (auto it = erodedDataJoin.begin(); it != erodedDataJoin.end(); it++)
		it->join();
	erodedDataJoin.clear();

	//release memory
	delete[] array21; delete[] array22; delete[] array23; delete[] array24;
}

void ml3DThinning::StartOpenning()
{
	//
	//first : eroded
	//
	array21 = new unsigned char[size3D]; memset(array21, 0, sizeof(unsigned char)*size3D);
	array22 = new unsigned char[size3D]; memset(array22, 0, sizeof(unsigned char)*size3D);
	array23 = new unsigned char[size3D]; memset(array23, 0, sizeof(unsigned char)*size3D);
	array24 = new unsigned char[size3D]; memset(array24, 0, sizeof(unsigned char)*size3D);

	//pararell eroded
	std::vector<std::thread> parallelEroded;
	parallelEroded.resize(4);
	parallelEroded[0] = std::thread([this] {ErosionBetween(openingKernel, 1, slice_num / 4 + 10, 1); });
	parallelEroded[1] = std::thread([this] {ErosionBetween(openingKernel, slice_num / 4 - 10, slice_num / 2 + 10, 2); });
	parallelEroded[2] = std::thread([this] {ErosionBetween(openingKernel, slice_num / 2 - 10, 3 * slice_num / 4 + 10, 3); });
	parallelEroded[3] = std::thread([this] {ErosionBetween(openingKernel, 3 * slice_num / 4 - 10, slice_num - 1, 4); });
	for (auto it = parallelEroded.begin(); it != parallelEroded.end(); it++)
		it->join();
	parallelEroded.clear();

	//pararell joint
	std::vector<std::thread> erodedDataJoin;
	erodedDataJoin.resize(4);
	erodedDataJoin[0] = std::thread([this] {Join2VesselArrayBetween(1, slice_num / 4, 1); });
	erodedDataJoin[1] = std::thread([this] {Join2VesselArrayBetween(slice_num / 4, slice_num / 2, 2); });
	erodedDataJoin[2] = std::thread([this] {Join2VesselArrayBetween(slice_num / 2, 3 * slice_num / 4, 3); });
	erodedDataJoin[3] = std::thread([this] {Join2VesselArrayBetween(3 * slice_num / 4, slice_num - 1, 4); });
	for (auto it = erodedDataJoin.begin(); it != erodedDataJoin.end(); it++)
		it->join();
	erodedDataJoin.clear();

	//release memory
	delete[] array21; delete[] array22; delete[] array23; delete[] array24;

	//
	//second : dilated
	//
	array21 = new unsigned char[size3D]; memset(array21, 0, sizeof(unsigned char)*size3D);
	array22 = new unsigned char[size3D]; memset(array22, 0, sizeof(unsigned char)*size3D);
	array23 = new unsigned char[size3D]; memset(array23, 0, sizeof(unsigned char)*size3D);
	array24 = new unsigned char[size3D]; memset(array24, 0, sizeof(unsigned char)*size3D);

	//pararell dilated
	std::vector<std::thread> parallelDilated;
	parallelDilated.resize(4);
	parallelDilated[0] = std::thread([this] {DilationBetween(openingKernel, 1, slice_num / 4 + 10, 1); });
	parallelDilated[1] = std::thread([this] {DilationBetween(openingKernel, slice_num / 4 - 10, slice_num / 2 + 10, 2); });
	parallelDilated[2] = std::thread([this] {DilationBetween(openingKernel, slice_num / 2 - 10, 3 * slice_num / 4 + 10, 3); });
	parallelDilated[3] = std::thread([this] {DilationBetween(openingKernel, 3 * slice_num / 4 - 10, slice_num - 1, 4); });
	for (auto it = parallelDilated.begin(); it != parallelDilated.end(); it++)
		it->join();
	parallelDilated.clear();

	//pararell joint
	std::vector<std::thread> dilatedDataJoin;
	dilatedDataJoin.resize(4);
	dilatedDataJoin[0] = std::thread([this] {Join2VesselArrayBetween(1, slice_num / 4, 1); });
	dilatedDataJoin[1] = std::thread([this] {Join2VesselArrayBetween(slice_num / 4, slice_num / 2, 2); });
	dilatedDataJoin[2] = std::thread([this] {Join2VesselArrayBetween(slice_num / 2, 3 * slice_num / 4, 3); });
	dilatedDataJoin[3] = std::thread([this] {Join2VesselArrayBetween(3 * slice_num / 4, slice_num - 1, 4); });
	for (auto it = dilatedDataJoin.begin(); it != dilatedDataJoin.end(); it++)
		it->join();
	dilatedDataJoin.clear();

	//release memory
	delete[] array21; delete[] array22; delete[] array23; delete[] array24;
}

void ml3DThinning::DilationBetween(const int _n, const int fromIdx, const int toIdx, const int type)
{
	for (int z = fromIdx; z<toIdx; z++){
		for (unsigned int y = 1; y<height - 1; y++){
			for (unsigned int x = 1; x<width - 1; x++){
				if (vesselArray[x + y*width + z*size2D] == 0){
					int count = 0;
					for (int vz = -_n; vz <= _n; vz++) {
						for (int vy = -_n; vy <= _n; vy++) {
							for (int vx = -_n; vx <= _n; vx++) {
								if (_n*_n < vx*vx + vy*vy + vz*vz) {
									continue;
								}
								const int xx = x + vx;
								const int yy = y + vy;
								const int zz = z + vz;
								if (xx<1 || xx>width - 1 || yy<1 || yy>height - 1 || zz<1 || zz>slice_num - 1){
									continue;
								}
								if (vesselArray[xx + yy*width + zz*size2D] == 1){
									count++;
								}
							}
						}
					}
					if (count>0){
						if (type == 1)
							array21[x + y*width + z*size2D] = 1;
						else if (type == 2)
							array22[x + y*width + z*size2D] = 1;
						else if (type == 3)
							array23[x + y*width + z*size2D] = 1;
						else if (type == 4)
							array24[x + y*width + z*size2D] = 1;
					}
					else{
						if (type == 1)
							array21[x + y*width + z*size2D] = 0;
						else if (type == 2)
							array22[x + y*width + z*size2D] = 0;
						else if (type == 3)
							array23[x + y*width + z*size2D] = 0;
						else if (type == 4)
							array24[x + y*width + z*size2D] = 0;
					}
				}
				else{
					if (type == 1)
						array21[x + y*width + z*size2D] = 1;
					else if (type == 2)
						array22[x + y*width + z*size2D] = 1;
					else if (type == 3)
						array23[x + y*width + z*size2D] = 1;
					else if (type == 4)
						array24[x + y*width + z*size2D] = 1;
				}
			}
		}
	}
}

void ml3DThinning::ErosionBetween(const int _n, const int fromIdx, const int toIdx, const int type)
{
	for (int z = fromIdx; z<toIdx; z++){
		for (unsigned int y = 1; y<height - 1; y++){
			for (unsigned int x = 1; x<width - 1; x++){
				if (vesselArray[x + y*width + z*size2D] == 1){
					int count = 0;
					for (int vz = -_n; vz <= _n; vz++) {
						for (int vy = -_n; vy <= _n; vy++) {
							for (int vx = -_n; vx <= _n; vx++) {
								if (_n*_n < vx*vx + vy*vy + vz*vz) {
									continue;
								}
								const int xx = x + vx;
								const int yy = y + vy;
								const int zz = z + vz;
								if (xx<1 || xx>width - 1 || yy<1 || yy>height - 1 || zz<1 || zz>slice_num - 1){
									continue;
								}
								if (vesselArray[xx + yy*width + zz*size2D] == 0){
									count++;
								}
							}
						}
					}
					if (count>0){
						if (type == 1)
							array21[x + y*width + z*size2D] = 0;
						else if (type == 2)
							array22[x + y*width + z*size2D] = 0;
						else if (type == 3)
							array23[x + y*width + z*size2D] = 0;
						else if (type == 4)
							array24[x + y*width + z*size2D] = 0;
					}
					else{
						if (type == 1)
							array21[x + y*width + z*size2D] = 1;
						else if (type == 2)
							array22[x + y*width + z*size2D] = 1;
						else if (type == 3)
							array23[x + y*width + z*size2D] = 1;
						else if (type == 4)
							array24[x + y*width + z*size2D] = 1;
					}
				}
				else{
					if (type == 1)
						array21[x + y*width + z*size2D] = 0;
					else if (type == 2)
						array22[x + y*width + z*size2D] = 0;
					else if (type == 3)
						array23[x + y*width + z*size2D] = 0;
					else if (type == 4)
						array24[x + y*width + z*size2D] = 0;
				}
			}
		}
	}
}

void ml3DThinning::Join2VesselArrayBetween(const int fromIdx, const int toIdx, const int type)
{
	for (int z = fromIdx; z<toIdx; z++){
		for (unsigned int y = 1; y<height - 1; y++){
			for (unsigned int x = 1; x<width - 1; x++){
				if (type == 1)
					vesselArray[x + y*width + z*size2D] = array21[x + y*width + z*size2D];
				else if (type == 2)
					vesselArray[x + y*width + z*size2D] = array22[x + y*width + z*size2D];
				else if (type == 3)
					vesselArray[x + y*width + z*size2D] = array23[x + y*width + z*size2D];
				else if (type == 4)
					vesselArray[x + y*width + z*size2D] = array24[x + y*width + z*size2D];
			}
		}
	}
}

// Dilation & Erosion
void ml3DThinning::SetArraySize(unsigned int _w, unsigned int _h, unsigned int _sn)
{
	width = _w; height = _h; slice_num = _sn;
	size2D = width*height;
	size3D = width*height*slice_num;
}

void ml3DThinning::ErosionOnArray(const int erosionSize, unsigned char * inputArray)
{
	//delete[] m_pTemp;
	unsigned char * m_pTemp = new unsigned char[size3D];
	memset(m_pTemp, 0, sizeof(unsigned char)*size3D);

	for (unsigned int z = 1; z<slice_num - 1; z++){
		for (unsigned int y = 1; y<height - 1; y++){
			for (unsigned int x = 1; x<width - 1; x++){
				if (inputArray[x + y*width + z*size2D] == 1){
					int count = 0;
					for (int vz = -erosionSize; vz <= erosionSize; vz++) {
						for (int vy = -erosionSize; vy <= erosionSize; vy++) {
							for (int vx = -erosionSize; vx <= erosionSize; vx++) {
								if (erosionSize*erosionSize < vx*vx + vy*vy + vz*vz) {
									continue;
								}
								const int xx = x + vx;
								const int yy = y + vy;
								const int zz = z + vz;
								if (xx<1 || xx>width - 1 || yy<1 || yy>height - 1 || zz<1 || zz>slice_num - 1){
									continue;
								}
								if (inputArray[xx + yy*width + zz*size2D] == 0){
									count++;
								}
							}
						}
					}
					if (count>0){
						m_pTemp[x + y*width + z*size2D] = 0;
					}
					else{
						m_pTemp[x + y*width + z*size2D] = 1;
					}
				}
				else{
					m_pTemp[x + y*width + z*size2D] = 0;
				}
			}
		}
	}
	for (unsigned int z = 1; z<slice_num - 1; z++){
		for (unsigned int y = 1; y<height - 1; y++){
			for (unsigned int x = 1; x<width - 1; x++){
				inputArray[x + y*width + z*size2D] = m_pTemp[x + y*width + z*size2D];
			}
		}
	}
	delete[] m_pTemp;
}

void ml3DThinning::DilationOnArray(const int dilationSize, unsigned char * inputArray)
{
	//delete[] m_pTemp;
	unsigned char *m_pTemp = new unsigned char[size3D];
	memset(m_pTemp, 0, sizeof(unsigned char)*size3D);

	for (unsigned int z = 1; z<slice_num - 1; z++){
		for (unsigned int y = 1; y<height - 1; y++){
			for (unsigned int x = 1; x<width - 1; x++){
				if (inputArray[x + y*width + z*size2D] == 0){
					int count = 0;
					for (int vz = -dilationSize; vz <= dilationSize; vz++) {
						for (int vy = -dilationSize; vy <= dilationSize; vy++) {
							for (int vx = -dilationSize; vx <= dilationSize; vx++) {
								if (dilationSize*dilationSize < vx*vx + vy*vy + vz*vz) {
									continue;
								}
								const int xx = x + vx;
								const int yy = y + vy;
								const int zz = z + vz;
								if (xx<1 || xx>width - 1 || yy<1 || yy>height - 1 || zz<1 || zz>slice_num - 1){
									continue;
								}
								if (inputArray[xx + yy*width + zz*size2D] == 1){
									count++;
								}
							}
						}
					}
					if (count>0){
						m_pTemp[x + y*width + z*size2D] = 1;
					}
					else{
						m_pTemp[x + y*width + z*size2D] = 0;
					}
				}
				else{
					m_pTemp[x + y*width + z*size2D] = 1;
				}
			}
		}
	}
	for (unsigned int z = 1; z<slice_num - 1; z++){
		for (unsigned int y = 1; y<height - 1; y++){
			for (unsigned int x = 1; x<width - 1; x++){
				inputArray[x + y*width + z*size2D] = m_pTemp[x + y*width + z*size2D];
			}
		}
	}
	delete[] m_pTemp;
}

void ml3DThinning::StartThinning()
{
	//Step1
	//ユークリッド２乗距離変換
	//2013/12/27 Modified by PHAN ( no need to compute distance by each x,y,z dimension, we can compute in one times )
	//vtkSmartPointer<vtkImageData> combineMHD = vtkSmartPointer<vtkImageData>::New();
	vtkImageData* combineMHD = vtkImageData::New();
	combineMHD->SetSpacing(1, 1, 1);
	combineMHD->SetExtent(0, width - 1, 0, height - 1, 0, slice_num - 1);
	combineMHD->AllocateScalars(VTK_UNSIGNED_SHORT, 1);
	for (unsigned int z = 0; z < slice_num; ++z){
		for (unsigned int y = 0; y < height; ++y){
			for (unsigned int x = 0; x < width; ++x){
				if (vesselArray[x + y*width + z*size2D] == 1){
					unsigned short* combine_volume_scalar = static_cast<unsigned short*>(combineMHD->GetScalarPointer(x, y, z));
					*combine_volume_scalar = 11;
				}
			}
		}
	}
	combineMHD->Modified();

	//vtkSmartPointer<vtkImageEuclideanDistance> imageEuclideanDistance = vtkSmartPointer<vtkImageEuclideanDistance>::New();
	vtkImageEuclideanDistance* imageEuclideanDistance = vtkImageEuclideanDistance::New();
	imageEuclideanDistance->SetInputData(combineMHD);
	imageEuclideanDistance->Update();
	convertedImg = imageEuclideanDistance->GetOutput();

	array21 = new unsigned char[size3D]; memset(array21, 0, sizeof(unsigned char)*size3D);
	array22 = new unsigned char[size3D]; memset(array22, 0, sizeof(unsigned char)*size3D);
	array23 = new unsigned char[size3D]; memset(array23, 0, sizeof(unsigned char)*size3D);
	array24 = new unsigned char[size3D]; memset(array24, 0, sizeof(unsigned char)*size3D);

	//pararell thinning
	std::vector<std::thread> parallelThinning;
	parallelThinning.resize(4);
	parallelThinning[0] = std::thread([this] {ThinningBetween(1, slice_num / 4 + 10, 1); });
	parallelThinning[1] = std::thread([this] {ThinningBetween(slice_num / 4 - 10, slice_num / 2 + 10, 2); });
	parallelThinning[2] = std::thread([this] {ThinningBetween(slice_num / 2 - 10, 3 * slice_num / 4 + 10, 3); });
	parallelThinning[3] = std::thread([this] {ThinningBetween(3 * slice_num / 4 - 10, slice_num - 1, 4); });
	for (auto it = parallelThinning.begin(); it != parallelThinning.end(); it++)
		it->join();
	parallelThinning.clear();

	//pararell joint to vessel array
	std::vector<std::thread> thinnedDataJoin;
	thinnedDataJoin.resize(4);
	thinnedDataJoin[0] = std::thread([this] {Join2VesselArrayBetween(1, slice_num / 4, 1); });
	thinnedDataJoin[1] = std::thread([this] {Join2VesselArrayBetween(slice_num / 4, slice_num / 2, 2); });
	thinnedDataJoin[2] = std::thread([this] {Join2VesselArrayBetween(slice_num / 2, 3 * slice_num / 4, 3); });
	thinnedDataJoin[3] = std::thread([this] {Join2VesselArrayBetween(3 * slice_num / 4, slice_num - 1, 4); });
	for (auto it = thinnedDataJoin.begin(); it != thinnedDataJoin.end(); it++)
		it->join();
	thinnedDataJoin.clear();

	//release memory
	delete[] array21; delete[] array22; delete[] array23; delete[] array24;
	combineMHD->Delete(); imageEuclideanDistance->Delete();
}

void ml3DThinning::ThinningBetween(const int fromSlice, const int toSlice, const int type)
{
	unsigned char * array2ForParallel = new unsigned char[size3D];
	memset(array2ForParallel, 0, sizeof(unsigned char)*size3D);

	int * arrayLabeling = new int[100];
	memset(arrayLabeling, 0, sizeof(int)* 27);

	std::list<List> l_Obj;
	std::list<List>::iterator it;
	std::list<List>::iterator it2;
	std::list<List>::iterator it3;
	l_Obj.clear();

	int nmax = 0, nmin = 100000;
	for (auto z = fromSlice; z<toSlice; ++z){
		for (unsigned int y = 0; y < height; ++y){
			for (unsigned int x = 0; x < width; ++x){
				int mainIndex = x + y*width + z*size2D;
				double *scalar = static_cast<double*>(convertedImg->GetScalarPointer(x, y, z));

				if (scalar[0] != 0){
					//array2ForParallel[x + y*w + z*wh] = 21;
					array2ForParallel[mainIndex] = (int)scalar[0] + 20;
					//printf("<%d> ",(int)scalar[0]+20);
				}

				nmax = array2ForParallel[mainIndex] > nmax ? array2ForParallel[mainIndex] : nmax; //ユークリッド距離値の最大値→nmax
				nmin = array2ForParallel[mainIndex] < nmin ? array2ForParallel[mainIndex] : nmin; //ユークリッド距離値の最小値→d
			}
		}
	}

	//Step2 初期境界画素の検出
	int fcount = 0;
	for (int z = fromSlice; z<toSlice; z++){
		for (unsigned int y = 1; y<height - 1; y++){
			for (unsigned int x = 1; x<width - 1; x++){
				fcount = 0;
				for (int k = 0; k<6; k++){
					int xx = x + dir6[k * 3];
					int yy = y + dir6[k * 3 + 1];
					int zz = z + dir6[k * 3 + 2];
					if (array2ForParallel[xx + yy*width + zz*size2D] == 0){
						fcount++;
					}
				}
				if (fcount != 0 && array2ForParallel[x + y*width + z*size2D]>20){
					struct List p;
					p.i = x;
					p.j = y;
					p.k = z;
					p.l = array2ForParallel[x + y*width + z*size2D];
					l_Obj.push_back(p);//リストに加える
					array2ForParallel[x + y*width + z*size2D] = 1;
				}
			}
		}
	}

	int nList = 0;
	int mm = 0;
	//ループ
	do{
		//Step3　永久保存点の検出と境界画素の分類
		for (it = l_Obj.begin(); it != l_Obj.end();) {
			if (it->l <= nmin){
				if (Deletable(it->i, it->j, it->k, array2ForParallel, arrayLabeling, type) == 0){
					it->l = 16;//一次保存
					it++;
				}
				else{
					int count26 = 0;
					for (int k = 0; k<26; k++){
						int xx = it->i + dir26[k * 3];
						int yy = it->j + dir26[k * 3 + 1];
						int zz = it->k + dir26[k * 3 + 2];
						if (array2ForParallel[xx + yy*width + zz*size2D] != 0){
							count26++;
						}
					}
					if (count26 == 1){//端点
						it = l_Obj.erase(it);//リストから消去
					}
					else{
						int m = count26;
						it->l = (int)(m / 3) + 7;
						it++;
					}
				}
			}
			else it++;
		}

		//Step4　画像の消去
		for (int bordertype = 7; bordertype<16; bordertype++){
			for (it2 = l_Obj.begin(); it2 != l_Obj.end();){
				if (it2->l == bordertype){
					if (Deletable(it2->i, it2->j, it2->k, array2ForParallel, arrayLabeling, type) == 0){
						it2->l = 16;//一次保存
						it2++;
					}
					else{
						int count26 = 0;
						for (int m = 0; m<26; m++){
							int xx = it2->i + dir26[m * 3];
							int yy = it2->j + dir26[m * 3 + 1];
							int zz = it2->k + dir26[m * 3 + 2];
							if (array2ForParallel[xx + yy*width + zz*size2D] != 0){
								count26++;
							}
						}
						if (count26 == 1){//端点
							it2 = l_Obj.erase(it2);//リストから消去
						}
						else{
							int ix = it2->i;
							int iy = it2->j;
							int iz = it2->k;
							array2ForParallel[ix + iy*width + iz*size2D] = 0;
							it2 = l_Obj.erase(it2);//リストから消去
							for (int m = 0; m<6; m++){
								int xx = ix + dir6[m * 3];
								int yy = iy + dir6[m * 3 + 1];
								int zz = iz + dir6[m * 3 + 2];
								if (array2ForParallel[xx + yy*width + zz*size2D]>20){
									struct List p;
									p.i = xx;
									p.j = yy;
									p.k = zz;
									p.l = array2ForParallel[xx + yy*width + zz*size2D];
									l_Obj.push_back(p);//リストに加える
									array2ForParallel[xx + yy*width + zz*size2D] = 1;
								}
							}
						}
					}
				}
				else it2++;
			}
		}

		//Step5　終了判定
		nList = 0;
		mm = 0;
		nmin = 10000000;
		for (it3 = l_Obj.begin(); it3 != l_Obj.end(); it3++){
			if (it3->l == 16){
				mm++;
			}
			nList++;
			if (it3->l > 20){
				if (it3->l < nmin){
					nmin = it3->l;
				}
			}
		}
	} while (nmin<nmax || mm != nList);

	delete[] array2ForParallel;
	delete[] arrayLabeling;
}

int ml3DThinning::Deletable(int _x, int _y, int _z, unsigned char *Array2ForParalell, int *ArrayLabelling, const int type)
{
	int count6 = 0;
	int count18 = 0;
	int count26 = 0;

	for (int k = 0; k<26; k++){
		int xx = _x + dir26[k * 3];
		int yy = _y + dir26[k * 3 + 1];
		int zz = _z + dir26[k * 3 + 2];
		if (Array2ForParalell[xx + yy*width + zz*size2D] != 0){
			if (type == 1) {
				array21[xx + yy*width + zz*size2D] = 1;
			}
			else if (type == 2) {
				array22[xx + yy*width + zz*size2D] = 1;
			}
			else if (type == 3) {
				array23[xx + yy*width + zz*size2D] = 1;
			}
			else if (type == 4) {
				array24[xx + yy*width + zz*size2D] = 1;
			}
		}
		else{
			if (type == 1) {
				array21[xx + yy*width + zz*size2D] = 0;
			}
			else if (type == 2) {
				array22[xx + yy*width + zz*size2D] = 0;
			}
			else if (type == 3) {
				array23[xx + yy*width + zz*size2D] = 0;
			}
			else if (type == 4) {
				array24[xx + yy*width + zz*size2D] = 0;
			}
		}
	}
	if (type == 1) {
		array21[_x + _y*width + _z*size2D] = 1;
	}
	else if (type == 2) {
		array22[_x + _y*width + _z*size2D] = 1;
	}
	else if (type == 3) {
		array23[_x + _y*width + _z*size2D] = 1;
	}
	else if (type == 4) {
		array24[_x + _y*width + _z*size2D] = 1;
	}

	for (int k = 0; k<6; k++){
		int xx = _x + dir6[k * 3];
		int yy = _y + dir6[k * 3 + 1];
		int zz = _z + dir6[k * 3 + 2];

		if (type == 1) {
			if (array21[xx + yy*width + zz*size2D] == 1)
				count6++;
		}
		else if (type == 2) {
			if (array22[xx + yy*width + zz*size2D] == 1)
				count6++;
		}
		else if (type == 3) {
			if (array23[xx + yy*width + zz*size2D] == 1)
				count6++;
		}
		else if (type == 4) {
			if (array24[xx + yy*width + zz*size2D] == 1)
				count6++;
		}
	}
	for (int k = 0; k<18; k++){
		int xx = _x + dir18[k * 3];
		int yy = _y + dir18[k * 3 + 1];
		int zz = _z + dir18[k * 3 + 2];

		if (type == 1) {
			if (array21[xx + yy*width + zz*size2D] == 1)
				count18++;
		}
		else if (type == 2) {
			if (array22[xx + yy*width + zz*size2D] == 1)
				count18++;
		}
		else if (type == 3) {
			if (array23[xx + yy*width + zz*size2D] == 1)
				count18++;
		}
		else if (type == 4) {
			if (array24[xx + yy*width + zz*size2D] == 1)
				count18++;
		}
	}
	for (int k = 0; k<26; k++){
		int xx = _x + dir26[k * 3];
		int yy = _y + dir26[k * 3 + 1];
		int zz = _z + dir26[k * 3 + 2];

		if (type == 1) {
			if (array21[xx + yy*width + zz*size2D] == 1)
				count26++;
		}
		else if (type == 2) {
			if (array22[xx + yy*width + zz*size2D] == 1)
				count26++;
		}
		else if (type == 3) {
			if (array23[xx + yy*width + zz*size2D] == 1)
				count26++;
		}
		else if (type == 4) {
			if (array24[xx + yy*width + zz*size2D] == 1)
				count26++;
		}
	}
	////端点かどうか
	//if(count26==1){
	//	Array3[((z*m_nSizey+y)*m_nSizex+x)]=2;//永久保存点
	//	return 0;
	//}	
	//空洞の個数Ym
	if (count6 == 6){
		//Array3[((z*m_nSizey+y)*m_nSizex+x)]=2;//永久保存点
		return 0;
	}

	//連結成分の個数Rm
	int* piTemp1 = new int[27];
	int* piTemp2 = new int[27];
	memset(piTemp1, 0, sizeof(int)* 27);
	memset(piTemp2, 0, sizeof(int)* 27);
	for (int k = 0; k<26; k++){
		if (type == 1) {
			piTemp1[((dir26[k * 3 + 2] + 1) * 3 + (dir26[k * 3 + 1] + 1)) * 3 + (dir26[k * 3] + 1)] = array21[((dir26[k * 3 + 2] + _z)*height + (dir26[k * 3 + 1] + _y))*width + (dir26[k * 3] + _x)];
		}
		else if (type == 2) {
			piTemp1[((dir26[k * 3 + 2] + 1) * 3 + (dir26[k * 3 + 1] + 1)) * 3 + (dir26[k * 3] + 1)] = array22[((dir26[k * 3 + 2] + _z)*height + (dir26[k * 3 + 1] + _y))*width + (dir26[k * 3] + _x)];
		}
		else if (type == 3) {
			piTemp1[((dir26[k * 3 + 2] + 1) * 3 + (dir26[k * 3 + 1] + 1)) * 3 + (dir26[k * 3] + 1)] = array23[((dir26[k * 3 + 2] + _z)*height + (dir26[k * 3 + 1] + _y))*width + (dir26[k * 3] + _x)];
		}
		else if (type == 4) {
			piTemp1[((dir26[k * 3 + 2] + 1) * 3 + (dir26[k * 3 + 1] + 1)) * 3 + (dir26[k * 3] + 1)] = array24[((dir26[k * 3 + 2] + _z)*height + (dir26[k * 3 + 1] + _y))*width + (dir26[k * 3] + _x)];
		}
	}
	piTemp1[13] = 1;
	//memset(ArrayLabelling,0,sizeof(int)*27);
	LiLabeling3D(3, 3, 3, piTemp1, ArrayLabelling);
	for (int k = 0; k<26; k++){
		if (ArrayLabelling[((dir26[k * 3 + 2] + 1) * 3 + (dir26[k * 3 + 1] + 1)) * 3 + (dir26[k * 3] + 1)] != ArrayLabelling[13]){
			continue;
		}
		else{
			piTemp2[((dir26[k * 3 + 2] + 1) * 3 + (dir26[k * 3 + 1] + 1)) * 3 + (dir26[k * 3] + 1)] = 1;
		}
	}
	int Rm = LiLabeling3D(3, 3, 3, piTemp2, piTemp1);
	delete[] piTemp1;
	delete[] piTemp2;
	if (Rm != 1){
		//Array3[((z*m_nSizey+y)*m_nSizex+x)]=2;//永久保存点
		return 0;
	}
	//穴の個数Hm
	//連結数の計算
	int Ncount = 2;
	const int dir[] = { 1, 0, 1, 1, 0, 1, -1, 1, -1, 0, -1, -1, 0, -1, 1, -1, 1, 0 };

	for (int vz = -1; vz <= 1; vz += 2){
		if (type == 1) {
			Ncount -= (1 - array21[(((_z + vz)*height + _y)*width + _x)]);//1項
		}
		else if (type == 2) {
			Ncount -= (1 - array22[(((_z + vz)*height + _y)*width + _x)]);//1項
		}
		else if (type == 3) {
			Ncount -= (1 - array23[(((_z + vz)*height + _y)*width + _x)]);//1項
		}
		else if (type == 4) {
			Ncount -= (1 - array24[(((_z + vz)*height + _y)*width + _x)]);//1項
		}
		for (int k = 0; k<4; k++){
			int xx = _x + dir[k * 4];
			int yy = _y + dir[k * 4 + 1];

			if (type == 1) {
				Ncount += (1 - array21[(((_z + vz)*height + _y)*width + _x)])*
					(1 - array21[(((_z + vz)*height + yy)*width + xx)])*
					(1 - array21[((_z*height + yy)*width + xx)]);//2項
			}
			else if (type == 2) {
				Ncount += (1 - array22[(((_z + vz)*height + _y)*width + _x)])*
					(1 - array22[(((_z + vz)*height + yy)*width + xx)])*
					(1 - array22[((_z*height + yy)*width + xx)]);//2項
			}
			else if (type == 3) {
				Ncount += (1 - array23[(((_z + vz)*height + _y)*width + _x)])*
					(1 - array23[(((_z + vz)*height + yy)*width + xx)])*
					(1 - array23[((_z*height + yy)*width + xx)]);//2項
			}
			else if (type == 4) {
				Ncount += (1 - array24[(((_z + vz)*height + _y)*width + _x)])*
					(1 - array24[(((_z + vz)*height + yy)*width + xx)])*
					(1 - array24[((_z*height + yy)*width + xx)]);//2項
			}
		}
	}

	for (int k = 0; k<4; k++){
		int dx1 = _x + dir[k * 4];
		int dy1 = _y + dir[k * 4 + 1];

		if (type == 1) {
			Ncount -= (1 - array21[((_z*height + dy1)*width + dx1)]);//3項
		}
		else if (type == 2) {
			Ncount -= (1 - array22[((_z*height + dy1)*width + dx1)]);//3項
		}
		else if (type == 3) {
			Ncount -= (1 - array23[((_z*height + dy1)*width + dx1)]);//3項
		}
		else if (type == 4) {
			Ncount -= (1 - array24[((_z*height + dy1)*width + dx1)]);//3項
		}

		int dx2 = _x + dir[k * 4 + 2];
		int dy2 = _y + dir[k * 4 + 3];
		int dx3 = _x + dir[k * 4 + 4];
		int dy3 = _y + dir[k * 4 + 5];

		if (type == 1) {
			Ncount += (1 - array21[((_z*height + dy1)*width + dx1)])*
				(1 - array21[((_z*height + dy2)*width + dx2)])*
				(1 - array21[((_z*height + dy3)*width + dx3)]);//4項
		}
		else if (type == 2) {
			Ncount += (1 - array22[((_z*height + dy1)*width + dx1)])*
				(1 - array22[((_z*height + dy2)*width + dx2)])*
				(1 - array22[((_z*height + dy3)*width + dx3)]);//4項
		}
		else if (type == 3) {
			Ncount += (1 - array23[((_z*height + dy1)*width + dx1)])*
				(1 - array23[((_z*height + dy2)*width + dx2)])*
				(1 - array23[((_z*height + dy3)*width + dx3)]);//4項
		}
		else if (type == 4) {
			Ncount += (1 - array24[((_z*height + dy1)*width + dx1)])*
				(1 - array24[((_z*height + dy2)*width + dx2)])*
				(1 - array24[((_z*height + dy3)*width + dx3)]);//4項
		}

		for (int vz = -1; vz <= 1; vz += 2){
			if (type == 1) {
				Ncount -= (1 - array21[((_z*height + dy1)*width + dx1)])*
					(1 - array21[((_z*height + dy2)*width + dx2)])*
					(1 - array21[((_z*height + dy3)*width + dx3)])*
					(1 - array21[(((_z + vz)*height + _y)*width + _x)])*
					(1 - array21[(((_z + vz)*height + dy1)*width + dx1)])*
					(1 - array21[(((_z + vz)*height + dy2)*width + dx2)])*
					(1 - array21[(((_z + vz)*height + dy3)*width + dx3)]);
			}
			else if (type == 2) {
				Ncount -= (1 - array22[((_z*height + dy1)*width + dx1)])*
					(1 - array22[((_z*height + dy2)*width + dx2)])*
					(1 - array22[((_z*height + dy3)*width + dx3)])*
					(1 - array22[(((_z + vz)*height + _y)*width + _x)])*
					(1 - array22[(((_z + vz)*height + dy1)*width + dx1)])*
					(1 - array22[(((_z + vz)*height + dy2)*width + dx2)])*
					(1 - array22[(((_z + vz)*height + dy3)*width + dx3)]);
			}
			else if (type == 3) {
				Ncount -= (1 - array23[((_z*height + dy1)*width + dx1)])*
					(1 - array23[((_z*height + dy2)*width + dx2)])*
					(1 - array23[((_z*height + dy3)*width + dx3)])*
					(1 - array23[(((_z + vz)*height + _y)*width + _x)])*
					(1 - array23[(((_z + vz)*height + dy1)*width + dx1)])*
					(1 - array23[(((_z + vz)*height + dy2)*width + dx2)])*
					(1 - array23[(((_z + vz)*height + dy3)*width + dx3)]);
			}
			else if (type == 4) {
				Ncount -= (1 - array24[((_z*height + dy1)*width + dx1)])*
					(1 - array24[((_z*height + dy2)*width + dx2)])*
					(1 - array24[((_z*height + dy3)*width + dx3)])*
					(1 - array24[(((_z + vz)*height + _y)*width + _x)])*
					(1 - array24[(((_z + vz)*height + dy1)*width + dx1)])*
					(1 - array24[(((_z + vz)*height + dy2)*width + dx2)])*
					(1 - array24[(((_z + vz)*height + dy3)*width + dx3)]);
			}
		}
	}

	//int Nm=Ncount;
	//int Hm=1-Nm;
	if (Ncount == 1){
		//m_pArray[((z*m_nSizey+y)*m_nSizex+x)]=0;//消去可能
		//bThinning=true;//ループ
		return 1;
	}
	else{
		//Array3[((z*m_nSizey+y)*m_nSizex+x)]=2;//永久保存点
		return 0;
	}
}

int ml3DThinning::LiLabeling3D(const int _sizex, const int _sizey, const int _sizez, const int *pBin, int *pObj)
{
	int* checked = new int[_sizex*_sizey*_sizez];
	memset(checked, 0, sizeof(int)*_sizex*_sizey*_sizez);
	std::list<int> lObj;
	std::list<int>::iterator itrYet;

	int objNumber = 1;
	memset(pObj, 0, sizeof(int)*_sizex*_sizey*_sizez);

	for (int z = 0; z<_sizez; z++){
		for (int y = 0; y<_sizey; y++) {
			for (int x = 0; x<_sizex; x++) {
				// オブジェクト判定
				if (pBin[(z*_sizey + y)*_sizex + x] == 0) {
					continue;
				}
				// 評価済み判定
				if (checked[(z*_sizey + y)*_sizex + x] == 1) {
					continue;
				}
				// Li探索開始
				lObj.clear();
				lObj.push_back((z*_sizey + y)*_sizex + x);
				// 評価済みフラグ
				checked[(z*_sizey + y)*_sizex + x] = 1;
				for (itrYet = lObj.begin(); itrYet != lObj.end(); itrYet++) {
					const int ind = *itrYet;
					const int x0 = ind%_sizex;
					const int y0 = (ind / _sizex) % _sizey;
					const int z0 = ind / _sizex / _sizey;
					for (int k = 0; k<26; k++){
						int xx = x0 + dir26[k * 3];
						int yy = y0 + dir26[k * 3 + 1];
						int zz = z0 + dir26[k * 3 + 2];
						if (xx<0 || yy<0 || zz<0 || xx>_sizex - 1 || yy>_sizey - 1 || zz>_sizez - 1){
							continue;
						}
						// オブジェクト判定
						if (pBin[(zz*_sizey + yy)*_sizex + xx] == 0) {
							continue;
						}
						// 評価済み判定
						if (checked[(zz*_sizey + yy)*_sizex + xx] == 1) {
							continue;
						}
						// 新規追加
						lObj.push_back((zz*_sizey + yy)*_sizex + xx);
						// 評価済みフラグ
						checked[(zz*_sizey + yy)*_sizex + xx] = 1;
					}
				}
				// オブジェクト番号を登録
				for (itrYet = lObj.begin(); itrYet != lObj.end(); itrYet++) {
					const int ind = *itrYet;
					const int xx = ind%_sizex;
					const int yy = (ind / _sizex) % _sizey;
					const int zz = ind / _sizex / _sizey;
					pObj[(zz*_sizey + yy)*_sizex + xx] = objNumber;
				}
				lObj.clear();
				objNumber++;
			}
		}
	}
	delete[] checked;
	return objNumber - 1;
}

int ml3DThinning::LiLabeling3D_2(const int _sizex, const int _sizey, const int _sizez, const int* pBin, int* pObj)
{
	int* checked = new int[_sizex*_sizey*_sizez];
	memset(checked, 0, sizeof(int)*_sizex*_sizey*_sizez);
	std::list<int> lObj;
	std::list<int>::iterator itrYet;

	int objNumber = 1;
	memset(pObj, 0, sizeof(int)*_sizex*_sizey*_sizez);

	for (int z = 0; z<_sizez; z++){
		for (int y = 0; y<_sizey; y++) {
			for (int x = 0; x<_sizex; x++) {
				// オブジェクト判定
				if (pBin[(z*_sizey + y)*_sizex + x] == 0) {
					continue;
				}
				// 評価済み判定
				if (checked[(z*_sizey + y)*_sizex + x] == 1) {
					continue;
				}
				// Li探索開始
				lObj.clear();
				lObj.push_back((z*_sizey + y)*_sizex + x);
				// 評価済みフラグ
				checked[(z*_sizey + y)*_sizex + x] = 1;
				for (itrYet = lObj.begin(); itrYet != lObj.end(); itrYet++) {
					const int ind = *itrYet;
					const int x0 = ind%_sizex;
					const int y0 = (ind / _sizex) % _sizey;
					const int z0 = ind / _sizex / _sizey;
					for (int k = 0; k<6; k++){
						int xx = x0 + dir26[k * 3];
						int yy = y0 + dir26[k * 3 + 1];
						int zz = z0 + dir26[k * 3 + 2];
						if (xx<0 || yy<0 || zz<0 || xx>_sizex - 1 || yy>_sizey - 1 || zz>_sizez - 1){
							continue;
						}
						// オブジェクト判定
						if (pBin[(zz*_sizey + yy)*_sizex + xx] == 0) {
							continue;
						}
						// 評価済み判定
						if (checked[(zz*_sizey + yy)*_sizex + xx] == 1) {
							continue;
						}
						// 新規追加
						lObj.push_back((zz*_sizey + yy)*_sizex + xx);
						// 評価済みフラグ
						checked[(zz*_sizey + yy)*_sizex + xx] = 1;
					}
				}
				// オブジェクト番号を登録
				for (itrYet = lObj.begin(); itrYet != lObj.end(); itrYet++) {
					const int ind = *itrYet;
					const int xx = ind%_sizex;
					const int yy = (ind / _sizex) % _sizey;
					const int zz = ind / _sizex / _sizey;
					pObj[(zz*_sizey + yy)*_sizex + xx] = objNumber;
				}
				lObj.clear();
				objNumber++;
			}
		}
	}
	delete[] checked;
	return objNumber - 1;
}

void ml3DThinning::DetectionBifurcationPoints()
{
	bifurcationPoints.clear();
	bifurcationArray = new unsigned char[size3D];
	memset(bifurcationArray, 0, sizeof(unsigned char)*size3D);

	vtkSmartPointer<vtkPoints> thinnedPointsWithSpacing = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> thinnedPointsWithoutSpacing = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkUnsignedCharArray> thinnedPointsColour = vtkSmartPointer<vtkUnsignedCharArray>::New();
	thinnedPointsColour->SetName("colors");
	thinnedPointsColour->SetNumberOfComponents(3);

	vtkSmartPointer<vtkStringArray> thinnedPointsType = vtkSmartPointer<vtkStringArray>::New();
	thinnedPointsType->SetName("type");

	for (unsigned int z = 1; z<slice_num - 1; z++){
		for (unsigned int y = 1; y<height - 1; y++){
			for (unsigned int x = 1; x<width - 1; x++){
				if (vesselArray[x + y*width + z*size2D] == 1){
					//連結成分の個数Rm
					int* piTemp1 = new int[27];
					int* piTemp2 = new int[27];
					memset(piTemp1, 0, sizeof(int)* 27);
					memset(piTemp2, 0, sizeof(int)* 27);
					for (int k = 0; k<26; k++){
						piTemp1[((dir26[k * 3 + 2] + 1) * 3 + (dir26[k * 3 + 1] + 1)) * 3 + (dir26[k * 3] + 1)] = vesselArray[((dir26[k * 3 + 2] + z)*height + (dir26[k * 3 + 1] + y))*width + (dir26[k * 3] + x)];
					}
					piTemp1[13] = 1;

					int *arrayForLabeling = new int[27]; memset(arrayForLabeling, 0, sizeof(int)* 27);
					memset(arrayForLabeling, 0, sizeof(int)* 27);
					LiLabeling3D(3, 3, 3, piTemp1, arrayForLabeling);

					for (int k = 0; k<26; k++){
						if (arrayForLabeling[((dir26[k * 3 + 2] + 1) * 3 + (dir26[k * 3 + 1] + 1)) * 3 + (dir26[k * 3] + 1)] != arrayForLabeling[13]){
							continue;
						}
						else{
							piTemp2[((dir26[k * 3 + 2] + 1) * 3 + (dir26[k * 3 + 1] + 1)) * 3 + (dir26[k * 3] + 1)] = 1;
						}
					} delete[] arrayForLabeling;

					int Rm = LiLabeling3D_2(3, 3, 3, piTemp2, piTemp1);
					delete[] piTemp1;
					delete[] piTemp2;
					if (Rm >= 3){
						bifurcationPoints.emplace_back(x*spacing[0], y*spacing[1], z*spacing[2]);
						thinnedPointsWithSpacing->InsertNextPoint(x*spacing[0], y*spacing[1], z*spacing[2]);
						thinnedPointsWithoutSpacing->InsertNextPoint(x, y, z);
						unsigned char col[3] = { 255, 255, 0 };
						thinnedPointsColour->InsertNextTupleValue(col);
						thinnedPointsType->InsertNextValue("bifurcation");
						bifurcationArray[x + y*width + z*size2D] = 1;

					}
					else if (Rm == 1/* || Rm == 0*/) {
						thinnedPointsWithSpacing->InsertNextPoint(x*spacing[0], y*spacing[1], z*spacing[2]);
						thinnedPointsWithoutSpacing->InsertNextPoint(x, y, z);
						unsigned char col[3] = { 0, 255, 0 };
						thinnedPointsColour->InsertNextTupleValue(col);
						thinnedPointsType->InsertNextValue("end");
						bifurcationArray[x + y*width + z*size2D] = 2;

					}
					else {
						thinnedPointsWithSpacing->InsertNextPoint(x*spacing[0], y*spacing[1], z*spacing[2]);
						thinnedPointsWithoutSpacing->InsertNextPoint(x, y, z);
						unsigned char col[3] = { 255, 0, 0 };
						thinnedPointsColour->InsertNextTupleValue(col);
						thinnedPointsType->InsertNextValue("vessel");
					}
				}
			}
		}
	}

	thinnedDataWithSpacing = vtkSmartPointer<vtkPolyData>::New();
	thinnedDataWithSpacing->SetPoints(thinnedPointsWithSpacing);
	thinnedDataWithSpacing->GetPointData()->SetScalars(thinnedPointsColour);
	thinnedDataWithSpacing->GetPointData()->AddArray(thinnedPointsType);
	thinnedDataWithSpacing->Modified();

	thinnedDataWithoutSpacing = vtkSmartPointer<vtkPolyData>::New();
	thinnedDataWithoutSpacing->SetPoints(thinnedPointsWithoutSpacing);
	thinnedDataWithoutSpacing->GetPointData()->SetScalars(thinnedPointsColour);
	thinnedDataWithoutSpacing->GetPointData()->AddArray(thinnedPointsType);
	thinnedDataWithoutSpacing->Modified();
}

void ml3DThinning::ReDetectBifurcation()
{
	for (int z = 0; z<slice_num; ++z){
		for (int y = 0; y<height; ++y){
			for (int x = 0; x<width; ++x){
				if (vesselArray[x + y*width + z*size2D] == 1){
					if (neigborNumberFrom(ml::vector3d(x, y, z)) == 3){
						bifurcationArray[x + y*width + z*size2D] = 1;
					}
					else if (neigborNumberFrom(ml::vector3d(x, y, z)) == 1){
						bifurcationArray[x + y*width + z*size2D] = 2;
					}
				}
			}
		}
	}

	errorPaternLenght = 3;
	ml3DThinning::DeleteErrorFromBifurcation();
	ml3DThinning::DetectionBifurcationPoints();
}

int ml3DThinning::neigborNumberFrom(ml::vector3d p)
{
	const int size1D[3] = { width, height, slice_num };
	for (int i = 0; i<3; ++i){
		if (p[i] <= 0){
			return false;
		}
	}
	for (int i = 0; i<3; ++i){
		if (p[i] >= size1D[i] - 1){
			return false;
		}
	}

	int num = 0;
	for (int z = p[2] - 1; z <= p[2] + 1; ++z){
		for (int y = p[1] - 1; y <= p[1] + 1; ++y){
			for (int x = p[0] - 1; x <= p[0] + 1; ++x){
				//printf ("<%d, %d, %d> ",x,y,z);
				if (vesselArray[x + y*size1D[0] + z*size2D] == 0){
					continue;
				}
				if (x == p[0] && y == p[1] && z == p[2]){
					continue;
				}
				num++;
			}
		}
	}
	return num;
}

void ml3DThinning::DetectBifurcationFromArray(unsigned char *srcArray, unsigned char* dstArray, unsigned int _w, unsigned int _h, unsigned int _sn)
{
	width = _w;
	height = _h;
	slice_num = _sn;
	size2D = _w*_h;
	size3D = _w*_h*_sn;

	vesselArray = srcArray;
	bifurcationArray = dstArray;

	for (unsigned int z = 2; z<_sn - 2; z++){
		for (unsigned int y = 2; y<_h - 2; y++){
			for (unsigned int x = 2; x<_w - 2; x++){
				if (srcArray[x + y*_w + z*size2D] == 1){
					//連結成分の個数Rm
					int* piTemp1 = new int[27];
					int* piTemp2 = new int[27];
					memset(piTemp1, 0, sizeof(int)* 27);
					memset(piTemp2, 0, sizeof(int)* 27);
					for (int k = 0; k<26; k++){
						piTemp1[((dir26[k * 3 + 2] + 1) * 3 + (dir26[k * 3 + 1] + 1)) * 3 + (dir26[k * 3] + 1)] = vesselArray[((dir26[k * 3 + 2] + z)*height + (dir26[k * 3 + 1] + y))*width + (dir26[k * 3] + x)];
					}
					piTemp1[13] = 1;

					int *arrayForLabeling = new int[27]; memset(arrayForLabeling, 0, sizeof(int)* 27);
					memset(arrayForLabeling, 0, sizeof(int)* 27);
					LiLabeling3D(3, 3, 3, piTemp1, arrayForLabeling);

					for (int k = 0; k<26; k++){
						if (arrayForLabeling[((dir26[k * 3 + 2] + 1) * 3 + (dir26[k * 3 + 1] + 1)) * 3 + (dir26[k * 3] + 1)] != arrayForLabeling[13]){
							continue;
						}
						else{
							piTemp2[((dir26[k * 3 + 2] + 1) * 3 + (dir26[k * 3 + 1] + 1)) * 3 + (dir26[k * 3] + 1)] = 1;
						}
					} delete[] arrayForLabeling;

					int Rm = LiLabeling3D_2(3, 3, 3, piTemp2, piTemp1);
					delete[] piTemp1;
					delete[] piTemp2;
					if (Rm >= 3){
						dstArray[x + y*width + z*size2D] = 1;
					}
					else if (Rm == 1/* || Rm == 0*/) {
						dstArray[x + y*width + z*size2D] = 2;
					}
					else {
						dstArray[x + y*width + z*size2D] = 0;
					}
				}
			}
		}
	}
}

// Second thinning process
void ml3DThinning::SetDeleteErrorPaternOn()
{
	isNeedDeleteError = true;
}

void ml3DThinning::SetDeleteErrorPaternOff()
{
	isNeedDeleteError = false;
}

void ml3DThinning::SetDeleteErrorPatern(bool _isNeedDeleteError) // Default is ON
{
	isNeedDeleteError = _isNeedDeleteError;
}

void ml3DThinning::SetErrorPaternLenght(unsigned int _errorPaternLenght) // Default is 3
{
	if (_errorPaternLenght > 0)
		errorPaternLenght = _errorPaternLenght;
}

void ml3DThinning::DetectErrorFromArray(unsigned char *_binariArray, unsigned char* _vesselArray, unsigned int _w, unsigned int _h, unsigned int _sn)
{
	std::cout << "DetectErrorFromArray in\n";
	width = _w;
	height = _h;
	slice_num = _sn;
	size2D = _w*_h;
	size3D = _w*_h*_sn;

	vesselArray = _binariArray;
	bifurcationArray = _vesselArray;

	vtkSmartPointer<vtkPoints> thinnedPointsWithoutSpacing = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkUnsignedCharArray> thinnedPointsColour = vtkSmartPointer<vtkUnsignedCharArray>::New();
	thinnedPointsColour->SetName("colors");
	thinnedPointsColour->SetNumberOfComponents(3);

	vtkSmartPointer<vtkStringArray> thinnedPointsType = vtkSmartPointer<vtkStringArray>::New();
	thinnedPointsType->SetName("type");

	for (unsigned int z = 1; z<slice_num - 1; z++){
		for (unsigned int y = 1; y<height - 1; y++){
			for (unsigned int x = 1; x<width - 1; x++){
				if (!IsValidPoint(x, y, z)) { continue; }
				if (bifurcationArray[x + y*width + z*size2D] == 1){
					thinnedPointsWithoutSpacing->InsertNextPoint(x, y, z);
					unsigned char col[3] = { 255, 255, 0 };
					thinnedPointsColour->InsertNextTupleValue(col);
					thinnedPointsType->InsertNextValue("bifurcation");
				}
				else if (bifurcationArray[x + y*width + z*size2D] == 2){
					thinnedPointsWithoutSpacing->InsertNextPoint(x, y, z);
					unsigned char col[3] = { 0, 255, 0 };
					thinnedPointsColour->InsertNextTupleValue(col);
					thinnedPointsType->InsertNextValue("end");
				}
				else{
					if (vesselArray[x + y*width + z*size2D] == 1){
						thinnedPointsWithoutSpacing->InsertNextPoint(x, y, z);
						unsigned char col[3] = { 255, 0, 0 };
						thinnedPointsColour->InsertNextTupleValue(col);
						thinnedPointsType->InsertNextValue("vessel");
					}
				}
			}
		}
	}

	thinnedDataWithoutSpacing = vtkSmartPointer<vtkPolyData>::New();
	thinnedDataWithoutSpacing->SetPoints(thinnedPointsWithoutSpacing);
	thinnedDataWithoutSpacing->GetPointData()->SetScalars(thinnedPointsColour);
	thinnedDataWithoutSpacing->GetPointData()->AddArray(thinnedPointsType);
	thinnedDataWithoutSpacing->Modified();

	//std::cout << "DetectErrorFromArray out\n";
	ml3DThinning::ErrorDeleteProcess();
}

void ml3DThinning::NeighborPointsFromPoint(const ml::vector3d startPoint, ml::points3d &foundedPoints, ml::points3d &edgePoints)
{
	for (int i = 0; i<3; ++i){
		if (startPoint[i] <= 0)
			return;
	}

	unsigned int size1D[3] = { width, height, slice_num };
	for (int i = 0; i<3; ++i){
		if (startPoint[i] >= size1D[i] - 1)
			return;
	}

	bool avail = false;
	for (int i = 0; i<edgePoints.size(); ++i){
		if (startPoint[0] == edgePoints[i][0] && startPoint[1] == edgePoints[i][1] && startPoint[2] == edgePoints[i][2]){
			avail = true;
			break;
		}
	}
	if (!avail)
		edgePoints.push_back(startPoint);

	for (int z = startPoint[2] - 1; z <= startPoint[2] + 1; ++z){
		for (int y = startPoint[1] - 1; y <= startPoint[1] + 1; ++y){
			for (int x = startPoint[0] - 1; x <= startPoint[0] + 1; ++x){

				if (x == startPoint[0] && y == startPoint[1] && z == startPoint[2])
					continue;

				avail = false;
				for (int i = 0; i<edgePoints.size(); i++){
					if (x == edgePoints[i][0] && y == edgePoints[i][1] && z == edgePoints[i][2]){
						avail = true;
						break;
					}
				}
				if (!avail) {
					if (vesselArray[x + y*width + z*size2D] == 1)
						foundedPoints.push_back(ml::vector3d(x, y, z));
				}
			}
		}
	}
}

void ml3DThinning::FindNeighborPointsFromPoint(const ml::vector3d fromPoint, ml::points3d &edgePoints)
{
	ml::points3d foundedPoints;
	ml3DThinning::NeighborPointsFromPoint(fromPoint, foundedPoints, edgePoints);
	for (int i = 0; i<foundedPoints.size(); ++i)
		ml3DThinning::FindNeighborPointsFromPoint(foundedPoints[i], edgePoints);
}

void ml3DThinning::DeleteEdgesFromPoint(unsigned char *_binariArray, unsigned char* _vesselArray, unsigned int _w, unsigned int _h, unsigned int _sn, ml::vector3d startP)
{
	width = _w;
	height = _h;
	slice_num = _sn;
	size2D = _w*_h;
	size3D = _w*_h*_sn;

	vesselArray = _binariArray;
	bifurcationArray = _vesselArray;

	ml::points3d edgePoints;
	ml3DThinning::FindNeighborPointsFromPoint(startP, edgePoints);
	for (int i = 0; i<edgePoints.size(); ++i){
		ml::vector3d p = edgePoints[i];
		vesselArray[(int)p[0] + (int)p[1] * width + (int)p[2] * size2D] = 0;
		bifurcationArray[(int)p[0] + (int)p[1] * width + (int)p[2] * size2D] = 0;
	}
}

void ml3DThinning::DeletePointsInArray(unsigned char *_binariArray, unsigned char* _vesselArray, unsigned int _w, unsigned int _h, unsigned int _sn, ml::points3d deletePoints)
{
	if (deletePoints.size() <= 0) { return; }
	unsigned int _size2D = _w*_h;
	for (int i = 0; i<deletePoints.size(); ++i){
		ml::vector3d p = deletePoints[i];
		_vesselArray[(int)p[0] + (int)p[1] * _w + (int)p[2] * _size2D] = 0;
		_binariArray[(int)p[0] + (int)p[1] * _w + (int)p[2] * _size2D] = 0;
	}
}

bool ml3DThinning::IsValidPoint(const int x, const int y, const int z)
{
	int size1D[3] = { width, height, slice_num };
	int aPoint[3] = { x, y, z };
	for (int i = 0; i<3; ++i){
		if (aPoint[i] <= 2){ return false; }
	}
	for (int i = 0; i<3; ++i){
		if (aPoint[i] >= size1D[i] - 2){ return false; }
	}
	return true;
}

void ml3DThinning::ErrorDeleteProcess()
{
	//std::cout << "ErrorDeleteProcess in\n";
	vtkPoints *thinnedPointsWithoutSpacing = thinnedDataWithoutSpacing->GetPoints();
	vtkStringArray *thinnedPointsType = (vtkStringArray*)thinnedDataWithoutSpacing->GetPointData()->GetAbstractArray("type");

	for (vtkIdType id = 0; id < thinnedPointsWithoutSpacing->GetNumberOfPoints(); id++){
		if (thinnedPointsType->GetValue(id) == "end"){

			oneEdgePoints4Detect = vtkSmartPointer<vtkPoints>::New();
			double aPoint[3]; thinnedPointsWithoutSpacing->GetPoint(id, aPoint);

			int endPoint[3] = { (int)aPoint[0], (int)aPoint[1], (int)aPoint[2] };
			bool av = false; ml3DThinning::DetectError(endPoint, av);
			if (av){
				//printf("Avaible error with lenght %d \n",oneEdgeLength4Detect);
				for (vtkIdType _id = 0; _id<oneEdgePoints4Detect->GetNumberOfPoints(); _id++)
				{
					double deletePoint[3]; oneEdgePoints4Detect->GetPoint(_id, deletePoint);
					if (IsValidPoint((int)deletePoint[0], (int)deletePoint[1], (int)deletePoint[2])){
						vesselArray[(int)deletePoint[0] + (int)deletePoint[1] * width + (int)deletePoint[2] * size2D] = 0;
					}
				}
			}
		}
	}

	//bifurcation detection again
	delete[] bifurcationArray;
	ml3DThinning::DetectionBifurcationPoints();

	thinnedPointsWithoutSpacing = thinnedDataWithoutSpacing->GetPoints();
	thinnedPointsType = (vtkStringArray*)thinnedDataWithoutSpacing->GetPointData()->GetAbstractArray("type");

	// check for need continue deletion or not
	bool needContinue = false;
	for (vtkIdType id = 0; id < thinnedPointsWithoutSpacing->GetNumberOfPoints(); id++){
		if (thinnedPointsType->GetValue(id) == "end"){

			oneEdgePoints4Detect = vtkSmartPointer<vtkPoints>::New();
			double aPoint[3]; thinnedPointsWithoutSpacing->GetPoint(id, aPoint);

			int endPoint[3] = { (int)aPoint[0], (int)aPoint[1], (int)aPoint[2] };
			bool av = false; ml3DThinning::DetectError(endPoint, av);
			if (av){
				needContinue = true;
				id = thinnedPointsWithoutSpacing->GetNumberOfPoints() + 10;
			}
		}
	}
	if (needContinue){
		ml3DThinning::ErrorDeleteProcess();
	}
	else if (!needContinue){
		ml3DThinning::DeleteErrorFromBifurcation();
	}
	//std::cout << "ErrorDeleteProcess out\n";
}

void ml3DThinning::ErrorDeleteProcess2()
{
	// 2 : Detect error pattern in converted graph, and insert into the error_list
	vtkPoints *thinnedPointsWithoutSpacing = thinnedDataWithoutSpacing->GetPoints();
	vtkStringArray *thinnedPointsType = (vtkStringArray*)thinnedDataWithoutSpacing->GetPointData()->GetAbstractArray("type");
	for (vtkIdType id = 0; id < thinnedPointsWithoutSpacing->GetNumberOfPoints(); id++){
		if (thinnedPointsType->GetValue(id) == "end"){

			oneEdgePoints4Detect = vtkSmartPointer<vtkPoints>::New();
			double aPoint[3]; thinnedPointsWithoutSpacing->GetPoint(id, aPoint);

			int endPoint[3] = { (int)aPoint[0], (int)aPoint[1], (int)aPoint[2] };
			bool av = false; ml3DThinning::DetectError(endPoint, av);
			if (av){
				Error anError = Error(oneEdgePoints4Detect);
				if (!IsAvaibleErrorInArray(anError, errorList)){
					errorList.push_back(anError);
				}
			}
		}
	}


	// 3 : Check each error in the error_list is trully error or not
	ml3DThinning::CheckTrullyError();


	// 4 : Delete trully error in the error_list
	delete[] bifurcationArray;
	ml3DThinning::DetectionBifurcationPoints();


	// 5 : Check still avaible error or not. If avaible rerun this function
	/*bool needContinueDetect = false;
	thinnedPointsWithoutSpacing = thinnedDataWithoutSpacing->GetPoints();
	thinnedPointsType = (vtkStringArray*)thinnedDataWithoutSpacing->GetPointData()->GetAbstractArray("type");
	for (vtkIdType id = 0; id < thinnedPointsWithoutSpacing->GetNumberOfPoints(); id++){
	if (thinnedPointsType->GetValue(id) == "end"){

	oneEdgePoints4Detect = vtkSmartPointer<vtkPoints>::New();
	double aPoint[3]; thinnedPointsWithoutSpacing->GetPoint(id, aPoint);

	int endPoint[3] = { (int)aPoint[0], (int)aPoint[1], (int)aPoint[2] };
	bool av = false; ml3DThinning::DetectError(endPoint, av);
	if (av){
	Error anError = Error(oneEdgePoints4Detect);
	if (!IsAvaibleErrorInArray(anError, errorList)){
	needContinueDetect = true;
	}
	}
	}
	}
	if (needContinueDetect){
	ml3DThinning::ErrorDeleteProcess2();
	} else {*/
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputConnection(appendPoly->GetOutputPort());
	writer->SetFileName("vesselCircle.vtp");
	writer->SetCompressorTypeToNone();
	writer->Write();
	//}
}

void ml3DThinning::CheckTrullyError()
{
	//cout << "\nError list num : " << errorList.size() << endl;
	vector<Error> _errors = errorList;
	errorList.clear();


	for (size_t index = 0; index < _errors.size(); index++)
	{
		Error anError = _errors[index];
		currentErrorPointList = anError.pointList;
		if (currentErrorPointList->GetNumberOfPoints() <= 5){
			anError.needDelete = true;
			for (vtkIdType id = 0; id < anError.pointList->GetNumberOfPoints(); id++)
			{
				double p[3]; anError.pointList->GetPoint(id, p);
				//if (!IsValidPoint((int)p[0], (int)p[1], (int)p[2])) { continue; }
				vesselArray[(int)p[0] + (int)p[1] * width + (int)p[2] * size2D] = 0;
			}
			continue;
		}

		/*for (vtkIdType id = 0; id < currentErrorPointList->GetNumberOfPoints(); id++){
		const double *p = currentErrorPointList->GetPoint(id);
		printf("currentErrorPointList %.0f %.0f %.0f \n", p[0], p[1], p[2]);
		}*/

		vtkSmartPointer<vtkPoints> detectedPointList = vtkSmartPointer<vtkPoints>::New();
		for (vtkIdType id = 0; id < currentErrorPointList->GetNumberOfPoints() - 2; id++)
		{
			//double vesselRadius = ml3DThinning::RadiusWithCenter(currentErrorPointList->GetPoint(id));
			//double vesselRadius = ml3DThinning::RadiusWithCenter2(currentErrorPointList->GetPoint(id), detectedPointList);
			double vesselRadius = ml3DThinning::RadiusWithCenter3(currentErrorPointList->GetPoint(id), detectedPointList);
			anError.averageRadius += vesselRadius;
			anError.radiusList.push_back(vesselRadius);
			//cout << id << " : " << vesselRadius << endl;
		}
		anError.averageRadius = anError.averageRadius / anError.radiusList.size();
		//cout << "Average : " << anError.averageRadius << endl;


		if (anError.radiusList.size() > 2){//3){
			for (size_t idx = 1; idx < anError.radiusList.size() - 1; idx++)
			{
				if (anError.averageRadius <= 0){
					for (vtkIdType id = 0; id < anError.pointList->GetNumberOfPoints(); id++)
					{
						double p[3]; anError.pointList->GetPoint(id, p);
						vesselArray[(int)p[0] + (int)p[1] * width + (int)p[2] * size2D] = 0;
					}
					anError.needDelete = true;
					//cout << "Going delete error ....\n";
					break;
				}
				else{
					double previousRadius = anError.radiusList[idx - 1];
					double currentRadius = anError.radiusList[idx];
					//double nextRadius = anError.radiusList[idx+1];

					if (fabs(currentRadius - previousRadius) > radiusThreshold){//0.31){//&& nextRadius - currentRadius > 0.2){
						for (vtkIdType id = 0; id < anError.pointList->GetNumberOfPoints(); id++)
						{
							double p[3]; anError.pointList->GetPoint(id, p);
							vesselArray[(int)p[0] + (int)p[1] * width + (int)p[2] * size2D] = 0;
						}
						anError.needDelete = true;
						//cout << "Going delete error ....\n";
						break;
					}
				}
				//if (anError.averageRadius <= 0 || (fabs(vesselRadius - anError.averageRadius) >= 0.4)){
				//	anError.needDelete = true;
				//	for (vtkIdType id = 0; id < anError.pointList->GetNumberOfPoints(); id++)
				//	{
				//		double p[3]; anError.pointList->GetPoint(id, p);
				//		//if (!IsValidPoint((int)p[0], (int)p[1], (int)p[2])) { continue; }
				//		vesselArray[(int)p[0] + (int)p[1] * width + (int)p[2] * size2D] = 0;
				//	}
				//	break;
				//}
			}
		}

		if (!anError.needDelete){
			errorList.push_back(anError);
		}
		//cout << endl;
	}
}

double ml3DThinning::RadiusWithCenter(const double center[3])
{
	double returnRadius = 0;
	vector3d direction(0, 0, 0);
	vtkSmartPointer<vtkPoints> pointList = vtkSmartPointer<vtkPoints>::New();
	ml3DThinning::DirectionVectorFromPoint(center, pointList, direction);

	if (!isZeroPoint(direction)){
		unsigned int boundIdx = 0;
		double _center[3] = { 0, 0, 0 };
		for (int i = 0; i<3; ++i){
			_center[i] = center[i] * spacing[i];
			_center[i] = _center[i] + bound[boundIdx];
			boundIdx += 2;
		}
		returnRadius = VesselRadius(direction, vector3d(_center[0], _center[1], _center[2]), bspTree);

		vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
		polygonSource->SetNumberOfSides(20);
		polygonSource->SetCenter(_center[0], _center[1], _center[2]);
		polygonSource->SetRadius(returnRadius);
		polygonSource->SetNormal(direction[0], direction[1], direction[2]);
		polygonSource->Update();

		appendPoly->AddInputConnection(polygonSource->GetOutputPort());
		appendPoly->Modified();
		appendPoly->Update();
	}

	return returnRadius;
}

double ml3DThinning::RadiusWithCenter2(const double center[3], vtkSmartPointer<vtkPoints> detectedPointList)
{
	double returnRadius = 0;
	vector3d direction(0, 0, 0);
	vtkSmartPointer<vtkPoints> pointList = vtkSmartPointer<vtkPoints>::New();
	ml3DThinning::DirectionVectorFromPoint2(center, pointList, detectedPointList, direction);

	if (!isZeroPoint(direction)){
		unsigned int boundIdx = 0;
		double _center[3] = { 0, 0, 0 };
		for (int i = 0; i<3; ++i){
			_center[i] = center[i] * spacing[i];
			_center[i] = _center[i] + bound[boundIdx];
			boundIdx += 2;
		}
		returnRadius = VesselRadius(direction, vector3d(_center[0], _center[1], _center[2]), bspTree);

		vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
		polygonSource->SetNumberOfSides(20);
		polygonSource->SetCenter(_center[0], _center[1], _center[2]);
		polygonSource->SetRadius(returnRadius);
		polygonSource->SetNormal(direction[0], direction[1], direction[2]);
		polygonSource->Update();

		appendPoly->AddInputConnection(polygonSource->GetOutputPort());
		appendPoly->Modified();
		appendPoly->Update();
	}
	else {
		cout << "in\n";
		::system("PAUSE");
	}

	return returnRadius;
}

double ml3DThinning::RadiusWithCenter3(const double center[3], vtkSmartPointer<vtkPoints> detectedPointList)
{
	double returnRadius = 0;
	vector3d direction(0, 0, 0), foot(center[0], center[1], center[2]);
	unsigned int boundIdx = 0;
	for (int i = 0; i<3; ++i){
		foot[i] = foot[i] * spacing[i];
		foot[i] = foot[i] + bound[boundIdx];
		boundIdx += 2;
	}

	vtkSmartPointer<vtkPoints> pointList = vtkSmartPointer<vtkPoints>::New();
	ml3DThinning::DirectionVectorFromPoint3(center, pointList, detectedPointList, foot, direction);

	if (!isZeroPoint(direction)){
		returnRadius = VesselRadius(direction, foot, bspTree);

		vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
		polygonSource->SetNumberOfSides(20);
		polygonSource->SetCenter(foot[0], foot[1], foot[2]);
		polygonSource->SetRadius(returnRadius);
		polygonSource->SetNormal(direction[0], direction[1], direction[2]);
		polygonSource->Update();

		appendPoly->AddInputConnection(polygonSource->GetOutputPort());
		appendPoly->Modified();
		appendPoly->Update();
	}
	else {
		cout << "in\n";
		//::system("PAUSE");
	}

	return returnRadius;
}

void ml3DThinning::DirectionVectorFromPoint(const double startP[3], vtkSmartPointer<vtkPoints> pointList, ml::vector3d &direction)
{
	if (!IsValidPoint((int)startP[0], (int)startP[1], (int)startP[2])) { return; }

	int length = (int)currentErrorPointList->GetNumberOfPoints();
	if (length > 3)
		length = 3;

	if (pointList->GetNumberOfPoints() >= length){
		ml::points3d points4Fitting;
		for (vtkIdType id = 0; id < pointList->GetNumberOfPoints(); id++){
			double point[3]; pointList->GetPoint(id, point);
			unsigned int boundIdx = 0;
			for (int i = 0; i<3; ++i){
				point[i] = point[i] * spacing[i];
				point[i] = point[i] + bound[boundIdx];
				boundIdx += 2;
			}
			points4Fitting.push_back(ml::vector3d(point[0], point[1], point[2]));
		}
		ml::line fittingLine;
		fittingLine.fit(points4Fitting);
		direction = fittingLine.direction();
		return;
	}

	if (!ExistedPointInSet(ml::vector3d(startP[0], startP[1], startP[2]), pointList)){
		pointList->InsertNextPoint(startP[0], startP[1], startP[2]);
	}

	if (pointList->GetNumberOfPoints() >= length){
		ml::points3d points4Fitting;
		for (vtkIdType id = 0; id < pointList->GetNumberOfPoints(); id++){
			double point[3]; pointList->GetPoint(id, point);
			unsigned int boundIdx = 0;
			for (int i = 0; i<3; ++i){
				point[i] = (point[i] * spacing[i]);
				point[i] = point[i] + bound[boundIdx];
				boundIdx += 2;
			}
			points4Fitting.push_back(ml::vector3d(point[0], point[1], point[2]));
		}
		ml::line fittingLine;
		fittingLine.fit(points4Fitting);
		direction = fittingLine.direction();
		return;
	}

	ml::points3d neighborhoodPoints;
	for (int z = startP[2] - 1; z <= startP[2] + 1; z++){
		for (int y = startP[1] - 1; y <= startP[1] + 1; y++){
			for (int x = startP[0] - 1; x <= startP[0] + 1; x++){

				if (!IsValidPoint(x, y, z)) { continue; }
				if (x == startP[0] && y == startP[1] && z == startP[2]){ continue; }
				if (vesselArray[x + y*width + z*size2D] == 1) {
					if (ExistedPointInSet(vector3d(x, y, z), currentErrorPointList) && !ExistedPointInSet(vector3d(x, y, z), pointList)){
						neighborhoodPoints.push_back(vector3d(x, y, z));
					}
				}
			}
		}
	}

	for (size_t index = 0; index < neighborhoodPoints.size(); index++){
		vector3d _p = neighborhoodPoints[index];
		const double p[3] = { _p[0], _p[1], _p[2] };
		ml3DThinning::DirectionVectorFromPoint(p, pointList, direction);
	}
}

void ml3DThinning::DirectionVectorFromPoint2(const double startP[3], vtkSmartPointer<vtkPoints> pointList, vtkSmartPointer<vtkPoints> detectedPointList, vector3d &direction)
{
	if (!IsValidPoint((int)startP[0], (int)startP[1], (int)startP[2])) { return; }

	int length = (int)currentErrorPointList->GetNumberOfPoints();
	if (length > 3)
		length = 3;

	if (pointList->GetNumberOfPoints() >= length){
		ml::points3d points4Fitting;
		for (vtkIdType id = 0; id < pointList->GetNumberOfPoints(); id++){
			double point[3]; pointList->GetPoint(id, point);
			unsigned int boundIdx = 0;
			for (int i = 0; i<3; ++i){
				point[i] = point[i] * spacing[i];
				point[i] = point[i] + bound[boundIdx];
				boundIdx += 2;
			}
			points4Fitting.push_back(ml::vector3d(point[0], point[1], point[2]));
		}
		ml::line fittingLine;
		fittingLine.fit(points4Fitting);
		direction = fittingLine.direction();
		return;
	}

	if (!ExistedPointInSet(ml::vector3d(startP[0], startP[1], startP[2]), pointList)){
		pointList->InsertNextPoint(startP[0], startP[1], startP[2]);
	}

	if (pointList->GetNumberOfPoints() >= length){
		ml::points3d points4Fitting;
		for (vtkIdType id = 0; id < pointList->GetNumberOfPoints(); id++){
			double point[3]; pointList->GetPoint(id, point);
			unsigned int boundIdx = 0;
			for (int i = 0; i<3; ++i){
				point[i] = (point[i] * spacing[i]);
				point[i] = point[i] + bound[boundIdx];
				boundIdx += 2;
			}
			points4Fitting.push_back(ml::vector3d(point[0], point[1], point[2]));
		}
		ml::line fittingLine;
		fittingLine.fit(points4Fitting);
		direction = fittingLine.direction();
		return;
	}

	ml::points3d neighborhoodPoints;
	int x1 = startP[0] - 1, y1 = startP[1] - 1, z1 = startP[2] - 1;
	int x2 = startP[0] + 1, y2 = startP[1] + 1, z2 = startP[2] + 1;

	for (int z = z1; z <= z2; z++){
		for (int y = y1; y <= y2; y++){
			for (int x = x1; x <= x2; x++){
				if (!IsValidPoint(x, y, z)) { continue; }
				if (x == startP[0] && y == startP[1] && z == startP[2]){ continue; }
				if (vesselArray[x + y*width + z*size2D] == 1) {
					if (ExistedPointInSet(vector3d(x, y, z), currentErrorPointList) && !ExistedPointInSet(vector3d(x, y, z), pointList) && !ExistedPointInSet(vector3d(x, y, z), detectedPointList)){
						neighborhoodPoints.push_back(vector3d(x, y, z));
					}
				}
			}
		}
	}

	/*if (neighborhoodPoints.size() == 0){
	cout << "size 0 ";
	}*/

	if (pointList->GetNumberOfPoints() == 1){
		detectedPointList->InsertNextPoint(pointList->GetPoint(0));
	}

	for (size_t index = 0; index < neighborhoodPoints.size(); index++){
		vector3d _p = neighborhoodPoints[index];
		const double p[3] = { _p[0], _p[1], _p[2] };
		ml3DThinning::DirectionVectorFromPoint2(p, pointList, detectedPointList, direction);
	}
}

void ml3DThinning::DirectionVectorFromPoint3(const double startP[3], vtkSmartPointer<vtkPoints> pointList, vtkSmartPointer<vtkPoints> detectedPointList, vector3d &foot, vector3d &direction)
{
	if (!IsValidPoint((int)startP[0], (int)startP[1], (int)startP[2])) { return; }

	int length = (int)currentErrorPointList->GetNumberOfPoints();
	if (length > 3)
		length = 3;

	if (pointList->GetNumberOfPoints() >= length){
		ml::points3d points4Fitting;
		for (vtkIdType id = 0; id < pointList->GetNumberOfPoints(); id++){
			double point[3]; pointList->GetPoint(id, point);
			unsigned int boundIdx = 0;
			for (int i = 0; i<3; ++i){
				point[i] = point[i] * spacing[i];
				point[i] = point[i] + bound[boundIdx];
				boundIdx += 2;
			}
			points4Fitting.push_back(ml::vector3d(point[0], point[1], point[2]));
		}
		ml::line fittingLine;
		fittingLine.fit(points4Fitting);
		direction = fittingLine.direction();
		foot = fittingLine.perpendicular_foot(foot);
		return;
	}

	if (!ExistedPointInSet(ml::vector3d(startP[0], startP[1], startP[2]), pointList)){
		pointList->InsertNextPoint(startP[0], startP[1], startP[2]);
	}

	if (pointList->GetNumberOfPoints() >= length){
		ml::points3d points4Fitting;
		for (vtkIdType id = 0; id < pointList->GetNumberOfPoints(); id++){
			double point[3]; pointList->GetPoint(id, point);
			unsigned int boundIdx = 0;
			for (int i = 0; i<3; ++i){
				point[i] = (point[i] * spacing[i]);
				point[i] = point[i] + bound[boundIdx];
				boundIdx += 2;
			}
			points4Fitting.push_back(ml::vector3d(point[0], point[1], point[2]));
		}
		ml::line fittingLine;
		fittingLine.fit(points4Fitting);
		direction = fittingLine.direction();
		foot = fittingLine.perpendicular_foot(foot);
		return;
	}

	ml::points3d neighborhoodPoints;
	int x1 = startP[0] - 1, y1 = startP[1] - 1, z1 = startP[2] - 1;
	int x2 = startP[0] + 1, y2 = startP[1] + 1, z2 = startP[2] + 1;

	for (int z = z1; z <= z2; z++){
		for (int y = y1; y <= y2; y++){
			for (int x = x1; x <= x2; x++){
				if (!IsValidPoint(x, y, z)) { continue; }
				if (x == startP[0] && y == startP[1] && z == startP[2]){ continue; }
				if (vesselArray[x + y*width + z*size2D] == 1) {
					if (ExistedPointInSet(vector3d(x, y, z), currentErrorPointList) && !ExistedPointInSet(vector3d(x, y, z), pointList) && !ExistedPointInSet(vector3d(x, y, z), detectedPointList)){
						neighborhoodPoints.push_back(vector3d(x, y, z));
					}
				}
			}
		}
	}

	if (pointList->GetNumberOfPoints() == 1){
		detectedPointList->InsertNextPoint(pointList->GetPoint(0));
	}

	for (size_t index = 0; index < neighborhoodPoints.size(); index++){
		vector3d _p = neighborhoodPoints[index];
		const double p[3] = { _p[0], _p[1], _p[2] };
		ml3DThinning::DirectionVectorFromPoint3(p, pointList, detectedPointList, foot, direction);
	}
}

void ml3DThinning::DetectError(const int aPoint[3], bool &av)
{
	if (!IsValidPoint((int)aPoint[0], (int)aPoint[1], (int)aPoint[2])) { return; }

	if (!ExistedPointInSet(ml::vector3d(aPoint[0], aPoint[1], aPoint[2]), oneEdgePoints4Detect)){
		oneEdgePoints4Detect->InsertNextPoint(aPoint[0], aPoint[1], aPoint[2]);
	}

	ml::points3d neighborhoodPoints;
	for (int z = aPoint[2] - 1; z <= aPoint[2] + 1; z++){
		for (int y = aPoint[1] - 1; y <= aPoint[1] + 1; y++){
			for (int x = aPoint[0] - 1; x <= aPoint[0] + 1; x++){

				if (!IsValidPoint(x, y, z)) { continue; }
				if (x == aPoint[0] && y == aPoint[1] && z == aPoint[2]){ continue; }

				if (vesselArray[x + y*width + z*size2D] == 1) {
					if (!ExistedPointInSet(ml::vector3d(x, y, z), oneEdgePoints4Detect)){
						ml::vector3d aPoint(x, y, z);
						neighborhoodPoints.emplace_back(aPoint);
					}
				}
			}
		}
	}

	if (neighborhoodPoints.size() == 1)
	{
		ml::vector3d aPoint = neighborhoodPoints[0];
		int x = (int)aPoint[0], y = (int)aPoint[1], z = (int)aPoint[2];
		if (!IsValidPoint(x, y, z)) { return; }

		if (bifurcationArray[x + y*width + z*size2D] == 1 || bifurcationArray[x + y*width + z*size2D] == 2){
			if (oneEdgePoints4Detect->GetNumberOfPoints() <= errorPaternLenght){
				if (bifurcationArray[x + y*width + z*size2D] == 2){
					oneEdgePoints4Detect->InsertNextPoint(x, y, z);
				}
				av = true;
				return;
			}
		}
		else{
			if (oneEdgePoints4Detect->GetNumberOfPoints() < errorPaternLenght){
				int nextPoint[3] = { x, y, z };
				ml3DThinning::DetectError(nextPoint, av);
			}
		}

	}
	else
	{
		for (unsigned int idx = 0; idx<neighborhoodPoints.size(); idx++){
			ml::vector3d aPoint = neighborhoodPoints[idx];
			int x = (int)aPoint[0], y = (int)aPoint[1], z = (int)aPoint[2];
			if (!IsValidPoint(x, y, z)) { continue; }

			if (bifurcationArray[x + y*width + z*size2D] == 1){
				if (oneEdgePoints4Detect->GetNumberOfPoints() <= errorPaternLenght){
					av = true;
					return;
				}
			}
		}
	}

}

void ml3DThinning::DeleteErrorFromBifurcation()
{
	vtkPoints *thinnedPointsWithoutSpacing = thinnedDataWithoutSpacing->GetPoints();
	vtkStringArray *thinnedPointsType = (vtkStringArray*)thinnedDataWithoutSpacing->GetPointData()->GetAbstractArray("type");

	for (vtkIdType id = 0; id < thinnedPointsWithoutSpacing->GetNumberOfPoints(); id++){
		if (thinnedPointsType->GetValue(id) == "bifurcation"){

			ml::points3d vessel;
			double aPoint[3]; thinnedPointsWithoutSpacing->GetPoint(id, aPoint);
			for (int z = aPoint[2] - 1; z <= aPoint[2] + 1; z++){
				for (int y = aPoint[1] - 1; y <= aPoint[1] + 1; y++){
					for (int x = aPoint[0] - 1; x <= aPoint[0] + 1; x++){
						if (x == aPoint[0] && y == aPoint[1] && z == aPoint[2]){ continue; }
						if (!IsValidPoint(x, y, z)) { continue; }
						if (vesselArray[x + y*width + z*size2D] == 1) { vessel.emplace_back(ml::vector3d(x, y, z)); }
					}
				}
			}

			if (vessel.size() >= 3)
			{
				detectedPoints = vtkSmartPointer<vtkPoints>::New();
				detectedPoints->InsertNextPoint(aPoint[0], aPoint[1], aPoint[2]);
				for (unsigned int idx = 0; idx<vessel.size(); idx++)
				{
					ml::vector3d point = vessel[idx];
					detectedPoints->InsertNextPoint(point[0], point[1], point[2]);
				}
				for (unsigned int idx = 0; idx<vessel.size(); idx++)
				{
					ml::vector3d point = vessel[idx];
					oneEdgePoints4Detect = vtkSmartPointer<vtkPoints>::New();

					int endPoint[3] = { (int)point[0], (int)point[1], (int)point[2] };
					int length = 0; ml3DThinning::DetectErrorFromBifurcation(endPoint, length);
					//printf("Avaible vessel error with lenght %d \n",oneEdgeLength4Detect);

					if (length == 1){
						for (vtkIdType _id = 0; _id<oneEdgePoints4Detect->GetNumberOfPoints(); _id++)
						{
							double deletePoint[3]; oneEdgePoints4Detect->GetPoint(_id, deletePoint);
							if (!IsValidPoint((int)deletePoint[0], (int)deletePoint[1], (int)deletePoint[2])) { continue; }
							vesselArray[(int)deletePoint[0] + (int)deletePoint[1] * width + (int)deletePoint[2] * size2D] = 0;
						}
					}
				}
			}
		}
	}

	//bifurcation detection again
	delete[] bifurcationArray;
	ml3DThinning::DetectionBifurcationPoints();

	thinnedPointsWithoutSpacing = thinnedDataWithoutSpacing->GetPoints();
	thinnedPointsType = (vtkStringArray*)thinnedDataWithoutSpacing->GetPointData()->GetAbstractArray("type");

	// check for need continue deletion or not
	bool needContinue = false;
	for (vtkIdType id = 0; id < thinnedPointsWithoutSpacing->GetNumberOfPoints(); id++){
		if (thinnedPointsType->GetValue(id) == "end"){

			oneEdgePoints4Detect = vtkSmartPointer<vtkPoints>::New();
			double aPoint[3]; thinnedPointsWithoutSpacing->GetPoint(id, aPoint);

			int endPoint[3] = { (int)aPoint[0], (int)aPoint[1], (int)aPoint[2] };
			bool av = false; ml3DThinning::DetectError(endPoint, av);
			if (av){
				needContinue = true;
				id = thinnedPointsWithoutSpacing->GetNumberOfPoints() + 10;
			}
		}
	}
	if (needContinue){
		ml3DThinning::ErrorDeleteProcess();
	}
}

void ml3DThinning::DetectErrorFromBifurcation(const int aPoint[3], int &lenght)
{
	if (!ExistedPointInSet(ml::vector3d(aPoint[0], aPoint[1], aPoint[2]), oneEdgePoints4Detect)){
		oneEdgePoints4Detect->InsertNextPoint(aPoint[0], aPoint[1], aPoint[2]);
		lenght++;
	}

	if (!ExistedPointInSet(ml::vector3d(aPoint[0], aPoint[1], aPoint[2]), detectedPoints)){
		detectedPoints->InsertNextPoint(aPoint[0], aPoint[1], aPoint[2]);
	}

	for (int z = aPoint[2] - 1; z <= aPoint[2] + 1; z++){
		for (int y = aPoint[1] - 1; y <= aPoint[1] + 1; y++){
			for (int x = aPoint[0] - 1; x <= aPoint[0] + 1; x++){

				if (x == aPoint[0] && y == aPoint[1] && z == aPoint[2]){ continue; }

				if (vesselArray[x + y*width + z*size2D] == 1) {
					if (!ExistedPointInSet(ml::vector3d(x, y, z), oneEdgePoints4Detect)
						&& !ExistedPointInSet(ml::vector3d(x, y, z), detectedPoints)
						&& oneEdgePoints4Detect->GetNumberOfPoints() < errorPaternLenght){

						int nextPoint[3] = { x, y, z };
						ml3DThinning::DetectErrorFromBifurcation(nextPoint, lenght);
					}
				}
			}
		}
	}
}

ml::points3d ml3DThinning::GetBifurcationPoints()
{

	return this->bifurcationPoints;
}

vtkSmartPointer<vtkPolyData> ml3DThinning::GetThinnedDataWithSpacing()
{
	return this->thinnedDataWithSpacing;
}

vtkSmartPointer<vtkPolyData> ml3DThinning::GetThinnedDataWithoutSpacing()
{
	return this->thinnedDataWithoutSpacing;
}

void ml3DThinning::GetDimensions(unsigned int &_w, unsigned int &_h, unsigned int &_sn)
{
	_w = width;
	_h = height;
	_sn = slice_num;
}

void ml3DThinning::GetSpacing(double &_xScale, double &_yScale, double &_zScale)
{
	_xScale = spacing[0];
	_yScale = spacing[1];
	_zScale = spacing[2];
}

unsigned int ml3DThinning::GetMaxVoxelValue()
{
	return this->maxVoxelValue;
}

unsigned char * ml3DThinning::GetThinBinaryArray()
{
	return vesselArray;
}

unsigned char * ml3DThinning::GetThinVessselArray()
{
	return bifurcationArray;
}

void ml3DThinning::Delete()
{
	if (!this->thinnedDataWithSpacing)
		return;
	if (this->thinnedDataWithSpacing->GetNumberOfPoints() == 0)
		return;

	delete[] vesselArray;
	delete[] bifurcationArray;
	this->bifurcationPoints.clear();
}