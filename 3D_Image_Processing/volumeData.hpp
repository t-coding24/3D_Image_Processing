#ifndef VOLUME_DATA_H
#define VOLUME_DATA_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;


class volumeData
{
public:
	volumeData();
	volumeData(string _volumeName, string _volumeMatrix);
	volumeData(string _volumeName, string _surfacePath, string _volumeMatrix);
	~volumeData();

	// Set functions
	void SetData(string _volumeName, string _surfacePath, string _volumeMatrix);
	void SetArray2(unsigned char * _array2);
	void SetArray3(unsigned char * _array3);
	void SetInputFiles(vector<InputFile> files);

	// Compute functions
	void GenerateBifurcationEndGraph();
	void GenerateVolumeData();
	void GenerateThinnedSurface();

	// Get functions
	vtkSmartPointer<vtkImageData> imageDataOfFile(string mhdPath);
	int* GetVolumeSize1D();
	double* GetVolumeSpacing();
	unsigned int* GetVolumeDataRGB();

	Matrix4x4 GetTranslationMatrix();
	vtkMutableDirectedGraph* GetThinnedGraph();
	vtkPolyData* GetThinnedSurface();
	vtkDataSetSurfaceFilter* GetVolumeSurface();
	vtkActor* VolumeSurfaceActor();

	// Correct graph
	void SaveThinDataToMHD();
	Vector3D TranslatedPoint(double beforeTranslatePositionX, double beforeTranslatePositionY, double beforeTranslatePositionZ, Matrix4x4 matrix);

	// Auto compute registration matrix between two input volume and matrix
	void AutoCombineMHD(vector<InputFile> files);
	void AutoCombineMHD2(vector<string> vtpFiles, vector<string> mtxFiles);




	//äÓëbãZèpââèKIIÇÃÇΩÇﬂÇ≈Ç∑
	vtkSmartPointer<vtkAppendPolyData> appendPoly;
	void StartComputeVesselRadius();
	void EdgeFromVertexID(vtkIdType fromId, vtkSmartPointer<vtkIdList> resultIds, ml::vector3d &foot, ml::vector3d &direction);
	//END

	//To calcrate angle between beam and vessel
	double GetVectorLength(Vector3D v);
	double GetInnerProduct(Vector3D v1, Vector3D v2);
	double GetAngleOf2Vector(Vector3D v1, Vector3D v2);
	void CreateLineFromBeamToVessel(double v1[3], double v2[3]);



private:
	Graph boostGraph; //Graph for detect avaible path from one vertex to another vertex
	vector<Edge> edges; //Edges array for boost graph
	vector<InputFile> inputFiles;

	// Parameters created when constructor was called 
	string		volumeName, surfacePath, volumeMatrix;
	vtkImageData	*volumeImg;
	Matrix4x4		translationMatrix;

	int				size1D[3];			//width, height, slice number
	int				size2D;				//size2D = width * height
	int				size3D;				//size3D = width * height * sliceNumber
	int				ext[6];
	double			bound[6];
	double			spacing[3];			//volume spacing

	// Parameters created when make full graph
	vtkSmartPointer<vtkMutableDirectedGraph>	graph;
	vtkSmartPointer<vtkPolyData>				thinnedSurface;				//Thinned data surface, use to display
	vtkSmartPointer<vtkPoints>					thinnedPoints;				//Thinned data points
	vtkSmartPointer<vtkUnsignedCharArray>		thinnedColour;				//Thinned data points color
	vtkSmartPointer<vtkStringArray>				thinnedType;				//Thinned data points type
	vtkSmartPointer<vtkPoints>					vertexsPoint;				//Array include vertex position
	vtkSmartPointer<vtkUnsignedCharArray>		vertexsColour;				//Array include vertex colour (red, green, yellow)
	vtkSmartPointer<vtkStringArray>				vertexsType;				//Array include vertex name (bifurcation, end, vessel)
	vtkSmartPointer<vtkModifiedBSPTree>			bspTree;					//Volume surface map, use to find cross point between line and surface

	bool										isSetArray2, isSetArray3;	//If true, we can start make graph function
	bool										needComputeVolumeSurface;
	unsigned char								*Array2, *Array3;			//The array include thinned data information

	Graph VTKGraphToBoostGraph();
	void InitParameters4Graph();

	Vector3D TranslatedPoint(double beforeTranslatePositionX, double beforeTranslatePositionY, double beforeTranslatePositionZ);
	Vector3D TranslatedPoint(double beforeTranslatePosition[3]){ return VolumeData::TranslatedPoint(beforeTranslatePosition[0], beforeTranslatePosition[1], beforeTranslatePosition[2]); };
	Vector3D TranslatedPoint(Vector3D beforeTranslatePosition){ return VolumeData::TranslatedPoint(beforeTranslatePosition[0], beforeTranslatePosition[1], beforeTranslatePosition[2]); };
	void addNoiseToPoint(Vector3D &aPoint);

	vtkSmartPointer<vtkActor> ShortPathBetWeen(Vector3D fromVertex, Vector3D toVertex, vtkMutableDirectedGraph* _graph);
	bool avaiblePathBetween(vtkIdType fromVertexID, vtkIdType toVertexID);

};

#endif