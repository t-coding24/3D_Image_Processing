#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

#include "VolumeData.h"
#include "Global.hpp"

#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>
#include <vtkObjectFactory.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkLink.h>
#include <vtkPropPicker.h>
#include <vtkAbstractPropPicker.h>
#include <vtkInteractorStyleTrackballActor.h>
#include <vtkActorCollection.h>
#include <vtkLineSource.h>
#include "vtkGraphGenerator.hpp"
#include "LoopDetection.h"
#include "vtkGraphFixer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkPointData.h"
#include <ml\mlVolumeFusion.hpp>

#include <Windows.h>
#include <wincon.h>


vtkSmartPointer<vtkImageData> VolumeData::imageDataOfFile(std::string mhdPath)
{
	vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(mhdPath.c_str());
	reader->Update();
	return reader->GetOutput();
}

VolumeData::VolumeData()
{
	printf("/ VolumeDataのデフォルトコンストラクタが呼ばれました\n");
}

VolumeData::VolumeData(std::string _volumeName, std::string _volumeMatrix)
{
	printf("/ VolumeDataのコンストラクタが呼ばれました\n");

	std::cout << _volumeName << std::endl;
	std::cout << _volumeMatrix << std::endl << std::endl;

	this->volumeName = _volumeName;
	this->volumeMatrix = _volumeMatrix;
	this->needComputeVolumeSurface = true;
	this->translationMatrix = ml::transform_matrix(this->volumeMatrix);

	vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(_volumeName.c_str());
	reader->Update();

	this->volumeImg = reader->GetOutput();
	this->volumeImg->GetBounds(this->bound);
	this->volumeImg->GetExtent(this->ext);
	this->volumeImg->GetSpacing(this->spacing);
	this->volumeImg->GetDimensions(this->size1D);
	this->size2D = this->size1D[0] * this->size1D[1];
	this->size3D = this->size1D[0] * this->size1D[1] * this->size1D[2];

	this->isSetArray2 = false;
	this->isSetArray3 = false;
}

VolumeData::VolumeData(std::string _volumeName, std::string _surfacePath, std::string _volumeMatrix)
{
	VolumeData::SetData(_volumeName, _surfacePath, _volumeMatrix);
}

void VolumeData::SetData(std::string _volumeName, std::string _surfacePath, std::string _volumeMatrix)
{
	printf("/ VolumeDataのコンストラクタが呼ばれました\n");

	// If volume avaible
	if (_volumeName == "nothing")
		return;
	else if (_volumeName != "nothing"){
		std::cout << _volumeName << std::endl;
		this->volumeName = _volumeName;

		vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
		reader->SetFileName(_volumeName.c_str());
		reader->Update();

		this->volumeImg = reader->GetOutput();
		this->volumeImg->GetBounds(this->bound);
		this->volumeImg->GetExtent(this->ext);
		this->volumeImg->GetSpacing(this->spacing);
		this->volumeImg->GetDimensions(this->size1D);
		this->size2D = this->size1D[0] * this->size1D[1];
		this->size3D = this->size1D[0] * this->size1D[1] * this->size1D[2];
		this->isSetArray2 = false;
		this->isSetArray3 = false;
	}

	// If surface avaible
	this->surfacePath = _surfacePath;
	if (_surfacePath != "nothing"){
		std::cout << _surfacePath << std::endl;
		this->needComputeVolumeSurface = false;
	}
	else if (_surfacePath == "nothing")
		this->needComputeVolumeSurface = true;


	// If matrix avaible
	this->volumeMatrix = _volumeMatrix;
	if (_volumeMatrix != "nothing"){
		std::cout << _volumeMatrix << std::endl << std::endl;
		this->translationMatrix = ml::transform_matrix(this->volumeMatrix);
	}
	else{
		printf("No input volume matrix\n");
	}

}

VolumeData::~VolumeData()
{
	printf("/ VolumeDataのデスコンストラクタが呼ばれました\n\n");
}

void VolumeData::SetArray2(unsigned char * _array2)
{
	Array2 = _array2;
	isSetArray2 = true;
}

void VolumeData::SetArray3(unsigned char * _array3)
{
	Array3 = _array3;
	isSetArray3 = true;
}

void VolumeData::SetInputFiles(std::vector<InputFile> files)
{
	this->inputFiles = files;
}

void VolumeData::GenerateThinnedSurface()
{
	if (isSetArray2 && isSetArray3)
	{
		this->thinnedSurface = vtkSmartPointer<vtkPolyData>::New();
		this->thinnedPoints = vtkSmartPointer<vtkPoints>::New();
		this->thinnedColour = vtkSmartPointer<vtkUnsignedCharArray>::New();
		this->thinnedColour->SetName("colors");
		this->thinnedColour->SetNumberOfComponents(3);
		this->thinnedType = vtkSmartPointer<vtkStringArray>::New();
		this->thinnedType->SetName("type");

		for (auto z = 0; z < this->size1D[2]; ++z){
			for (auto y = 0; y < this->size1D[1]; y++){
				for (auto x = 0; x < this->size1D[0]; ++x){

					if (this->Array2[(int)x + (int)y*this->size1D[0] + (int)z*this->size2D] == 1){

						double pp[3] = { x*this->spacing[0], y*this->spacing[1], z*this->spacing[2] };
						Vector3D translated = VolumeData::TranslatedPoint(pp[0], pp[1], pp[2]);
						this->thinnedPoints->InsertNextPoint(translated[0], translated[1], translated[2]);

						if (this->Array3[(int)x + (int)y*this->size1D[0] + (int)z*this->size2D] == 2){
							unsigned char col[3] = { 0, 255, 0 };
							this->thinnedColour->InsertNextTupleValue(col);
							this->thinnedType->InsertNextValue("end");
						}

						else if (this->Array3[(int)x + (int)y*this->size1D[0] + (int)z*this->size2D] == 1){
							unsigned char col[3] = { 255, 255, 0 };
							this->thinnedColour->InsertNextTupleValue(col);
							this->thinnedType->InsertNextValue("bifurcation");
						}

						else{ //vessel
							unsigned char col[3] = { 255, 0, 0 };
							this->thinnedColour->InsertNextTupleValue(col);
							this->thinnedType->InsertNextValue("vessel");
						}
					}
				}
			}
		}

		for (vtkIdType id = 0; id < this->thinnedPoints->GetNumberOfPoints(); id++)
		{
			double p[3]; this->thinnedPoints->GetPoint(id, p);
			p[0] += this->bound[0];
			p[1] += this->bound[2];
			p[2] += this->bound[4];
			this->thinnedPoints->SetPoint(id, p);
		}

		this->thinnedSurface->SetPoints(this->thinnedPoints);
		this->thinnedSurface->GetPointData()->SetScalars(this->thinnedColour);
		this->thinnedSurface->GetPointData()->AddArray(this->thinnedType);
		printf("/ GenerateThinnedSurface done\n");
	}
}

void VolumeData::GenerateVolumeData()
{
	if (isSetArray2 && isSetArray3)
	{
		// Compute spacing thinned surface
		VolumeData::GenerateThinnedSurface();

		// Compute spacing graph from thin array
		vtkGraphConverter converter;
		converter.SetDataInfos(Array2, Array3, this->inputFiles[0].volumePath);
		converter.SetNeedPlusBound(true);
		converter.Update();

		if (bvnSetting.fixGraph.on && !bvnSetting.fixGraph.loopDetectOn){
			vtkGraphFixer fixer;
			fixer.SetVolumeAndMatrix(this->inputFiles[0].volumePath, this->inputFiles[0].matrixPath);
			fixer.SetGraph(converter.GetResulGraph());
			fixer.SetFittingLineLenght(bvnSetting.fixGraph.fittingLineLenght);
			fixer.SetAngleOffset(bvnSetting.fixGraph.angleOffset);
			fixer.Update();

			this->graph = vtkSmartPointer<vtkMutableDirectedGraph>::New();
			this->graph->DeepCopy(fixer.GetResultGraph());
		}

		else if (!bvnSetting.fixGraph.on && bvnSetting.fixGraph.loopDetectOn){
			LoopDetection detection;
			detection.SetGraph(converter.GetResulGraph());
			detection.SetDeleteLoop(bvnSetting.fixGraph.loopDeleteOn);
			detection.Update();

			cout << "\n\nNumber of Loop : " << detection.GetNumberOfLoop() << endl;

			this->graph = vtkSmartPointer<vtkMutableDirectedGraph>::New();
			this->graph->DeepCopy(detection.GetResultGraph());
		}
		else if (bvnSetting.fixGraph.on && bvnSetting.fixGraph.loopDetectOn){
			vtkGraphFixer fixer;
			fixer.SetVolumeAndMatrix(this->inputFiles[0].volumePath, this->inputFiles[0].matrixPath);
			fixer.SetGraph(converter.GetResulGraph());
			fixer.SetFittingLineLenght(bvnSetting.fixGraph.fittingLineLenght);
			fixer.SetAngleOffset(bvnSetting.fixGraph.angleOffset);
			fixer.Update();

			LoopDetection detection;
			detection.SetGraph(fixer.GetResultGraph());
			detection.SetDeleteLoop(bvnSetting.fixGraph.loopDeleteOn);
			detection.Update();

			cout << "\n\nNumber of Loop : " << detection.GetNumberOfLoop() << endl;

			this->graph = vtkSmartPointer<vtkMutableDirectedGraph>::New();
			this->graph->DeepCopy(detection.GetResultGraph());
		}
		else {
			this->graph = vtkSmartPointer<vtkMutableDirectedGraph>::New();
			this->graph->DeepCopy(converter.GetResulGraph());
		}

		// 2014/11/11 Add for zemi
		StartComputeVesselRadius();
		// End

		printf("/ <FullGraph>の作成が完了しました\n");
	}
}

void VolumeData::StartComputeVesselRadius()
{
	vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(this->inputFiles[0].volumePath.c_str());
	reader->Update();

	vtkSmartPointer<vtkImageMarchingCubes> marching_cubes = vtkSmartPointer<vtkImageMarchingCubes>::New();
	marching_cubes->ComputeScalarsOff();
	marching_cubes->ComputeNormalsOn();
	marching_cubes->SetValue(0, bvnSetting.binarization);
	marching_cubes->SetInputData(reader->GetOutput());
	marching_cubes->Update();

	vtkSmartPointer<vtkDataSetSurfaceFilter> extractSurfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
	extractSurfaceFilter->SetInputConnection(marching_cubes->GetOutputPort());
	extractSurfaceFilter->Update();

	this->bspTree = vtkSmartPointer<vtkModifiedBSPTree>::New();
	this->bspTree->SetDataSet(extractSurfaceFilter->GetOutput());
	this->bspTree->BuildLocator();

	vtkSmartPointer<vtkGraphToPolyData> graphToPolyData = vtkSmartPointer<vtkGraphToPolyData>::New();
	graphToPolyData->SetInputData(this->graph);
	graphToPolyData->Update();

	this->vertexsPoint = graphToPolyData->GetOutput()->GetPoints();
	this->vertexsType = (vtkStringArray*)graphToPolyData->GetOutput()->GetPointData()->GetAbstractArray("type");
	this->vertexsColour = (vtkUnsignedCharArray*)graphToPolyData->GetOutput()->GetPointData()->GetScalars();

	this->appendPoly = vtkSmartPointer<vtkAppendPolyData>::New();
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	std::vector<vector3d> directionArray;
	std::vector<double> radiusArray;
	std::vector<double> vector_vessel_probeArray;
	ml::vector3d beamDirection;

	for (vtkIdType id = 0; id < this->vertexsType->GetNumberOfValues(); id++){
		if (this->vertexsType->GetValue(id) == "vessel")
		{
			const double *p = this->vertexsPoint->GetPoint(id);
			ml::vector3d foot(p[0], p[1], p[2]);
			ml::vector3d direction(0, 0, 0);
			vtkSmartPointer<vtkIdList> resultIds = vtkSmartPointer<vtkIdList>::New();
			EdgeFromVertexID(id, resultIds, foot, direction);
			if (!isZeroPoint(direction)){
				double radius = VesselRadius(direction, foot, this->bspTree);
				radiusArray.push_back(radius);
				directionArray.push_back(direction);
				vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
				polygonSource->SetNumberOfSides(20);
				polygonSource->SetCenter(foot[0], foot[1], foot[2]);
				polygonSource->SetRadius(radius);
				polygonSource->SetNormal(direction[0], direction[1], direction[2]);
				polygonSource->Update();
				this->appendPoly->AddInputConnection(polygonSource->GetOutputPort());
				this->appendPoly->Modified();
				this->appendPoly->Update();

				//2015_04_21_Yamashita_added


				ml::vector3d vesselDirection;
				std::vector<vector3d> vesselDirectionArray;
				vesselDirection[0] = direction[0];
				vesselDirection[1] = direction[1];
				vesselDirection[2] = direction[2];
				double _vesselDirection[3] = { static_cast<double>(vesselDirection[0], vesselDirection[1], vesselDirection[2]) };

				double _foot[3] = { static_cast<double>(foot[0], foot[1], foot[2]) };
				double _direction[3] = { static_cast<double>(direction[0], direction[1], direction[2]) };

				beamDirection[0] = direction[0];
				beamDirection[1] = 0;
				beamDirection[2] = direction[2];
				double _beamDirection[3] = { static_cast<double>(beamDirection[0], beamDirection[1], beamDirection[2]) };

				double degree_vessel_probe = GetAngleOf2Vector(vesselDirection, beamDirection);
				if (degree_vessel_probe > 90){
					degree_vessel_probe = degree_vessel_probe - 90;
				}
				vector_vessel_probeArray.push_back(degree_vessel_probe);

				vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
				lineSource->SetPoint1(_foot);
				lineSource->SetPoint2(_direction);

				// Visualize
				vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
				mapper->SetInputConnection(lineSource->GetOutputPort());
				actor->SetMapper(mapper);
				actor->GetProperty()->SetLineWidth(4);
				lineSource->Update();
			}
		}
	}


	//円描画
	textCombine.str(""); textCombine << volumeParentPath << "\\" << boost::filesystem::path(this->volumeName).stem().string() << "_RadiusSurface.vtp";
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(this->appendPoly->GetOutput());
	writer->SetFileName(textCombine.str().c_str());
	writer->SetCompressorTypeToNone();
	writer->Write();


	textCombine.str(""); textCombine << volumeParentPath << "\\" << boost::filesystem::path(this->volumeName).stem().string() << "_RadiusDatas.csv";
	std::ofstream log(textCombine.str(), ios::trunc);


	//Start compute real vesssel size
	double *doubleArray = new double[(int)radiusArray.size()];
	for (auto index = 0; index < radiusArray.size(); index++)
	{
		doubleArray[(int)index] = radiusArray[index] * radiusArray[index] * 3.1416;
	}
	std::sort(doubleArray, doubleArray + (int)radiusArray.size(), std::greater<double>());

	for (auto index = 0; index < radiusArray.size(); index++)
	{
		log << doubleArray[index] * 2 << "," << vector_vessel_probeArray[index] << "\n";
	}

	log.close();
}

void VolumeData::EdgeFromVertexID(vtkIdType fromId, vtkSmartPointer<vtkIdList> resultIds, ml::vector3d &foot, ml::vector3d &direction)
{
	if (resultIds->GetNumberOfIds() >= 7){
		Points3D points4Fitting;
		for (vtkIdType i = resultIds->GetNumberOfIds() - 1; i >= 0; i--){
			vtkIdType id = resultIds->GetId(i);
			const double *p = this->vertexsPoint->GetPoint(id);
			points4Fitting.push_back(Vector3D(p[0], p[1], p[2]));
		}
		ml::line lineFit;
		lineFit.fit(points4Fitting);
		foot = lineFit.perpendicular_foot(foot);
		direction = lineFit.direction();
		return;
	}

	resultIds->InsertUniqueId(fromId);

	if (resultIds->GetNumberOfIds() >= 7){
		Points3D points4Fitting;
		for (vtkIdType i = resultIds->GetNumberOfIds() - 1; i >= 0; i--){
			vtkIdType id = resultIds->GetId(i);
			const double *p = this->vertexsPoint->GetPoint(id);
			points4Fitting.push_back(Vector3D(p[0], p[1], p[2]));
		}
		ml::line lineFit;
		lineFit.fit(points4Fitting);
		foot = lineFit.perpendicular_foot(foot);
		direction = lineFit.direction();
		return;
	}

	vtkSmartPointer<vtkInEdgeIterator> inIt = vtkSmartPointer<vtkInEdgeIterator>::New();
	this->graph->GetInEdges(fromId, inIt);

	vtkSmartPointer<vtkOutEdgeIterator> outIt = vtkSmartPointer<vtkOutEdgeIterator>::New();
	this->graph->GetOutEdges(fromId, outIt);

	bool hasEndOrBifur = false;
	vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
	while (inIt->HasNext()){
		vtkInEdgeType type = inIt->Next();
		if (type.Source == fromId)
			continue;

		if (!IsAvaibleIDInIdList(resultIds, type.Source)){
			idList->InsertUniqueId(type.Source);
			if (this->vertexsType->GetValue(type.Source) == "bifurcation" || this->vertexsType->GetValue(type.Source) == "end"){
				hasEndOrBifur = true;
			}
		}
	}
	while (outIt->HasNext()){
		vtkOutEdgeType type = outIt->Next();
		if (type.Target == fromId)
			continue;

		if (!IsAvaibleIDInIdList(resultIds, type.Target)){
			idList->InsertUniqueId(type.Target);
			if (this->vertexsType->GetValue(type.Target) == "bifurcation" || this->vertexsType->GetValue(type.Target) == "end"){
				hasEndOrBifur = true;
			}
		}
	}

	if (hasEndOrBifur){
		if (resultIds->GetNumberOfIds() >= 3){
			Points3D points4Fitting;
			for (vtkIdType i = resultIds->GetNumberOfIds() - 1; i >= 0; i--){
				vtkIdType id = resultIds->GetId(i);
				const double *p = this->vertexsPoint->GetPoint(id);
				points4Fitting.push_back(Vector3D(p[0], p[1], p[2]));
			}
			ml::line lineFit;
			lineFit.fit(points4Fitting);
			foot = lineFit.perpendicular_foot(foot);
			direction = lineFit.direction();
		}
		return;
	}
	else{
		for (vtkIdType id = 0; id < idList->GetNumberOfIds(); id++)
			EdgeFromVertexID(idList->GetId(id), resultIds, foot, direction);
	}
}

void VolumeData::GenerateBifurcationEndGraph()
{
	if (isSetArray2 && isSetArray3)
	{
		// Compute spacing thinned surface
		VolumeData::GenerateThinnedSurface();

		// Compute spacing graph from thin array
		vtkGraphConverter converter;
		converter.SetDataInfos(Array2, Array3, this->inputFiles[0].volumePath);
		converter.SetNeedPlusBound(true);
		converter.SetGraphMode(BIFURCATION_END_GRAPH);
		converter.Update();

		this->graph = vtkSmartPointer<vtkMutableDirectedGraph>::New();
		this->graph->DeepCopy(converter.GetResulGraph());

		printf("/ <FullGraph>の作成は完了でした\n");
	}
}

int* VolumeData::GetVolumeSize1D()
{
	return this->size1D;
}

double* VolumeData::GetVolumeSpacing()
{
	return this->spacing;
}

unsigned int* VolumeData::GetVolumeDataRGB()
{
	unsigned int *volumeDataRGB = new unsigned int[this->size3D];
	const double *scalar_range = this->volumeImg->GetScalarRange();

	for (auto z = this->ext[4]; z < this->ext[5]; ++z){
		for (auto y = this->ext[2]; y < this->ext[3]; ++y){
			for (auto x = this->ext[0]; x < this->ext[1]; ++x){
				if (this->volumeImg->GetScalarType() == VTK_UNSIGNED_CHAR){
					unsigned char* voxel = static_cast<unsigned char*>(this->volumeImg->GetScalarPointer(x, y, z));
					if (voxel[0]>0){
						int _x = x - this->ext[0];
						int _y = y - this->ext[2];
						int _z = z - this->ext[4];
						volumeDataRGB[_x + _y*this->size1D[0] + _z*this->size2D] = (unsigned int)voxel[0];
					}
				}
				else{
					unsigned short* voxel = static_cast<unsigned short*>(this->volumeImg->GetScalarPointer(x, y, z));
					if (*voxel > 0){
						int _x = x - this->ext[0];
						int _y = y - this->ext[2];
						int _z = z - this->ext[4];
						unsigned int unsigned_char_value = 255 * (*voxel) / scalar_range[1];
						volumeDataRGB[_x + _y*this->size1D[0] + _z*this->size2D] = unsigned_char_value;
					}
				}
			}
		}
	}

	return volumeDataRGB;
}

Matrix4x4 VolumeData::GetTranslationMatrix()
{
	return this->translationMatrix;
}

vtkMutableDirectedGraph* VolumeData::GetThinnedGraph()
{
	return this->graph;
}

vtkPolyData* VolumeData::GetThinnedSurface()
{
	return this->thinnedSurface;
}

void VolumeData::InitParameters4Graph()
{
	this->graph = vtkSmartPointer<vtkMutableDirectedGraph>::New();
	this->thinnedSurface = vtkSmartPointer<vtkPolyData>::New();
	this->thinnedPoints = vtkSmartPointer<vtkPoints>::New();
	this->thinnedColour = vtkSmartPointer<vtkUnsignedCharArray>::New();
	this->thinnedColour->SetName("colors");
	this->thinnedColour->SetNumberOfComponents(3);
	this->thinnedType = vtkSmartPointer<vtkStringArray>::New();
	this->thinnedType->SetName("type");
	printf("/ InitParameters4Graphが終了しました\n");
}

Graph VolumeData::VTKGraphToBoostGraph()
{
	for (vtkIdType e = 0; e < this->graph->GetNumberOfEdges(); e++)
	{
		vtkIdType fromVertexID = this->graph->GetSourceVertex(e);
		vtkIdType toVertexID = this->graph->GetTargetVertex(e);

		Edge aEdge = std::make_pair((int)fromVertexID, (int)toVertexID);
		edges.push_back(aEdge);
	}
	Graph g(edges.begin(), edges.end(), (int)this->graph->GetNumberOfVertices());

	//std::ofstream file("test.dot");
	//boost::write_graphviz(file, g);//, boost::make_label_writer("testGraph"));

	return g;
}

Vector3D VolumeData::TranslatedPoint(double beforeTranslatePositionX, double beforeTranslatePositionY, double beforeTranslatePositionZ)
{
	Vector3D returnPoint(beforeTranslatePositionX, beforeTranslatePositionY, beforeTranslatePositionZ);
	if (this->volumeMatrix != "nothing"){
		returnPoint = this->translationMatrix * returnPoint;
	}
	return returnPoint;
}

Vector3D VolumeData::TranslatedPoint(double beforeTranslatePositionX, double beforeTranslatePositionY, double beforeTranslatePositionZ, Matrix4x4 matrix)
{
	Vector3D returnPoint(beforeTranslatePositionX, beforeTranslatePositionY, beforeTranslatePositionZ);
	returnPoint = matrix * returnPoint;
	return returnPoint;
}

bool VolumeData::avaiblePathBetween(vtkIdType fromVertexID, vtkIdType toVertexID)
{
	std::vector<boost::default_color_type> color(this->graph->GetNumberOfVertices(), boost::white_color);
	return boost::is_reachable(fromVertexID, toVertexID, boostGraph, color.data());
}

vtkSmartPointer<vtkActor> VolumeData::ShortPathBetWeen(Vector3D fromVertex, Vector3D toVertex, vtkMutableDirectedGraph* _graph)
{
	vtkSmartPointer<vtkActor> pathActor = vtkSmartPointer<vtkActor>::New();
	vtkSmartPointer<vtkPoints> vPoints = _graph->GetPoints();
	vtkIdType fromVertexID, toVertexID;

	for (vtkIdType id = 0; id < vPoints->GetNumberOfPoints(); id++){
		double *p = vPoints->GetPoint(id);
		if (p[0] == fromVertex[0] && p[1] == fromVertex[1] && p[2] == fromVertex[2])
			fromVertexID = id;
		else if (p[0] == toVertex[0] && p[1] == toVertex[1] && p[2] == toVertex[2])
			toVertexID = id;
	}

	if (fromVertexID >= 0 && toVertexID >= 0){
		if (fromVertexID < vPoints->GetNumberOfPoints() && toVertexID < vPoints->GetNumberOfPoints()){
			// Convert the graph to a polydata
			vtkSmartPointer<vtkGraphToPolyData> graphToPolyData = vtkSmartPointer<vtkGraphToPolyData>::New();
			graphToPolyData->SetInputData(_graph);
			graphToPolyData->Update();

			vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra = vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
			dijkstra->SetInputConnection(graphToPolyData->GetOutputPort());
			dijkstra->SetStartVertex(fromVertexID);
			dijkstra->SetEndVertex(toVertexID);
			dijkstra->Update();

			vtkSmartPointer<vtkVertexGlyphFilter> vertexfilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
			vertexfilter->SetInputConnection(dijkstra->GetOutputPort());
			vertexfilter->Update();

			//printf("number of points %d\n",vertexfilter->GetOutput()->GetPoints()->GetNumberOfPoints());

			//printf("in\n");
			std::cout << vertexfilter->GetOutput()->GetPoints()->GetNumberOfPoints() << std::endl;

			//// Create a mapper and actor
			//vtkSmartPointer<vtkPolyDataMapper> pathMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			//pathMapper->SetInputConnection(dijkstra->GetOutputPort());
			//
			//pathActor->SetMapper(pathMapper);
			//pathActor->GetProperty()->SetColor(1,0,0); // Red
			//pathActor->GetProperty()->SetLineWidth(4);
		}
	}

	return pathActor;
}

void VolumeData::SaveThinDataToMHD()
{
	double min[3] = { 10000, 10000, 10000 }, max[3] = { -10000, -10000, -10000 };
	for (vtkIdType id = 0; id < this->thinnedPoints->GetNumberOfPoints(); id++)
	{
		double p[3]; this->thinnedPoints->GetPoint(id, p);
		for (int idx = 0; idx < 3; ++idx)
			min[idx] = min[idx] < p[idx] ? min[idx] : p[idx];
		for (int idx = 0; idx<3; ++idx)
			max[idx] = max[idx] > p[idx] ? max[idx] : p[idx];
	}

	int minInINT[3], maxInINT[3];
	for (int idx = 0; idx < 3; ++idx){
		minInINT[idx] = (min[idx] / this->spacing[idx]) - 2;
		maxInINT[idx] = (max[idx] / this->spacing[idx]) + 2;
	}

	vtkSmartPointer<vtkImageData> thinMHD = vtkSmartPointer<vtkImageData>::New();
	thinMHD->SetSpacing(this->spacing);
	thinMHD->SetExtent(minInINT[0], maxInINT[0], minInINT[1], maxInINT[1], minInINT[2], maxInINT[2]);
	thinMHD->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
	for (vtkIdType id = 0; id < this->thinnedPoints->GetNumberOfPoints(); id++)
	{
		double p[3]; this->thinnedPoints->GetPoint(id, p);
		unsigned char* volume_scalar = reinterpret_cast<unsigned char*>(thinMHD->GetScalarPointer(p[0] / this->spacing[0], p[1] / this->spacing[1], p[2] / this->spacing[2]));
		*volume_scalar = 255;
	}
	thinMHD->Modified();

	/*vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();
	writer->SetInput(thinMHD);
	writer->SetFileName("test.mhd");
	writer->Write();*/
}

void VolumeData::AutoCombineMHD(std::vector<InputFile> files)
{
	std::vector<std::string> _volumesPaths;
	std::vector<vtkSmartPointer<vtkPoints>> _thinBifurPoints;
	//ml::transform_matrix file1Mtx;

	for (int idx = 0; idx < files.size(); ++idx){
		_volumesPaths.push_back(files[idx].volumePath);
		ml3DThinning thin = ThinOnFile(files[idx].volumePath);

		vtkSmartPointer<vtkImageData> image = imageDataOfFile(files[idx].volumePath);
		double ori[3]; image->GetOrigin(ori);
		ml::transform_matrix matrixxx = ml::transform_matrix(files[idx].matrixPath);

		//if (idx == 0){
		//	file1Mtx = matrixxx;
		//}
		matrixxx.translate_local(ori[0], ori[1], ori[2]);

		vtkGraphConverter converter;
		converter.SetDataInfos(thin.GetThinBinaryArray(), thin.GetThinVessselArray(), files[idx].volumePath);
		converter.SetTranformMatrix(matrixxx);
		converter.Update();

		if (bvnSetting.fixGraph.loopDetectOn && bvnSetting.fixGraph.loopDeleteOn){
			cout << "Error loop delete start..." << endl;

			LoopDetection detection;
			detection.SetGraph(converter.GetResulGraph());
			detection.SetDeleteLoop(bvnSetting.fixGraph.loopDeleteOn);
			detection.Update();

			_thinBifurPoints.push_back(detection.GetVtkBifurcationPoints());
		}
		else {
			_thinBifurPoints.push_back(converter.GetVtkBifurcationPoints());
		}
	}

	// 2014/12/26 Mod
	std::vector<std::string> volumesPaths;
	std::vector<vtkSmartPointer<vtkPoints>> thinBifurPoints;
	ml::transform_matrix file1Mtx;

	if (_thinBifurPoints[0]->GetNumberOfPoints() < _thinBifurPoints[1]->GetNumberOfPoints()){
		volumesPaths = _volumesPaths;
		thinBifurPoints = _thinBifurPoints;
		file1Mtx = ml::transform_matrix(files[0].matrixPath);
	}
	else {
		for (size_t idx = _volumesPaths.size() - 1; idx >= 0; --idx){
			if (idx < 0 || idx >= _volumesPaths.size())
				break;

			volumesPaths.push_back(_volumesPaths[idx]);
			thinBifurPoints.push_back(_thinBifurPoints[idx]);
		} file1Mtx = ml::transform_matrix(files[1].matrixPath);
	}
	// End mod

	ml::auto_registration registration;
	registration.SetFile1Matrix(file1Mtx);
	registration.Set2OriginVolumePath(volumesPaths);
	registration.Set2vtkBifurcationPoints(thinBifurPoints);
	registration.SetMaxCommonBifurcationPoint(bvnSetting.fusionSettinng.maxCommonBifurcationPoint);
	registration.SetMinRMSE4StopProcess(bvnSetting.fusionSettinng.minRMSE4StopRegistration);
	registration.SetMinDisBetween2PointsIn1Set(bvnSetting.fusionSettinng.minDisBetween2PointIn1Set);
	registration.SetMaxDisBetween2VolumeOnWorldCoordinate(bvnSetting.fusionSettinng.maxDistaceBetween2Volume);
	registration.SetDis4Check2SetPointIsSameShape(bvnSetting.fusionSettinng.dis4Check2SetPointIsoform);
	registration.SetMaxRegistrationMatrixAngle(bvnSetting.fusionSettinng.maxRegistrationMatrixAngle);
	registration.SetRegistrationProcessMode(bvnSetting.fusionSettinng.registrationMode);
	registration.StartCompute();
}

void VolumeData::AutoCombineMHD2(std::vector<std::string> vtpFiles, std::vector<std::string> mtxFiles)
{
	std::vector<std::string> volumesPaths;
	std::vector<vtkSmartPointer<vtkPoints>> thinBifurPoints;
	ml::transform_matrix file1Mtx;

	for (int idx = 0; idx < vtpFiles.size(); ++idx){
		if (idx == 0){
			file1Mtx = ml::transform_matrix(mtxFiles[idx]);
		}

		ml::transform_matrix matrix(mtxFiles[idx]);

		vtkSmartPointer<vtkXMLPolyDataReader> reader = VTK_NEW(vtkXMLPolyDataReader);
		reader->SetFileName(vtpFiles[idx].c_str());
		reader->Update();
		vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();

		vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
		vtkSmartPointer<vtkStringArray> pointsType = (vtkStringArray*)polyData->GetPointData()->GetAbstractArray("type");
		vtkSmartPointer<vtkPoints> bifurPoints = vtkSmartPointer<vtkPoints>::New();

		for (vtkIdType id = 0; id < pointsType->GetNumberOfValues(); id++)
		{
			vtkStdString value = pointsType->GetValue(id);
			if (value == "bifurcation"){
				const double *p = points->GetPoint(id);
				ml::vector3d transform_p = matrix * ml::vector3d(p[0], p[1], p[2]);
				if (idx == 1){
					VolumeData::addNoiseToPoint(transform_p);
				}
				bifurPoints->InsertNextPoint(transform_p[0], transform_p[1], transform_p[2]);
			}
		}

		volumesPaths.push_back(vtpFiles[idx]);
		thinBifurPoints.push_back(bifurPoints);
	}

	ml::auto_registration registration;
	registration.SetFile1Matrix(file1Mtx);
	registration.Set2OriginVolumePath(volumesPaths);
	registration.Set2vtkBifurcationPoints(thinBifurPoints);
	registration.SetMaxCommonBifurcationPoint(bvnSetting.fusionSettinng.maxCommonBifurcationPoint);
	registration.SetMinRMSE4StopProcess(bvnSetting.fusionSettinng.minRMSE4StopRegistration);
	registration.SetMinDisBetween2PointsIn1Set(bvnSetting.fusionSettinng.minDisBetween2PointIn1Set);
	registration.SetMaxDisBetween2VolumeOnWorldCoordinate(bvnSetting.fusionSettinng.maxDistaceBetween2Volume);
	registration.SetDis4Check2SetPointIsSameShape(bvnSetting.fusionSettinng.dis4Check2SetPointIsoform);
	registration.SetMaxRegistrationMatrixAngle(bvnSetting.fusionSettinng.maxRegistrationMatrixAngle);
	registration.SetRegistrationProcessMode(bvnSetting.fusionSettinng.registrationMode);
	registration.StartCompute();
}

void VolumeData::addNoiseToPoint(Vector3D &aPoint)
{
	for (size_t idx = 0; idx < 2; idx++){
		aPoint[idx] += ((double)rand() / (RAND_MAX + 1)) * 2 - 1;
	}
}


