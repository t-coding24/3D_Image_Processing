









double VolumeData::GetVectorLength(Vector3D v)
{
	return pow((v[0] * v[0]) + (v[0] * v[0]) + (v[0] * v[0]), 0.5);
}

double VolumeData::GetInnerProduct(Vector3D v1, Vector3D v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

double VolumeData::GetAngleOf2Vector(Vector3D v1, Vector3D v2)
{
	double length_v1 = GetVectorLength(v1);
	double length_v2 = GetVectorLength(v2);

	double cos_sita = GetInnerProduct(v1, v2) / (length_v1 * length_v2);

	double sita = acos(cos_sita);

	sita = sita * 180.0 / 3.1416;

	return sita;
}

void VolumeData::CreateBeamLineToVessel(double v1[3], double v2[3])
{

	vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
	lineSource->SetPoint1(v1);
	lineSource->SetPoint2(v2);
	lineSource->Update();

	// Visualize
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(lineSource->GetOutputPort());
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetLineWidth(4);

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);

	renderWindow->Render();
	renderWindowInteractor->Start();

}