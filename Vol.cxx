#include <vtkActor.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkStructuredPoints.h>
#include <vtkDataSetMapper.h>
#include <vtkLookupTable.h>
#include <vtkThreshold.h>
#include <vtkShrinkFilter.h>
#include <vtkUnstructuredGrid.h>
int main(int, char*[])
{
  vtkNew<vtkNamedColors> colors;

  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renWin;
  renWin->AddRenderer(renderer);
  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(renWin);

  vtkNew<vtkStructuredPoints> vol;
  int sizeCube = 101;
  vol->SetDimensions(sizeCube, sizeCube, sizeCube);
  vol->SetOrigin(-0.5, -0.5, -0.5);
  auto sp = 1.0 / float(sizeCube - 1);
  vol->SetSpacing(sp, sp, sp);

  vtkNew<vtkDoubleArray> scalars;
  scalars->SetNumberOfComponents(1);
  scalars->SetNumberOfTuples(sizeCube * sizeCube * sizeCube);
  for (auto k = 0; k < sizeCube; k++)
  {
    auto z = k * sp;
    auto kOffset = k * sizeCube * sizeCube;
    for (auto j = 0; j < sizeCube; j++)
    {
      auto y = -0.5 + j * sp;
      auto jOffset = j * sizeCube;
      for (auto i = 0; i < sizeCube; i++)
      {
        auto x = -0.5 + i * sp;
        auto s = x * x + y * y + z * z - (0.4 * 0.4);
        auto offset = i + jOffset + kOffset;
        auto value = z - (sinf(i / 5.0) + cosf(j / 5.0)) / 2.0f;
        scalars->InsertTuple(offset, &value);
      }
    }
  }
  vol->GetPointData()->SetScalars(scalars);


  // Threshold
  vtkNew<vtkThreshold> threshold;
  threshold->SetInputData(vol);
  // Criterion is cells whose scalars are greater or equal to threshold.
  threshold->SetLowerThreshold(0.3);
  threshold->SetThresholdFunction(vtkThreshold::THRESHOLD_LOWER);
  threshold->Update();



  vtkNew<vtkDataSetMapper> volMapper;
  volMapper->SetInputConnection(threshold->GetOutputPort());
  volMapper->ScalarVisibilityOn();
  vtkNew<vtkActor> volActor;
  volActor->SetMapper(volMapper);
  volActor->GetProperty()->EdgeVisibilityOn();
  volActor->GetProperty()->SetColor(colors->GetColor3d("Salmon").GetData());
  renderer->AddActor(volActor);
  renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());
  renWin->SetSize(812, 812);
  renWin->SetWindowName("Vol");

  // interact with data
  renWin->Render();

  iren->Start();

  return EXIT_SUCCESS;
}