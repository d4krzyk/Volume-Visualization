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
#include <vtkColorTransferFunction.h>
#include <QApplication>
#include <QSlider>
#include <QVBoxLayout>
#include <QWidget>
#include <QMainWindow>
#include <QDockWidget>
#include <QPointer>
#include <QPushButton>
#include <QLabel>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <QApplication>
#include <QDockWidget>
#include <QGridLayout>
#include <QLabel>
#include <QMainWindow>
#include <QPointer>
#include <QPushButton>
#include <QVBoxLayout>
#include <vtkCamera.h>
#include <qtimer.h>
#include <vtkCellPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkWarpScalar.h>
#include <vtkGeometryFilter.h>

class VolumeRenderer {
public:
    VolumeRenderer(vtkRenderer* renderer, vtkStructuredPoints* vol )
        : renderer(renderer), vol(vol), frequencyX(10.0f), frequencyY(10.0f), amplitude(10.0f) {}

    bool isPlaying = false;

    void updateVolume() {
        vtkNew<vtkDoubleArray> scalars;
        int sizeCube = vol->GetDimensions()[0];
        auto sp = vol->GetSpacing()[0];
        scalars->SetNumberOfComponents(1);
        scalars->SetNumberOfTuples(sizeCube * sizeCube * sizeCube);

        for (auto k = 0; k < sizeCube; k++) {
            auto z = k * sp;
            auto kOffset = k * sizeCube * sizeCube;
            for (auto j = 0; j < sizeCube; j++) {
                auto y = -0.5 + j * sp;
                auto jOffset = j * sizeCube;
                for (auto i = 0; i < sizeCube; i++) {
                    auto x = -0.5 + i * sp; 
                    auto value = z - ( sinf(i / ((sizeCube - 1) / frequencyX) + time) * cosf(j / ((sizeCube - 1) / frequencyY) + time)) / ((sizeCube - 1) / amplitude) ;
                    auto offset = i + jOffset + kOffset;
                    scalars->InsertTuple(offset, &value);
                }
            }
        }
        vol->GetPointData()->SetScalars(scalars);
        renderer->GetRenderWindow()->Render();
        
    }

    void setFrequencyX(double freqx) {
        frequencyX = freqx;
        updateVolume();
    }
    void setFrequencyY(double freqy) {
        frequencyY = freqy;
        updateVolume();
    }
    void setAmplitude(double amp) {
        amplitude = amp;
        updateVolume();
    }
	void setSizeCube(int size) {
        vol->SetDimensions(size, size, size);
        auto sp = 1.0 / float(size - 1);
        vol->SetSpacing(sp, sp, sp);
		updateVolume();
	}
    void updateTime(double deltaTime) {
        time += deltaTime; // Zwiększanie czasu
        updateVolume();
    }
    void pickPoint(int x, int y) {
        vtkNew<vtkCellPicker> picker;
        picker->SetTolerance(0.005);

        if (picker->Pick(x, y, 0, renderer)) {
            double pickedPosition[3];
            picker->GetPickPosition(pickedPosition);

            vtkIdType pointId = vol->FindPoint(pickedPosition);
            if (pointId != -1) {
                std::cout << "Picked Point ID: " << pointId << std::endl;
                std::cout << "Picked Position: (" << pickedPosition[0] << ", "
                    << pickedPosition[1] << ", " << pickedPosition[2] << ")" << std::endl;

                //changePointScalar(pointId, 100.0); // Nowa wartość skalarna
                increasePointHeight(pointId, -10.0);
                //updateVoxels();
                //applyWarpScalar(); // Zastosuj deformację
                // Pobierz wartość skalaru dla wybranego punktu
                vtkDataArray* scalars = vol->GetPointData()->GetScalars();
                if (scalars) {
                    double value = scalars->GetTuple1(pointId);
                    std::cout << "Scalar Value at Picked Point: " << value << std::endl;
                }
            }
            else {
                std::cout << "No valid point found at picked position." << std::endl;
            }
        }
        else {
            std::cout << "No point picked." << std::endl;
        }
    }

    void changePointScalar(vtkIdType pointId, double newValue) {
        vtkDataArray* scalars = vol->GetPointData()->GetScalars();
        if (scalars && pointId != -1) {
            scalars->SetTuple1(pointId, newValue); // Zmiana wartości skalarnej
            vol->Modified(); 
            renderer->GetRenderWindow()->Render(); 
        }
    }
    void updateVoxels() {
        vtkNew<vtkThreshold> threshold;
        threshold->SetInputData(vol);
        threshold->SetLowerThreshold(0.55);
        threshold->SetThresholdFunction(vtkThreshold::THRESHOLD_LOWER);
        threshold->Update();

        vtkNew<vtkDataSetMapper> mapper;
        mapper->SetInputConnection(threshold->GetOutputPort());
        mapper->ScalarVisibilityOn();
        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(1.0, 1.0, 1.0); // Ustaw kolor biały (lub dynamiczny)
        actor->GetProperty()->EdgeVisibilityOn();
        actor->GetProperty()->SetInterpolationToFlat();
        actor->GetProperty()->SetAmbient(1.0);  // Ustawienie maksymalnej jasności materiału (1.0 = maksymalna jasność)
        actor->GetProperty()->SetDiffuse(0.2);  // Zmniejszenie odbicia rozproszonego
        actor->GetProperty()->SetSpecular(0.1);  // Można dodać połysk (dla efektu refleksji)
        actor->GetProperty()->SetSpecularPower(1.0);  // Zwiększenie refleksów świetlnych
        renderer->RemoveAllViewProps(); // Usuń poprzednie
        renderer->AddActor(actor);
        renderer->GetRenderWindow()->Render();
    }

    void increasePointHeight(vtkIdType pointId, double increment) {
        vtkDataArray* scalars = vol->GetPointData()->GetScalars();
        if (scalars && pointId != -1) {
            double currentValue = scalars->GetTuple1(pointId);
            double newValue = currentValue + increment; // Zwiększenie wysokości
            scalars->SetTuple1(pointId, newValue);      
            vol->Modified();                            
            renderer->GetRenderWindow()->Render();      
        }
    }
    void applyWarpScalar() {
        vtkNew<vtkWarpScalar> warp;
        warp->SetInputData(vol);
        warp->SetScaleFactor(0.0); // Skala deformacji
        warp->Update();

        vtkNew<vtkDataSetMapper> mapper;
        mapper->SetInputConnection(warp->GetOutputPort());

        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);

        renderer->RemoveAllViewProps(); // Usuń poprzednie obiekty
        renderer->AddActor(actor);
        renderer->GetRenderWindow()->Render();
    }

private:
    vtkRenderer* renderer;
    vtkStructuredPoints* vol;
    double frequencyX;
    double frequencyY;
    double amplitude;
	
    double time; // Czas, który będzie animował falę
};
void rotateCameraRight(vtkNew<vtkCamera>& camera, vtkRenderer* renderer) {
    double* position = camera->GetPosition();
    double* focalPoint = camera->GetFocalPoint();
    camera->Azimuth(15);  // Obrót kamery o 10 stopni w prawo
    camera->SetViewUp(0.0, 0.0, 1.0);
    camera->SetFocalPoint(0.0, 0.0, 0.0);
    camera->SetPosition(position);
    renderer->GetActiveCamera()->Modified();  // Aktualizowanie widoku
}
// Funkcja obracania kamery w lewo
void rotateCameraLeft(vtkNew<vtkCamera>& camera, vtkRenderer* renderer) {
    double* position = camera->GetPosition();
    double* focalPoint = camera->GetFocalPoint();
    camera->Azimuth(-15);  // Obrót kamery o 10 stopni w lewo
    camera->SetViewUp(0.0, 0.0, 1.0);
    camera->SetFocalPoint(0.0, 0.0, 0.0);
    camera->SetPosition(position);
    renderer->GetActiveCamera()->Modified();  // Aktualizowanie widoku
}

class MouseInteractorStyle : public vtkInteractorStyleTrackballCamera {
public:
    static MouseInteractorStyle* New();
    vtkTypeMacro(MouseInteractorStyle, vtkInteractorStyleTrackballCamera);

    void OnRightButtonDown() override {
        int* clickPos = this->GetInteractor()->GetEventPosition();

        if (volumeRenderer) {
            volumeRenderer->pickPoint(clickPos[0], clickPos[1]);
        }

        vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }

    void setVolumeRenderer(VolumeRenderer* renderer) {
        volumeRenderer = renderer;
    }

private:
    VolumeRenderer* volumeRenderer = nullptr;
};
vtkStandardNewMacro(MouseInteractorStyle);
int main(int argc, char* argv[])
{
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());

    QApplication app(argc, argv);

    // main window
    QMainWindow mainWindow;
    mainWindow.resize(900, 900);
    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkDataSetMapper> volMapper;
    vtkNew<vtkRenderer> renderer;

    vtkNew<vtkCamera> camera;
    camera->SetPosition(2.0, -2.0, 2.0);  // Ustawienie kamery na pozycji (0, -5, 5)
    camera->SetFocalPoint(0.0, 0.0, 0.0); // Skierowanie kamery na środek sceny (0, 0, 0)
    camera->SetViewUp(0.0, 0.0, 1.0); // Określenie, który kierunek jest "górą"
    renderer->SetActiveCamera(camera);

    // Parameters
    QPushButton RotateCameraRight("Rotate Camera Right");
    QPushButton RotateCameraLeft("Rotate Camera Left");
    RotateCameraLeft.setMinimumWidth(50);
    RotateCameraRight.setMinimumWidth(50);
    RotateCameraLeft.setMaximumWidth(150);
    RotateCameraRight.setMaximumWidth(150);
	// Time offset button
    QPushButton playStopButton("Play/Stop Phase");

    QHBoxLayout* buttonLayout = new QHBoxLayout();
    buttonLayout->addWidget(&RotateCameraLeft);
    buttonLayout->addWidget(&RotateCameraRight);
    buttonLayout->addWidget(&playStopButton);
    QWidget layoutContainer;

    QVBoxLayout* sliderLayout = new QVBoxLayout();

    QSlider frequencyXSlider(Qt::Horizontal);
    frequencyXSlider.setRange(-100, 100);
    frequencyXSlider.setValue(0);

    QSlider frequencyYSlider(Qt::Horizontal);
    frequencyYSlider.setRange(-100, 100);
    frequencyYSlider.setValue(0);

    QSlider amplitudeSlider(Qt::Horizontal);
    amplitudeSlider.setRange(-100, 100);
    amplitudeSlider.setValue(0);

    QSlider sizeCubeSlider(Qt::Horizontal);
    sizeCubeSlider.setRange(3, 101);
    sizeCubeSlider.setValue(51);




    QLabel label_freq_x;
    label_freq_x.setText("Frequency X");
    sliderLayout->addWidget(&label_freq_x);
    sliderLayout->addWidget(&frequencyXSlider);
    QLabel label_freq_y;
    label_freq_y.setText("Frequency Y");
    sliderLayout->addWidget(&label_freq_y);
    sliderLayout->addWidget(&frequencyYSlider);
    QLabel label_Amp;
    label_Amp.setText("Amplitude");
    sliderLayout->addWidget(&label_Amp);
    sliderLayout->addWidget(&amplitudeSlider);
    QLabel label_SizeCube;
    label_SizeCube.setText("Cube Size");
    sliderLayout->addWidget(&label_SizeCube);
    sliderLayout->addWidget(&sizeCubeSlider);
    label_freq_x.setAlignment(Qt::AlignCenter);
    label_freq_y.setAlignment(Qt::AlignCenter);
    label_Amp.setAlignment(Qt::AlignCenter);
    label_SizeCube.setAlignment(Qt::AlignCenter);

    // control area
    QDockWidget controlDock;
    mainWindow.addDockWidget(Qt::BottomDockWidgetArea, &controlDock);

    QLabel controlDockTitle("Control Deck");
    controlDockTitle.setAlignment(Qt::AlignCenter);
    controlDockTitle.setMargin(5);
    controlDock.setTitleBarWidget(&controlDockTitle);

    QPointer<QVBoxLayout> dockLayout = new QVBoxLayout();

    dockLayout->addLayout(buttonLayout);
    dockLayout->addLayout(sliderLayout);

    layoutContainer.setLayout(dockLayout);
    controlDock.setWidget(&layoutContainer);
    // render area
    QPointer<QVTKOpenGLNativeWidget> vtkRenderWidget = new QVTKOpenGLNativeWidget();
    mainWindow.setCentralWidget(vtkRenderWidget);



    // VTK part
    vtkNew<vtkGenericOpenGLRenderWindow> window;
    vtkRenderWidget->setRenderWindow(window.Get());

    vtkNew<vtkStructuredPoints> vol;
    int sizeCube = 51;
    vol->SetDimensions(sizeCube, sizeCube, sizeCube);
    vol->SetOrigin(-0.5, -0.5, -0.5);
    auto sp = 1.0 / float(sizeCube - 1) ;
    vol->SetSpacing(sp, sp, sp);

    vtkNew<vtkDoubleArray> scalars;
    scalars->SetNumberOfComponents(1);
    scalars->SetNumberOfTuples(sizeCube * sizeCube * sizeCube);
    
    vol->GetPointData()->SetScalars(scalars);

    // Threshold
    vtkNew<vtkThreshold> threshold;
    threshold->SetInputData(vol);
    // Criterion is cells whose scalars are greater or equal to threshold.
    threshold->SetLowerThreshold(0.55);
    threshold->SetThresholdFunction(vtkThreshold::THRESHOLD_LOWER);
    threshold->Update();


    volMapper->SetInputConnection(threshold->GetOutputPort());

    //// Dopasowanie zakresu kolorów do wartości z
    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    double range[2];
    vol->GetScalarRange(range);
    
    volMapper->ScalarVisibilityOn();
    lut->SetTableRange(range[0], range[1]);
    lut->SetNumberOfTableValues(256);

    lut->Build();
    // Ustawienie kolorów w przestrzeni HLS
    volMapper->SetLookupTable(lut);

    vtkNew<vtkActor> volActor;
    volActor->SetMapper(volMapper);
    volActor->GetProperty()->EdgeVisibilityOn();
    volActor->GetProperty()->SetInterpolationToFlat();
    volActor->GetProperty()->SetAmbient(1.0);  // Ustawienie maksymalnej jasności materiału (1.0 = maksymalna jasność)
    volActor->GetProperty()->SetDiffuse(0.2);  // Zmniejszenie odbicia rozproszonego
    volActor->GetProperty()->SetSpecular(0.1);  // Można dodać połysk (dla efektu refleksji)
    volActor->GetProperty()->SetSpecularPower(1.0);  // Zwiększenie refleksów świetlnych
    volActor->GetProperty()->SetColor(colors->GetColor3d("Salmon").GetData());

    renderer->AddActor(volActor);
    renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());
    window->AddRenderer(renderer);

    
    VolumeRenderer volRenderer(renderer, vol);
    vtkNew<vtkRenderWindowInteractor> interactor;
    vtkNew<MouseInteractorStyle> style;

    style->setVolumeRenderer(&volRenderer); // Ustaw instancję VolumeRenderer
    interactor->SetInteractorStyle(style);


	volRenderer.updateVolume();
    QObject::connect(&frequencyXSlider, &QSlider::valueChanged, [&volRenderer](int value) {
        float frequencyX = (value / 10.0f ) + 10.0f;
        volRenderer.setFrequencyX(frequencyX);
        });
    QObject::connect(&frequencyYSlider, &QSlider::valueChanged, [&volRenderer](int value) {
        float frequencyY = (value / 10.0f) + 10.0f;  // (1-100) -> (1-50)
        volRenderer.setFrequencyY(frequencyY);
        });
    QObject::connect(&amplitudeSlider, &QSlider::valueChanged, [&volRenderer](int value) {
        float amplitude = (value / 10.0f) + 10.0f;
        volRenderer.setAmplitude(amplitude);

        });
    QObject::connect(&sizeCubeSlider, &QSlider::valueChanged, [&volRenderer](int value) {
        float sizeCube = (value);
        volRenderer.setSizeCube(sizeCube);

        });
    // Połączenie przycisków z funkcjami
    QObject::connect(&RotateCameraRight, &QPushButton::clicked, [&]() {
        rotateCameraRight(camera, renderer);
        window->Render();
        });

    QObject::connect(&RotateCameraLeft, &QPushButton::clicked, [&]() {
        rotateCameraLeft(camera, renderer);
        window->Render();
        });
    QObject::connect(&playStopButton, &QPushButton::clicked, [&volRenderer]() {
        volRenderer.isPlaying = !volRenderer.isPlaying;
        });

    // Animacja
    QTimer timer;
    QObject::connect(&timer, &QTimer::timeout, [&volRenderer]() {
		if (volRenderer.isPlaying)
            volRenderer.updateTime(0.10);  
        });
    timer.start(16);  // 60 fps (1000 ms / 60 ≈ 16 ms na klatkę)

    mainWindow.show();
    renderer->GetRenderWindow()->SetInteractor(interactor);
    interactor->Initialize();
    interactor->Start();
    

    return app.exec();
}

