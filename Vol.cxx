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
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkSphereSource.h>
#include <QGridLayout>
#include <vtkCamera.h>
#include <qtimer.h>
#include <vtkCellPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkWarpScalar.h>
#include <vtkGeometryFilter.h>
#include <vtkFloatArray.h>
#include <vtkGenericDataArray.h>
#include <vtkShortArray.h>
#include <vtkOpenGLGPUVolumeRayCastMapper.h>
#include <vtkVolumeProperty.h>
#include <vtkPiecewiseFunction.h>
#include <vtkOutlineFilter.h>
#include <QThreadPool>
#include <cmath>
#include <chrono>
#include <vtkCharArray.h>

#define M_PI 3.14159265358979323846
vtkNew<vtkActor> sphereActor;
//static int sizeCube = 64; // LOW QUALITY (FAST RENDER)
static int sizeCube = 128; // MEDIUM QUALITY (MEDIUM RENDER)
//static int sizeCube = 256; // HIGH QUALITY (LOW RENDER)
constexpr std::size_t LOOKUP_SIZE = 4096;
static float sp = 1.0 / float(sizeCube - 1);

float sinLookup[LOOKUP_SIZE];
float cosLookup[LOOKUP_SIZE];

float* zLookup = new float[sizeCube];
float* yLookup = new float[sizeCube];
float* xLookup = new float[sizeCube];

float pickedPointXY[2];


void setupLookup() {
    pickedPointXY[0] = 0;
    pickedPointXY[1] = 0;


    for (int i = 0; i < LOOKUP_SIZE; i++) {
        sinLookup[i] = sinf(i / float(LOOKUP_SIZE) * 2 * M_PI);
        cosLookup[i] = cosf(i / float(LOOKUP_SIZE) * 2 * M_PI);
    }
    for (int i = 0; i < sizeCube; i++) {
        zLookup[i] = i * sp;
        yLookup[i] = -0.5 + i * sp;
        xLookup[i] = -0.5 + i * sp;
    }
}
float lightX(int i) {
	return xLookup[i];
}
float lightY(int i) {
    return yLookup[i];
}
float lightZ(int i) {
    return zLookup[i];
}
float lightSin(float angle) {
    return sinLookup[int(angle / 2 * M_PI * LOOKUP_SIZE/4) % LOOKUP_SIZE];
}

float lightCos(float angle) {
    return cosLookup[int(angle / 2 * M_PI * LOOKUP_SIZE/4) % LOOKUP_SIZE];
}


#define BENCHMARK

void benchmark(std::function<void()> f) {
    auto start = std::chrono::high_resolution_clock::now();
    f();
    auto stop = std::chrono::high_resolution_clock::now();

#ifdef BENCHMARK
	auto currentTime = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    if (duration >= 1) {
        std::cout << "FPS: " << int( 1000/duration) << std::endl;
        start = stop;
    }
#endif
}


class VolumeRenderer {
public:
    VolumeRenderer(vtkRenderer* renderer, vtkStructuredPoints* volPoints )
        : renderer(renderer), volPoints(volPoints), frequencyX(5.0f), frequencyY(5.0f), amplitude(5.0f) {}

    bool isPlaying = false;

    float damping(float distance, float time, float decayRate) {
        return exp(-decayRate * distance) * exp(-decayRate *time);
    }
    float gaussian(float x, float y, float sigma) {
        return exp(-(x * x + y * y) / (2 * sigma * sigma));
    }

    void updateVolume() {
        QThreadPool pool;

        float* begin = (float*)volPoints->GetPointData()->GetScalars()->GetVoidPointer(0);
        auto scalars = volPoints->GetPointData()->GetScalars();
     
        std::vector<float> sinValues(sizeCube);
        std::vector<float> cosValues(sizeCube);

        float decayRate = 2.0f; // Szybkość tłumienia zakłócenia
        auto currentTime = std::chrono::high_resolution_clock::now();
        float elapsedTime = std::chrono::duration<float>(currentTime - lastClickTime).count(); // Czas od ostatniego kliknięcia

        for (int i = 0; i < sizeCube; ++i) {
            sinValues[i] = lightSin(i / ((sizeCube - 1) / frequencyX) + time);
            cosValues[i] = lightCos(i / ((sizeCube - 1) / frequencyY) + time);
            
            
        }

            for (int k = 0; k < sizeCube; k++) {
                float z = lightZ(k);
                int kOffset = k * sizeCube * sizeCube;
                pool.start([this, z, kOffset, begin, &sinValues, &cosValues, decayRate, elapsedTime ]() {
                    for (int j = 0; j < sizeCube; j++) {
                        float y = lightY(j);
                        int jOffset = j * sizeCube;
                        for (int i = 0; i < sizeCube; i++) {
                            float x = lightX(i);
                            float value = 0;
                            if (pickedPointXY[0] == 0 || pickedPointXY[1] == 0) {
                                value = z - ((sinValues[i] * cosValues[j]) / ((sizeCube - 1) / amplitude));
                            } else {
                                float dx = pickedPointXY[0] - x;
                                float dy = pickedPointXY[1] - y;
                                float distance = sqrt(dx * dx + dy * dy);
                                float disturbance = 0.0f;
                                for (int n = 1; n <= 3; ++n) { // kolejne fale
                                    float factor = 0.10f / n;
                                    float gaussianValue = gaussian(dx, dy, 0.05f * n);
                                    float sinValue = sin(elapsedTime * 10.0f - distance * 10.0f * n);
                                    float dampingValue = damping(distance, elapsedTime, decayRate);
                                    disturbance += factor * gaussianValue * sinValue * dampingValue;
                                }
                                value = z - ((sinValues[i] * cosValues[j]) / ((sizeCube - 1) / amplitude))  + disturbance;
                            }
                            int offset = i + jOffset + kOffset;
                            begin[offset] = value;
                        }
                    }
                });
        }
        pool.waitForDone();

        volPoints->Modified();
        //volPoints->GetPointData()->SetScalars(scalars);
        renderer->GetRenderWindow()->Render();
        
    }

    void setFrequencyX(float freqx) {
        frequencyX = freqx;
        updateVolume();
    }
    void setFrequencyY(float freqy) {
        frequencyY = freqy;
        updateVolume();
    }
    void setAmplitude(float amp) {
        amplitude = amp;
        updateVolume();
    }
	void setSizeCube(int size) {
        volPoints->SetDimensions(size, size, size);
        auto sp = 1.0 / float(size - 1);
        volPoints->SetSpacing(sp, sp, sp);
		updateVolume();
	}
    void updateTime(float deltaTime) {
        time += deltaTime; // Zwiększanie czasu
        updateVolume();
    }
    void pickPoint(int x, int y) {
        vtkNew<vtkCellPicker> picker;
        picker->SetTolerance(0.005);

        if (picker->Pick(x, y, 0, renderer)) {
            double pickedPosition[3];
            picker->GetPickPosition(pickedPosition);

            vtkIdType pointId = volPoints->FindPoint(pickedPosition);
            if (pointId != -1) {
                std::cout << "Picked Point ID: " << pointId << std::endl;
                std::cout << "Picked Position: (" << pickedPosition[0] << ", "
                    << pickedPosition[1] << ", " << pickedPosition[2] << ")" << std::endl;

                pickedPointXY[0] = pickedPosition[0];
                pickedPointXY[1] = pickedPosition[1];
                std::cout << "Amp: " << amplitude << " freqX: " << frequencyX << " freqY: " << frequencyY << std::endl;
                //changePointScalar(pointId, 100.0); // Nowa wartość skalarna
                //increasePointHeight(pointId, -10.0);
                //updateVoxels();
                //applyWarpScalar(); // Zastosuj deformację
                // Pobierz wartość skalaru dla wybranego punktu
                lastClickTime = std::chrono::high_resolution_clock::now(); // Resetowanie czasu ostatniego kliknięcia

                vtkDataArray* scalars = volPoints->GetPointData()->GetScalars();
                if (scalars) {
                    float value = scalars->GetTuple1(pointId);
					pickedPointScalar = value;
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
        vtkDataArray* scalars = volPoints->GetPointData()->GetScalars();
        if (scalars && pointId != -1) {
            scalars->SetTuple1(pointId, newValue); // Zmiana wartości skalarnej
            volPoints->Modified();
            renderer->GetRenderWindow()->Render(); 
        }
    }
    

private:
    vtkRenderer* renderer;
    vtkStructuredPoints* volPoints;
    float frequencyX;
    float frequencyY;
    float amplitude;
    float pickedPointScalar;
    std::chrono::time_point<std::chrono::high_resolution_clock> lastClickTime; // Czas ostatniego kliknięcia
    float time; // Czas, który będzie animował falę
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

        //vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
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
    //vtkNew<vtkDataSetMapper> volMapper;
    vtkNew<vtkOpenGLGPUVolumeRayCastMapper> volMapper;
    vtkNew<vtkVolumeProperty> volProperty;
    vtkNew<vtkRenderer> renderer;
    setupLookup();

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
    amplitudeSlider.setRange(-1000, 1000);
    amplitudeSlider.setValue(0);

    QSlider sizeCubeSlider(Qt::Horizontal);
    sizeCubeSlider.setRange(3, 256);
    sizeCubeSlider.setValue(256);




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
    //sliderLayout->addWidget(&label_SizeCube);
    //sliderLayout->addWidget(&sizeCubeSlider);
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

    vtkNew<vtkStructuredPoints> volPoints;

    volPoints->SetDimensions(sizeCube, sizeCube, sizeCube);
    volPoints->SetOrigin(-0.5, -0.5, -0.5);

    volPoints->SetSpacing(sp, sp, sp);

    vtkNew<vtkFloatArray> scalars;
    scalars->SetNumberOfComponents(1);
    scalars->SetNumberOfTuples(sizeCube * sizeCube * sizeCube);
    
    volPoints->GetPointData()->SetScalars(scalars);

    volMapper->SetInputData(volPoints);
    
    vtkNew<vtkPiecewiseFunction> opacityTransferFunction;
    opacityTransferFunction->AddPoint(0.1, 0.9999);
    opacityTransferFunction->AddPoint(0.3, 0.999);
    opacityTransferFunction->AddPoint(0.45, 0.9999);
    opacityTransferFunction->AddPoint(0.50, 1.0);
    opacityTransferFunction->AddPoint(0.51, 0.0);
    //opacityTransferFunction->AddPoint(0.5, 1.0);
    //opacityTransferFunction->AddPoint(0.51, 0.0);
	vtkSmartPointer<vtkColorTransferFunction> colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
	colorTransferFunction->AddRGBPoint(0.1, 0.45, 0.10, 0.0);
	colorTransferFunction->AddRGBPoint(0.8, 0.2, 0.35, 0.0);
	colorTransferFunction->AddRGBPoint(0.4, 0.7, 0.6, 1.0);

    volProperty->SetColor(colorTransferFunction);
    volProperty->SetScalarOpacity(opacityTransferFunction);
    volProperty->SetInterpolationTypeToLinear();
    volProperty->ShadeOn();
    volProperty->SetAmbient(0.4);
    volProperty->SetDiffuse(0.6);
    volProperty->SetSpecular(0.2);
    volProperty->SetSpecularPower(10.0);


    volMapper->SetAutoAdjustSampleDistances(true);
	vtkNew<vtkVolume> vol;
    vol->SetMapper(volMapper);
    vol->SetProperty(volProperty);

    renderer->AddVolume(vol);

    // Add wireframe outline
    vtkNew<vtkOutlineFilter> outlineFilter;
    outlineFilter->SetInputData(volPoints);
    outlineFilter->Update();


    vtkNew<vtkPolyDataMapper> outlineMapper;
    outlineMapper->SetInputConnection(outlineFilter->GetOutputPort());


    

    vtkNew<vtkActor> outlineActor;
    outlineActor->SetMapper(outlineMapper);
    outlineActor->GetProperty()->SetColor(1.0, 1.0, 1.0); // Ustaw kolor siatki na biały
    outlineActor->GetProperty()->SetRepresentationToWireframe();

    renderer->AddActor(outlineActor);


    //sphereActor->VisibilityOff();
    

    //renderer->AddActor(volActor);

    renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());
    window->AddRenderer(renderer);

    
    VolumeRenderer volRenderer(renderer, volPoints);
    vtkNew<vtkRenderWindowInteractor> interactor;
    vtkNew<MouseInteractorStyle> style;

    style->setVolumeRenderer(&volRenderer); // Ustaw instancję VolumeRenderer
    interactor->SetInteractorStyle(style);


	volRenderer.updateVolume();
    QObject::connect(&frequencyXSlider, &QSlider::valueChanged, [&volRenderer](int value) {
        float frequencyX = (value / 5.0f ) + 5.0f;
        volRenderer.setFrequencyX(frequencyX);
        });
    QObject::connect(&frequencyYSlider, &QSlider::valueChanged, [&volRenderer](int value) {
        float frequencyY = (value / 5.0f) + 5.0f;  // (1-100) -> (1-50)
        volRenderer.setFrequencyY(frequencyY);
        });
    QObject::connect(&amplitudeSlider, &QSlider::valueChanged, [&volRenderer](int value) {
        float amplitude = (value / 5.0f) + 5.0f;
        volRenderer.setAmplitude(amplitude);

        });
    /*QObject::connect(&sizeCubeSlider, &QSlider::valueChanged, [&volRenderer](int value) {
        float sizeCube = (value);
        volRenderer.setSizeCube(sizeCube);

        });*/
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
        benchmark([&]() {
            if (volRenderer.isPlaying)
                volRenderer.updateTime(0.01);
            });
        
     });
		
    timer.start(16);  // 60 fps (1000 ms / 60 ≈ 16 ms na klatkę)

    mainWindow.show();
    renderer->GetRenderWindow()->SetInteractor(interactor);
    interactor->Initialize();
    interactor->Start();
    

    return app.exec();
}

