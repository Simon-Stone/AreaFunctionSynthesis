cmake_minimum_required(VERSION 3.24)
project(AreaFunctionSynthesis)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(portaudio-19.7.0 EXCLUDE_FROM_ALL)
add_subdirectory(wxWidgets-3.2.1)

add_executable(AreaFunctionSynthesis src/Frontend/Application.cpp src/Frontend/MainWindow.cpp src/Frontend/MainPage.cpp src/Frontend/AreaFunctionPicture.cpp src/Frontend/Graph.cpp src/Frontend/BasicPicture.cpp src/Frontend/Data.h src/Backend/Synthesizer.cpp src/Backend/TdsModel.cpp src/Backend/TlModel.cpp src/Backend/LfPulse.cpp src/Backend/Tube.cpp src/Backend/TriangularGlottis.cpp src/Backend/Dsp.cpp src/Backend/IirFilter.cpp src/Backend/Signal.cpp src/Backend/TimeFunction.cpp src/Backend/OneDimAreaFunction.cpp src/Frontend/Data.cpp src/Backend/Glottis.cpp src/Backend/Matrix2x2.cpp src/Backend/XmlNode.cpp src/Frontend/SynthesisThread.cpp)
target_link_libraries(AreaFunctionSynthesis wx::net wx::core wx::base portaudio)
