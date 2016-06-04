set(EXPERIMENT_SOURCES 
  ORNL2016Processor.cpp

)

add_library(ExperimentObjects OBJECT ${EXPERIMENT_SOURCES})
add_library(ExperimentObjectsTOF OBJECT ORNL2016Processor-TOF.cpp)
