/** \file ORNL2016Processor.hpp
 * \brief A class to process data from the OLTF experiments in Feb/March 2016 using
 * VANDLE at ORNl
 *
 *\author S. V. Paulauskas
 *\date Feb, 2016
 */

#ifndef __ORNL2016PROCESSOR_HPP_
#define __ORNL2016PROCESSOR_HPP_
#include <fstream>


#include "EventProcessor.hpp"


#ifdef useroot
#include <TFile.h>
#include <TTree.h>
#endif

/// Class to process VANDLE analysis for ORNL2016 campaign  at the OLTF
class ORNL2016Processor : public EventProcessor { 
public:
    /** Default Constructor */
    ORNL2016Processor();
    /** Default Destructor */
    ~ORNL2016Processor() ;
    /** Declare the plots used in the analysis */
    virtual void DeclarePlots(void);


    /** Preprocess the ORNL2016 data
     * \param [in] event : the event to preprocess
     * \return true if successful */
    virtual bool PreProcess(RawEvent &event);

    /** Process the event
    * \param [in] event : the event to process
    * \return Returns true if the processing was successful */
    virtual bool Process(RawEvent &event);



private:  
  TTree *tree ;
  TBranch *calbranch;
  TBranch *rawbranch;
  TBranch *auxBranch;

  struct RAY { 
    double Hag[16];
    double NaI[10];
    double Ge[4];
    double beta;
    int cycle;
    int eventNum;
    int gMulti;
    int nMulti;
    int hMulti;
  } aux;


  TFile *rootFName_;
    //functions for root preocessing
  void rootArrayreset(double arrayName[], int arraySize);

  void rootGstrutInit(RAY &strutName);
};
#endif

