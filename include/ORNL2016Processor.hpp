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

#include "DoubleBetaProcessor.hpp"
#include "EventProcessor.hpp"
#include "VandleProcessor.hpp"

#ifdef useroot
#include <TFile.h>
#include <TTree.h>
#endif

/// Class to process VANDLE analysis for ORNL2016 campaign  at the OLTF
class ORNL2016Processor : public VandleProcessor {
public:
    /** Default Constructor */
    ORNL2016Processor();
    /** Default Destructor */
    ~ORNL2016Processor() ;
    /** Declare the plots used in the analysis */
    virtual void DeclarePlots(void);

    /** Constructor taking a list of detector types as an argument
     * \param [in] typeList : the list of bar types that are in the analysis
     * \param [in] res : The resolution of the DAMM histograms
     * \param [in] offset : The offset of the DAMM histograms */
    ORNL2016Processor(const std::vector<std::string> &typeList,
                    const double &res, const double &offset,
                    const double &numStarts);

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

  struct RAY { 
    double Hag[16];
    double NaI[10];
    double Ge[4];
    double beta;
    int cycle;
  } calgam, rawgam;

  TFile *rootFName_;
  
};
#endif
