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

class ScintAddBack {
    public:
        /** Default constructor setting things to zero */
        ScintAddBack() {
            energy = time = multiplicity = 0;
        }

        /** Default constructor setting default values
         * \param [in] ienergy : the initial energy
         * \param [in] itime : the initial time
         * \param [in] imultiplicity : multiplicity of the event */
        ScintAddBack(double ienergy, double itime, unsigned imultiplicity) {
            energy = ienergy;
            time = itime;
            multiplicity = imultiplicity;
        }

        double energy;//!< Energy of the addback event
        double time;//!< time of the addback event
        unsigned multiplicity;//!< multiplicity of the event
};

/// Class to process VANDLE analysis for ORNL2016 campaign  at the OLTF
class ORNL2016Processor : public EventProcessor { 
public:
    /** Default Constructor */
    ORNL2016Processor(double gamma_threshold_L, double sub_event_L, double gamma_threshold_N, double sub_event_N);
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

      /** Returns the events that were added to the LaBr3 addback_ */
  std::vector< std::vector<ScintAddBack> > GetLaddBackEvent(void) {return(LaddBack_);}
  std::vector< std::vector<ScintAddBack> > GetNaddBackEvent(void) {return(NaddBack_);}
  

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
    int bMulti;
  } aux;


  TFile *rootFName_;
    //functions for root preocessing
  void rootArrayreset(double arrayName[], int arraySize);

  void rootGstrutInit(RAY &strutName);
  
  //thresholds and windows for gamma addback for LaBr3 (L*) and NaI (N*)
  double LgammaThreshold_;
  double LsubEventWindow_;
  double NgammaThreshold_;
  double NsubEventWindow_;

  /** addbackEvents vector of vectors, where first vector
   * enumerates crystal, second events */
  std::vector< std::vector<ScintAddBack> > LaddBack_;
  std::vector< std::vector<ScintAddBack> > NaddBack_;



};
#endif


