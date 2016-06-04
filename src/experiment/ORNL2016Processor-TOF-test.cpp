/** \file ORNL2016Processor-TOF.cpp
 * \brief A class to process data from the ORNL 2016 OLTF experiment using
 * VANDLE. Using Root and Damm for histogram analysis. 
 * THIS CPP has TOF Filling
 *
 *\author S. V. Paulauskas
 *\date February 10, 2016
 *
 *\Edits by Thomas King 
 *\Starting April 2016
 * 
 *
*/
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include "BarBuilder.hpp"
#include "RawEvent.hpp"
#include "TimingMapBuilder.hpp"
#include "TreeCorrelator.hpp"
#include "LogicProcessor.hpp"
#include "VandleProcessor.hpp"

#include "DammPlotIds.hpp"
#include "DetectorDriver.hpp"
#include "GetArguments.hpp"
#include "Globals.hpp"


#include "ORNL2016Processor.hpp"
#include "DoubleBetaProcessor.hpp"


namespace dammIds {
    namespace vandle {
        const unsigned int ORNL2016_OFFSET = 70;
	const int D_VANDLEMULT  = 0+ORNL2016_OFFSET;

        const int DD_CTOFNOTAPE  = 1+ORNL2016_OFFSET;
        const int DD_QDCTOFNOGATE  = 2+ORNL2016_OFFSET;
        const int DD_CORTOFVSEGAM  = 3+ORNL2016_OFFSET;
        const int DD_QDCVSCORTOFMULT1  = 4+ORNL2016_OFFSET;
        const int DD_MULT2SYM  = 5+ORNL2016_OFFSET;
	const int DD_LIGLEN = 6+ORNL2016_OFFSET;
	const int DD_LIGLTOF = 7+ORNL2016_OFFSET;

	const int D_LABR3SUM = 8+ORNL2016_OFFSET;
	const int D_LABR3BETA = 9+ORNL2016_OFFSET;

	const int D_NAISUM = 10+ORNL2016_OFFSET;
	const int D_NAIBETA = 11+ORNL2016_OFFSET;

	const int DD_TOFVSNAI = 12+ORNL2016_OFFSET;
	const int DD_TOFVSHAGRID = 13+ORNL2016_OFFSET;
	const int DD_TOFVSGE = 14+ORNL2016_OFFSET;

    
    }
}//namespace dammIds

using namespace std;
using namespace dammIds::vandle;

void ORNL2016Processor::DeclarePlots(void) {
    VandleProcessor::DeclarePlots();
    DeclareHistogram1D(D_VANDLEMULT, S7, "Vandle Multiplicity");
    DeclareHistogram2D(DD_QDCTOFNOGATE, SC, SD, "QDC ToF Ungated");
    DeclareHistogram2D(DD_QDCVSCORTOFMULT1, SC, SC, "QDC vs Cor Tof Mult1");
    DeclareHistogram2D(DD_LIGLEN, SC, S5, "E - LiGlass");
    DeclareHistogram2D(DD_LIGLTOF, SC, SC, "E vs ToF - LiGlass");
    
    DeclareHistogram1D(D_LABR3SUM, SC, "HAGRiD summed");
    DeclareHistogram1D(D_LABR3BETA, SC, "HAGRiD summed - BETA GATED");
    DeclareHistogram1D(D_NAISUM, SC, "NaI summed");
    DeclareHistogram1D(D_NAIBETA, SC, "NaI summed - BETA GATED");

    DeclareHistogram2D(DD_TOFVSNAI, SC, SB, "ToF vs. NaI");
    DeclareHistogram2D(DD_TOFVSHAGRID, SC, SB, "ToF vs. HAGRiD");
    DeclareHistogram2D(DD_TOFVSGE, SC, SB, "ToF vs. Ge");

}


ORNL2016Processor::ORNL2016Processor(const std::vector<std::string> &typeList,
    const double &res, const double &offset, const double &numStarts) :
    VandleProcessor(typeList,res,offset,numStarts) {
    associatedTypes.insert("vandle");
    associatedTypes.insert("liglass");
    associatedTypes.insert("nai");
    associatedTypes.insert("labr3");
    associatedTypes.insert("beta");
}

bool ORNL2016Processor::PreProcess(RawEvent &event) {
    if (!EventProcessor::PreProcess(event))
        return(false);
    if (!VandleProcessor::PreProcess(event))
	return(false);

    EndProcess();
    return(true);
}

bool ORNL2016Processor::Process(RawEvent &event) {
    if (!EventProcessor::Process(event))
        return(false);
    if(!VandleProcessor::Process(event))
        return(false);

    static const vector<ChanEvent*> & betaEvts =
        event.GetSummary("beta:double")->GetList();
    static const vector<ChanEvent*> &liEvts =
        event.GetSummary("liglass")->GetList();
    static const vector<ChanEvent*> &labr3Evts =
        event.GetSummary("labr3")->GetList();
    static const vector<ChanEvent*> &naiEvts =
        event.GetSummary("nai")->GetList();
    static const vector<ChanEvent*> &geEvts =
        event.GetSummary("ge")->GetList();

    BarBuilder builder(betaEvts);
    BarMap betas = builder.GetBarMap();
    
    bool hasMultOne = bars_.size() == 1;

    if(bars_.size() != 0)
	plot(D_VANDLEMULT, bars_.size());
    
    for(BarMap::iterator itStart = betas.begin();
	itStart != betas.end(); itStart++) {
	
	unsigned int startLoc = (*itStart).first.first;
	BarDetector start = (*itStart).second;
	
	if(bars_.size() != 0) {
	    for (BarMap::iterator it = bars_.begin(); it !=  bars_.end(); it++) {
		TimingDefs::TimingIdentifier barId = (*it).first;
		BarDetector bar = (*it).second;
		
		if(!bar.GetHasEvent())
		    continue;
		
		TimingCalibration cal = bar.GetCalibration();
		
		double tofOffset = cal.GetTofOffset(startLoc);
		double tof = bar.GetWalkCorTimeAve() -
		    start.GetWalkCorTimeAve() + tofOffset;
		
		double corTof =
		    CorrectTOF(tof, bar.GetFlightPath(), cal.GetZ0());
		
		plot(DD_QDCTOFNOGATE, tof*plotMult_+plotOffset_, bar.GetQdc());
		if(hasMultOne)
		    plot(DD_QDCVSCORTOFMULT1, corTof*plotMult_+plotOffset_, bar.GetQdc());
		
		for(vector<ChanEvent*>::const_iterator it = labr3Evts.begin(); 
		    it != labr3Evts.end(); it++)
		    plot(DD_TOFVSHAGRID, (*it)->GetCalEnergy(), tof*plotMult_+200);
		
		for(vector<ChanEvent*>::const_iterator naiIt = naiEvts.begin(); 
		    naiIt != naiEvts.end(); naiIt++)
		    plot(DD_TOFVSNAI, (*naiIt)->GetCalEnergy(), tof*plotMult_+200);
		
		for (vector<ChanEvent *>::const_iterator itGe = geEvts.begin();
		     itGe != geEvts.end(); itGe++) 
		    plot(DD_TOFVSGE, (*itGe)->GetCalEnergy(), tof*plotMult_+200);
		//		for (vector<ChanEvent *>::const_iterator itGe = geEvts.begin();
		//   itGe != geEvts.end(); itGe++) 
		//  plot(DD_GAMMVSTIME, (*itGe)->GetCalEnergy(), tof*plotMult_+200);
	    } //(BarMap::iterator itBar
	} //if(bars_.size != 0)
	
	for(vector<ChanEvent*>::const_iterator liIt = liEvts.begin(); 
	    liIt != liEvts.end(); liIt++) {
	    unsigned int id = (*liIt)->GetID();
	    
	    //Beta gated energy
	    plot(DD_LIGLEN, (*liIt)->GetEnergy(), id);

	    double li_tof = ((*liIt)->GetTime() - 0.5*(start.GetLeftSide().GetPixieTime()+
						       start.GetRightSide().GetPixieTime()));
	    //cout << id << " " << li_tof+1000 << endl;
	    if((*liIt)->GetID() != 12)
		plot(DD_LIGLTOF, li_tof+1000, (*liIt)->GetEnergy());
	}
    } // for(TimingMap::iterator itStart
    
    ///PLOTTING THE SINGLES SPECTRA FOR THE LIGLASS
    for(vector<ChanEvent*>::const_iterator liIt = liEvts.begin(); 
	liIt != liEvts.end(); liIt++)
	plot(DD_LIGLEN, (*liIt)->GetEnergy(), (*liIt)->GetID()-11);
    
    ///PLOTTING THE SINGLES AND BETA GATED SPECTRA FOR HAGRiD
    for(vector<ChanEvent*>::const_iterator it = labr3Evts.begin(); 
	it != labr3Evts.end(); it++) {
	plot(D_LABR3SUM, (*it)->GetCalEnergy());
	for(vector<ChanEvent*>::const_iterator bIt = betaEvts.begin(); 
	    bIt != betaEvts.end(); bIt++) {
	    if((*bIt)->GetID() == 158 || (*bIt)->GetID() == 159)
		plot(D_LABR3BETA, (*it)->GetCalEnergy());

	}
    }
    
    ///PLOTTING THE SINGLES AND BETA GATED SPECTRA FOR NAI
    for(vector<ChanEvent*>::const_iterator naiIt = naiEvts.begin(); 
     	naiIt != naiEvts.end(); naiIt++) {
	plot(D_NAISUM, (*naiIt)->GetCalEnergy());
	for(vector<ChanEvent*>::const_iterator bIt = betaEvts.begin(); 
	    bIt != betaEvts.end(); bIt++) {
	    if((*bIt)->GetID() == 158 || (*bIt)->GetID() == 159)
		plot(D_NAIBETA, (*naiIt)->GetCalEnergy());
	    
	}
    }

   
   /// PLOT ANALYSIS HISTOGRAMS

    //Cycle timing
    static double cycleLast = 2;
    static int cycleNum = 0;
  if (TreeCorrelator::get()->place("Cycle")->status()){	  
      double cycleTime = TreeCorrelator::get()->place("Cycle")->last().time;
      cycleTime *= (Globals::get()->clockInSeconds()*1.e9);
      if ( cycleTime != cycleLast ){
	double tdiff = (cycleTime - cycleLast) / 1e7; //Outputs cycle length in msecs.
	if (cycleNum == 0){	
	  cout<<" #  There are some events at the beginning of the first segment missing from Histograms that use cycleNum."<<endl<<" #  This is a product of not starting the cycle After the LDF."<<endl<<" #  This First TDIFF is most likly nonsense"<<endl;
	}
	cycleLast = cycleTime;
	cycleNum = cycleNum + 1;
	cout<< "Cycle Change "<<endl<<"Tdiff (Cycle start and Now) (ms)= "<<tdiff<<endl<<"Starting on Cycle #"<<cycleNum<<endl; 
      }
  }  


    EndProcess();
    
    return(true);

}
