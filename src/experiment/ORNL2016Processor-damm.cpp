/** \file ORNL2016Processor.cpp
 * \brief A class to process data from the ORNL 2016 OLTF experiment using
 * VANDLE.
 *
 *\author S. V. Paulauskas
 *\date February 10, 2016
 *
 *\Edits by Thomas King 
 *\Starting April 2016
 */
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <cmath>

#include "BarBuilder.hpp"
#include "DammPlotIds.hpp"
#include "DetectorDriver.hpp"
#include "GetArguments.hpp"
#include "Globals.hpp"
#include "RawEvent.hpp"
#include "TimingMapBuilder.hpp"
#include "ORNL2016Processor.hpp"
#include "DoubleBetaProcessor.hpp"
#include "TreeCorrelator.hpp"
#include "LogicProcessor.hpp"

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

      //Toby Adds
      const int D_BETASCALARVSTIME= 20+ORNL2016_OFFSET;

      //Ge Vs Cycle both Raw and calibrate
      const int DD_GEXVSTIME = 30+ORNL2016_OFFSET; //Full Histogram # is 3300
      const int DD_CALGEXVSTIME = 34+ORNL2016_OFFSET;
      //RAW NaI vs cycle
            const int DD_RAWNAIXVSTIME = 38+ORNL2016_OFFSET; 
      //RAW HAGRiD vs cycle
           const int DD_RAWHAGXVSTIME = 48+ORNL2016_OFFSET;
      //63 is the next free id
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

    //TOBY ADDS
    static int cycleCount = SA; // Sets max ploted cycles for the "per cycle" histograms
   
    DeclareHistogram1D(D_BETASCALARVSTIME,cycleCount,"Beta scalar per cycle");

    //Declaring Ge vs Cycle
    for (unsigned int i=0; i < 4; i++){
      //static  int n=1;
      stringstream ss;
      stringstream sss;
      //      int odd=i+1;
      ss<< "Raw  Energy/2 VS Cycle Number Ge " << i ;
      sss<< "Cal  Energy/2 VS Cycle Number Ge " << i ;
      DeclareHistogram2D(DD_GEXVSTIME + i,SD,cycleCount,ss.str().c_str());
      DeclareHistogram2D(DD_CALGEXVSTIME + i,SD,cycleCount,sss.str().c_str());
      //n=n+1;
    }

    //Declaring NaI vs Cycle
    for (unsigned int i=0; i < 10; i++){
      static int n=1;
      stringstream ss;
      ss<< "Raw Energy/2 VS Cycle Number NaI# " << n ;
      DeclareHistogram2D(DD_RAWNAIXVSTIME + i,SD,cycleCount,ss.str().c_str());
      n=n+1;
    }
    //Declaring HAGRiD vs Cycle
    for (unsigned int i = 0; i < 16; i++){
      static int n=1;
      stringstream ss;
      ss<< "Raw Energy/2 VS Cycle Number HAGRiD# " << n;
      DeclareHistogram2D(DD_RAWHAGXVSTIME + i,SD,cycleCount,ss.str().c_str());
      n=n+1;
    }

 
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

    static bool hasbeta=false;
    for(vector<ChanEvent*>::const_iterator bIt = betaEvts.begin(); 
	bIt != betaEvts.end(); bIt++) {
      plot(D_BETASCALARVSTIME,cycleNum ); //PLOTTING BETA SCALAR RATE (HIS# 759) per CYCLE
      if((*bIt)->GetID() == 158 || (*bIt)->GetID() == 159)
	hasbeta=true;
    }

    for(vector<ChanEvent*>::const_iterator naiIt = naiEvts.begin();
	naiIt != naiEvts.end(); naiIt++) {
      int nainum= (*naiIt)->GetChanID().GetLocation();

      //Filling NaI vs Cycle
      if (TreeCorrelator::get()->place("Cycle")->status()){
	plot(DD_RAWNAIXVSTIME + nainum,( (*naiIt)->GetEnergy())/2,cycleNum);
      }
    } //NaI loop End
    
    for(vector<ChanEvent*>::const_iterator itGe = geEvts.begin();
	itGe != geEvts.end(); itGe++) {
      int genum = (*itGe)->GetChanID().GetLocation();
      if (TreeCorrelator::get()->place("Cycle")->status()){
	plot(DD_GEXVSTIME + genum,(*itGe)->GetEnergy()/2,cycleNum);
	plot(DD_CALGEXVSTIME + genum,(*itGe)->GetCalEnergy()/2,cycleNum);
      }

    } //GE loop end
    
    for(vector<ChanEvent*>::const_iterator itHag = labr3Evts.begin();
	itHag != labr3Evts.end(); itHag++){
      int hagnum = (*itHag)->GetChanID().GetLocation();
      if (TreeCorrelator::get()->place("Cycle")->status()){
	plot(DD_RAWHAGXVSTIME + hagnum, ((*itHag)->GetEnergy()/2),cycleNum);
      }
    } //Hagrid loop end	  
    
      


    EndProcess();
    hasbeta=false;  
    return(true);

}
