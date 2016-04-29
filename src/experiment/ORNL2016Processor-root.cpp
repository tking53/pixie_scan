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
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TTree.h"

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

      // const int DD_CYCVSGE0 = 30+ORNL2016_OFFSET;
      // const int DD_CYCVSGE0R = 31+ORNL2016_OFFSET;
      // const int D_GE0 = 34+ORNL2016_OFFSET;
      // const int D_GE0R = 35+ORNL2016_OFFSET;
      
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
    // DeclareHistogram2D(DD_CYCVSGE0 , SD ,S7, "Apples+apples GE[0]");
    // DeclareHistogram2D(DD_CYCVSGE0R , SD ,S7, "Apples+apples RAW GE[0]");
    // DeclareHistogram1D(D_GE0 , SD , "Apples+apples GE[0] remove of 1313");
    // DeclareHistogram1D(D_GE0R , SD , "Apples+apples GE[0] remove of 113");

}

void ORNL2016Processor::rootArrayreset(double arrayName[], int arraySize){ //refills arrayName with 0s 
  fill(arrayName,arrayName + arraySize,0);
}

void ORNL2016Processor::rootGstrutInit(RAY &strutName) { //Zeros the entire root gamma structure

  fill(strutName.Hag,strutName.Hag + 16,0);
  fill(strutName.NaI,strutName.NaI + 10,0);
  fill(strutName.Ge,strutName.Ge + 4,0);
  strutName.beta = 0;
  strutName.cycle=0;

}

ORNL2016Processor::ORNL2016Processor(const std::vector<std::string> &typeList,
    const double &res, const double &offset, const double &numStarts) :
    VandleProcessor(typeList,res,offset,numStarts) {
    associatedTypes.insert("vandle");
    associatedTypes.insert("liglass");
    associatedTypes.insert("nai");
    associatedTypes.insert("labr3");
    associatedTypes.insert("beta");
    char hisFileName[32];
    GetArgument(1,hisFileName,32);
    string tmp = hisFileName;
    tmp= tmp.substr(0,tmp.find_first_of(" "));
    stringstream rootname;
    rootname<<tmp<<".root";
    rootFName_ =  new TFile(rootname.str().c_str(),"RECREATE");
    tree = new TTree("gammas","Tree containing gamma events with cycle and betas");
    //    calbranch = tree->Branch("calgam",&calgam,"Hag[16]/D:NaI[10]/D:Ge[4]/D:beta/D:cycle/i");
    rawbranch = tree->Branch("rawgam",&rawgam,"Hag[16]/D:NaI[10]/D:Ge[4]/D:beta/D:cycle/i");
    tree->SetAutoFlush(3000);
    //    rootGstrutInit(calgam);
    rootGstrutInit(rawgam);
}


ORNL2016Processor::~ORNL2016Processor(){
  rootFName_->Write();
  rootFName_->Close();
  
  delete(rootFName_);


}
 

bool ORNL2016Processor::PreProcess(RawEvent &event) {
    if (!EventProcessor::PreProcess(event))
        return(false);
   

    EndProcess();
    return(true);
}

bool ORNL2016Processor::Process(RawEvent &event) {
    if (!EventProcessor::Process(event))
        return(false);
  

    static const vector<ChanEvent*> & betaEvts =
        event.GetSummary("beta:double")->GetList();
    static const vector<ChanEvent*> &labr3Evts =
        event.GetSummary("labr3")->GetList();
    static const vector<ChanEvent*> &naiEvts =
        event.GetSummary("nai")->GetList();
    static const vector<ChanEvent*> &geEvts =
        event.GetSummary("ge")->GetList();


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
    //resets root struct
    //rootGstrutInit(calgam);
    rootGstrutInit(rawgam);
      
   //Cycle timing
    static double cycleLast = 2;
    static int cycleNum = 0;
  if (TreeCorrelator::get()->place("Cycle")->status()){	  
      double cycleTime = TreeCorrelator::get()->place("Cycle")->last().time;
      cycleTime *= (Globals::get()->clockInSeconds()*1.e9);
      if ( cycleTime != cycleLast ){

	double tdiff = (cycleTime - cycleLast) / 1e6; //Outputs cycle length in msecs.
	if (cycleNum == 0){cout<<" #  There are some events at the beginning of the first segment missing from Histograms that use cycleNum."<<endl<<" #  This is a product of not starting the cycle After the LDF."<<endl<<" #  This First TDIFF is most likly nonsense"<<endl;}
	cycleLast = cycleTime;
	cycleNum++;
	cout<< "Cycle Change "<<endl<<"Tdiff (Cycle start and Now) (ms)= "<<tdiff<<endl<<"Starting on Cycle #"<<cycleNum<<endl; 
      }
  	
  }
  //  calgam.cycle = cycleNum; //ROOT Set
  rawgam.cycle = cycleNum;
  
  //Betas
    static bool hasbeta=false;
    for(vector<ChanEvent*>::const_iterator bIt = betaEvts.begin(); 
	bIt != betaEvts.end(); bIt++) {
      

      //calgam.beta = (*bIt)->GetCalEnergy();
      rawgam.beta = (*bIt)->GetEnergy();
      
      // if((*bIt)->GetID() == 158 || (*bIt)->GetID() == 159)
      // 	hasbeta=true;
      
    }
    
    for(vector<ChanEvent*>::const_iterator naiIt = naiEvts.begin();
	naiIt != naiEvts.end(); naiIt++) {
      int nainum= (*naiIt)->GetChanID().GetLocation();
     
      //calgam.NaI[nainum] = (*naiIt)->GetCalEnergy();
      rawgam.NaI[nainum] = (*naiIt)->GetEnergy();
   

      //if (TreeCorrelator::get()->place("Cycle")->status()){

      //}
    } //NaI loop End

    for(vector<ChanEvent*>::const_iterator itGe = geEvts.begin();
	itGe != geEvts.end(); itGe++) {
      int genum = (*itGe)->GetChanID().GetLocation();
    
      // if (genum == 0){
      // 	plot(DD_CYCVSGE0,(*itGe)->GetCalEnergy(),cycleNum);
      // 	plot(DD_CYCVSGE0R,(*itGe)->GetEnergy(),cycleNum);
      // 	plot(D_GE0,(*itGe)->GetCalEnergy());
      // 	plot(D_GE0R,(*itGe)->GetEnergy());
      // }
      // if (TreeCorrelator::get()->place("Cycle")->status()){
      // }
      rawgam.Ge[genum] = (*itGe)->GetEnergy();
      //      calgam.Ge[genum] = (*itGe)->GetCalEnergy();
      
     
    } //GE loop end

    for(vector<ChanEvent*>::const_iterator itHag = labr3Evts.begin();
	itHag != labr3Evts.end(); itHag++){
      int hagnum = (*itHag)->GetChanID().GetLocation();
   

      rawgam.Hag[hagnum]= (*itHag)->GetEnergy();
      //     calgam.Hag[hagnum]= (*itHag)->GetCalEnergy();
   
      // if (TreeCorrelator::get()->place("Cycle")->status()){
      // }

    } //Hagrid loop end	  
    
    tree->Fill();      


    EndProcess();
    hasbeta=false;  
    return(true);

}
