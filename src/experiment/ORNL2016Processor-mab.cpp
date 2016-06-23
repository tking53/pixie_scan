/** \file ORNL2016Processor-mab.cpp
 * \brief A class to process data from the ORNL 2016 OLTF experiment using
 * VANDLE. Using Root and Damm for histogram analysis. 
 * Generates aux root branch  with add backs and multiplicies 
 * And pab root branch with procesor addback.
 *
 *\author S. V. Paulauskas
 *\date February 10, 2016
 *
 *\Edits by Thomas King 
 *\Starting June 2016
 * 
 *
*/
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>


#include "DammPlotIds.hpp"
#include "DetectorDriver.hpp"
#include "GetArguments.hpp"
#include "Globals.hpp"

#include "ORNL2016Processor.hpp"
#include "DoubleBetaProcessor.hpp"
#include "GeProcessor.hpp"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"

 static unsigned int evtNum = 0;

namespace dammIds {
    namespace experiment {
      const unsigned int ORNL2016_OFFSET = 0;
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
      
      //seeds for Damm cycle his
      const int D_BETASCALARRATE= 29+ORNL2016_OFFSET; //6079 in his
      const int D_BETAENERGY=30+ORNL2016_OFFSET;
      const int D_IGEBETA=31+ORNL2016_OFFSET;
      const int D_INAIBETA=35+ORNL2016_OFFSET;
      const int D_IHAGBETA=45+ORNL2016_OFFSET;


     
    }
}//namespace dammIds

using namespace std;
using namespace dammIds::experiment;

void ORNL2016Processor::DeclarePlots(void) {
    
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
   
    DeclareHistogram1D(D_BETASCALARRATE,SB,"Beta scalar per cycle");
    DeclareHistogram1D(D_BETAENERGY,SD,"Beta Energy");
    
 
    // //Declaring beta gated
    // for (unsigned int i=0; i < 4; i++){
    //   stringstream ss;
    //   ss<< "HPGe " << i <<" - beta gated " ;
    //   DeclareHistogram1D(D_IGEBETA + i,SC,ss.str().c_str());     
    // }

    // //Declaring NaI beta gated
    // for (unsigned int i=0; i < 10; i++){
    //   stringstream ss;
    //   ss<< "NaI " << i << " - beta gated ";
    //   DeclareHistogram1D(D_INAIBETA + i ,SC,ss.str().c_str());     
    // }

    // //Declaring Hagrid beta gated
    // for (unsigned int i = 0; i < 16; i++){
    //   stringstream ss;
    //   ss<< "LaBr " << i << " - beta gated ";
    //   DeclareHistogram1D(D_IHAGBETA + i ,SC,ss.str().c_str());     

    // }   

}

void ORNL2016Processor::rootArrayreset(double arrayName[], int arraySize){ //refills arrayName with 0s 
  fill(arrayName,arrayName + arraySize,0);
}

void ORNL2016Processor::rootGstrutInit(RAY &strutName) { //Zeros the entire aux  structure

  fill(strutName.LaBr,strutName.LaBr + 16,0);
  fill(strutName.NaI,strutName.NaI + 10,0);
  fill(strutName.Ge,strutName.Ge + 4,0);
  strutName.beta = -9999;
  strutName.cycle=-9999;
  strutName.eventNum=-9999;
  strutName.gMulti=-9999;
  strutName.nMulti=-9999;
  strutName.lMulti=-9999;
  strutName.bMulti=-9999;

}
void ORNL2016Processor::rootGstrutInit2(adbk &strutName) { //Zeros the entire pab  structure
  strutName.LabE=-999;
  strutName.NabE=-999;
  strutName.GabE=-999;
  strutName.LabEvtNum=-9999;
  strutName.NabEvtNum=-9999;
  strutName.GabEvtNum=-9999;
  strutName.Lmulti=-999;
  strutName.Nmulti=-999;  
  strutName.Gmulti=-999;
}

ORNL2016Processor::ORNL2016Processor(double gamma_threshold_L, double sub_event_L, double gamma_threshold_N, double sub_event_N, double gamma_threshold_G, double sub_event_G) :EventProcessor(OFFSET,RANGE,"ORNL2016Processor"){

 

  associatedTypes.insert("ge");
  associatedTypes.insert("nai");
  associatedTypes.insert("labr3");
  associatedTypes.insert("beta");

  LgammaThreshold_ = gamma_threshold_L;
  LsubEventWindow_ = sub_event_L;
  NgammaThreshold_ = gamma_threshold_N;
  NsubEventWindow_ = sub_event_N;
  GgammaThreshold_ = gamma_threshold_G;
  GsubEventWindow_ = sub_event_G;

  // //initalize addback vecs
  LaddBack_.push_back(ScintAddBack());
  NaddBack_.push_back(ScintAddBack());
  GaddBack_.push_back(ScintAddBack());




  char hisFileName[32];
  GetArgument(1,hisFileName,32);
  string tmp = hisFileName;
  tmp= tmp.substr(0,tmp.find_first_of(" "));
  stringstream rootname;
  rootname<<tmp<<".root";
  rootFName_ =  new TFile(rootname.str().c_str(),"RECREATE");
  Taux = new TTree("Taux","Tree containing ancillary detector events with cycle");
  auxBranch = Taux->Branch("aux",&aux,"LaBr[16]/D:NaI[10]/D:Ge[4]/D:beta/D:eventNum/D:cycle/i:gMulti/i:nMulti/i:hMulti/i:bMulti/i");
  Taux->SetAutoFlush(3000);
  Tab = new TTree("Tab","Tree containing Processor add back information");
  adbkBranch = Tab->Branch("pab",&pab,"LabE/D:NabE/D:GabE/D:LabEvtNum/D:NabEvtNum/D:GabEvtNum/D:Lmulti/i:Nmulti/i:Gmulti/i");
  Tab->SetAutoFlush(3000);

  
    rootGstrutInit(aux);
    rootGstrutInit2(pab);
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
  
    map<unsigned int, pair<double,double> > lrtBetas;
    bool hasBeta = TreeCorrelator::get()->place("Beta")->status(); //might need a static initalize to false + reset at the end

   if(event.GetSummary("beta:double")->GetList().size() != 0) {
     lrtBetas = ((DoubleBetaProcessor*)DetectorDriver::get()->
		 GetProcessor("DoubleBetaProcessor"))->GetLowResBars();
   }
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
	if(hasBeta)
	plot(D_LABR3BETA, (*it)->GetCalEnergy());

    }
    
    
    ///PLOTTING THE SINGLES AND BETA GATED SPECTRA FOR NAI
    for(vector<ChanEvent*>::const_iterator naiIt = naiEvts.begin(); 
     	naiIt != naiEvts.end(); naiIt++) {
	plot(D_NAISUM, (*naiIt)->GetCalEnergy());
	if (hasBeta)
	  plot(D_NAIBETA, (*naiIt)->GetCalEnergy());
    }
    
   
   /// PLOT ANALYSIS HISTOGRAMS-------------------------------------------------------------------------------------------------------------------------------------

    rootGstrutInit(aux);
    rootGstrutInit2(pab);
      
    //Seting vars for addback
    double LrefTime = -2.0 * LsubEventWindow_;
    double NrefTime = -2.0 * NsubEventWindow_;
    double GrefTime = -2.0 * GsubEventWindow_;

    


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
  aux.cycle = cycleNum;
    
  //Multiplicitys
  aux.gMulti=geEvts.size();
  aux.nMulti=naiEvts.size();
  aux.lMulti=labr3Evts.size();
  aux.bMulti=lrtBetas.size();


  //Betas
  for(map<unsigned int, pair<double,double> >::iterator bIt = lrtBetas.begin();
      bIt != lrtBetas.end(); bIt++){
    plot(D_BETASCALARRATE,cycleNum );//PLOTTING BETA SCALAR SUM per CYCLE (LIKE 759 but per cycle vs per second
    plot(D_BETAENERGY,bIt->second.second);
    aux.beta = (bIt->second.second);
  }
    
  //NaI    
  for(vector<ChanEvent*>::const_iterator itNai = naiEvts.begin();
      itNai != naiEvts.end(); itNai++) {
    int nainum= (*itNai)->GetChanID().GetLocation();
    aux.NaI[nainum] = (*itNai)->GetCalEnergy();
        
    if (hasBeta){  //Beta Gated to Remove La Contamination
      //begin addback calulations for NaI

      double energy = (*itNai)->GetCalEnergy();
      double time = (*itNai)->GetCorrectedTime();
      
      if (energy < NgammaThreshold_){
	continue;
      }//end energy comp if statment
      
      double t1 = Globals::get()->clockInSeconds();
      double dtime = abs(time-NrefTime)*t1;
      
      if (dtime > GsubEventWindow_){ //if event time is outside sub event window start new addback after filling tree
	
	pab.NabEvtNum = evtNum;
	pab.NabE = NaddBack_.back().energy;
	pab.Nmulti = NaddBack_.back().multiplicity;
	Tab->Fill();
	NaddBack_.push_back(ScintAddBack());
      }
      
      NaddBack_.back().energy += energy;
      NaddBack_.back().time = time;
      NaddBack_.back().multiplicity +=1;
      NrefTime = time;
    }
       
  } //NaI loop End

  //HPGe
  for(vector<ChanEvent*>::const_iterator itGe = geEvts.begin();
      itGe != geEvts.end(); itGe++) {
    int genum = (*itGe)->GetChanID().GetLocation();
    aux.Ge[genum] = (*itGe)->GetCalEnergy();
    
    if (hasBeta){ //betagated addback to cut LaBr contamination out
           
     //begin addback calulations for clover
    //ChanEvent *chan =(*itGe);
    double energy = (*itGe)->GetCalEnergy();
    double time = (*itGe)->GetCorrectedTime();
    
    if (energy < GgammaThreshold_){
      continue;
    }//end energy comp if statment
    
    double t1 = Globals::get()->clockInSeconds();
    double dtime = abs(time-GrefTime)*t1;
       
    if (dtime > GsubEventWindow_){ //if event time is outside sub event window start new addback after filling tree
      
      pab.GabEvtNum = evtNum;
      pab.GabE = GaddBack_.back().energy;
      pab.Gmulti = GaddBack_.back().multiplicity;
      Tab->Fill();
      GaddBack_.push_back(ScintAddBack());
    }
    
    GaddBack_.back().energy += energy;
    GaddBack_.back().time = time;
    GaddBack_.back().multiplicity +=1;
    GrefTime = time;
    }
 } //GE loop end

    //Hagrid 
    for(vector<ChanEvent*>::const_iterator itLabr = labr3Evts.begin();
	itLabr != labr3Evts.end(); itLabr++){
      int LaBrnum = (*itLabr)->GetChanID().GetLocation();
      
      if (hasBeta){
	//begin addback calulations for LaBr Beta Gated to Remove La Contamination
	
	double energy = (*itLabr)->GetCalEnergy();
	double time = (*itLabr)->GetCorrectedTime();
	
	if (energy < LgammaThreshold_){
	  continue;
	}//end energy comp if statment
	
      double t1 = Globals::get()->clockInSeconds();
      double dtime = abs(time-NrefTime)*t1;
      
      if (dtime > GsubEventWindow_){ //if event time is outside sub event window start new addback after filling tree
	
	pab.LabEvtNum = evtNum;
	pab.LabE = LaddBack_.back().energy;
	pab.Lmulti = LaddBack_.back().multiplicity;
	Tab->Fill();
	LaddBack_.push_back(ScintAddBack());
      }// end if for new entry in vector
      
      LaddBack_.back().energy += energy;
      LaddBack_.back().time = time;
      LaddBack_.back().multiplicity +=1;
      LrefTime = time;


      }//end beta gate


      aux.LaBr[LaBrnum]= (*itLabr)->GetCalEnergy();
    } //Hagrid loop end	  
    
   





    aux.eventNum=evtNum;
    

    Taux->Fill();      
    
    evtNum++;
    EndProcess();
    return(true);
}
