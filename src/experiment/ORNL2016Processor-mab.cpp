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

#include "ORNL2016Processor.hpp"
#include "DoubleBetaProcessor.hpp"
#include "GeProcessor.hpp"

#include "TH1.h"

static unsigned int evtNum = 0;

namespace dammIds {
    namespace experiment {
        const unsigned int ORNL2016_OFFSET = 0;
        const int D_VANDLEMULT = 0 + ORNL2016_OFFSET;

        const int DD_CTOFNOTAPE = 1 + ORNL2016_OFFSET;
        const int DD_QDCTOFNOGATE = 2 + ORNL2016_OFFSET;
        const int DD_CORTOFVSEGAM = 3 + ORNL2016_OFFSET;
        const int DD_QDCVSCORTOFMULT1 = 4 + ORNL2016_OFFSET;
        const int DD_MULT2SYM = 5 + ORNL2016_OFFSET;
        const int DD_LIGLEN = 6 + ORNL2016_OFFSET;
        const int DD_LIGLTOF = 7 + ORNL2016_OFFSET;

        const int D_LABR3SUM = 8 + ORNL2016_OFFSET;
        const int D_LABR3BETA = 9 + ORNL2016_OFFSET;

        const int D_NAISUM = 10 + ORNL2016_OFFSET;
        const int D_NAIBETA = 11 + ORNL2016_OFFSET;

        const int DD_TOFVSNAI = 12 + ORNL2016_OFFSET;
        const int DD_TOFVSHAGRID = 13 + ORNL2016_OFFSET;
        const int DD_TOFVSGE = 14 + ORNL2016_OFFSET;

        //seeds for Damm cycle his
        const int D_BETASCALARRATE = 29 + ORNL2016_OFFSET; //6079 in his
        const int D_BETAENERGY = 30 + ORNL2016_OFFSET;
        const int D_IGEBETA = 31 + ORNL2016_OFFSET;
        const int D_INAIBETA = 35 + ORNL2016_OFFSET;
        const int D_IHAGBETA = 45 + ORNL2016_OFFSET;


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

    DeclareHistogram1D(D_BETASCALARRATE, SB, "Beta scalar per cycle");
    DeclareHistogram1D(D_BETAENERGY, SD, "Beta Energy");


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


void ORNL2016Processor::rootGstrutInit(RAY &strutName) { //Zeros the entire aux  structure

    fill(strutName.LaBr, strutName.LaBr + 16, 0);
    fill(strutName.NaI, strutName.NaI + 10, 0);
    fill(strutName.Ge, strutName.Ge + 4, 0);
    strutName.beta = -9999;
    strutName.cycle = -9999;
    strutName.eventNum = -9999;
    strutName.gMulti = -9999;
    strutName.nMulti = -9999;
    strutName.lMulti = -9999;
    strutName.bMulti = -9999;

}

void ORNL2016Processor::rootGstrutInit2(PROSS &strutName) { //Zeros the entire processed structure
    strutName.AbE = -999;
    strutName.AbEvtNum = -9999;
    strutName.Multi = 0;
    strutName.SymX = -999;
    strutName.SymY = -999;
}

ORNL2016Processor::ORNL2016Processor(double gamma_threshold_L, double sub_event_L, double gamma_threshold_N,
                                     double sub_event_N, double gamma_threshold_G, double sub_event_G) : EventProcessor(
        OFFSET, RANGE, "ORNL2016Processor") {


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
    LaddBack_.push_back(ScintAddBack(0, 0, 0));
    NaddBack_.push_back(ScintAddBack(0, 0, 0));
    GaddBack_.push_back(ScintAddBack(0, 0, 0));


    char hisFileName[32];
    GetArgument(1, hisFileName, 32);
    string tmp = hisFileName;
    tmp = tmp.substr(0, tmp.find_first_of(" "));
    stringstream rootname;
    rootname << tmp << ".root";
    rootFName_ = new TFile(rootname.str().c_str(), "RECREATE");
    Taux = new TTree("Taux", "Tree containing ancillary detector events with cycle");
    singBranch = Taux->Branch("sing", &sing,
                             "LaBr[16]/D:NaI[10]/D:Ge[4]/D:beta/D:eventNum/D:cycle/i:gMulti/i:nMulti/i:hMulti/i:bMulti/i");

    gProcBranch = Taux->Branch("Gpro", &Gpro,"AbE/D:AbEvtNum/D:Multi/D:SymX/D:SymY/D");
    lProcBranch = Taux->Branch("Lpro", &Lpro,"AbE/D:AbEvtNum/D:Multi/D:SymX/D:SymY/D");
    nProcBranch = Taux->Branch("Npro", &Npro,"AbE/D:AbEvtNum/D:Multi/D:SymX/D:SymY/D");

    Taux->SetAutoFlush(3000);
    rootGstrutInit(sing);
    rootGstrutInit2(Gpro);
    rootGstrutInit2(Lpro);
    rootGstrutInit2(Npro);
}


ORNL2016Processor::~ORNL2016Processor() {

    rootFName_->Write();
    rootFName_->Close();

    delete (rootFName_);


}


bool ORNL2016Processor::PreProcess(RawEvent &event) {
    if (!EventProcessor::PreProcess(event))
        return (false);


    EndProcess();
    return (true);
}

bool ORNL2016Processor::Process(RawEvent &event) {
    if (!EventProcessor::Process(event))
        return (false);

    map<unsigned int, pair<double, double> > lrtBetas;
    bool hasBeta = TreeCorrelator::get()->place(
            "Beta")->status(); //might need a static initalize to false + reset at the end

    if (event.GetSummary("beta:double")->GetList().size() != 0) {
        lrtBetas = ((DoubleBetaProcessor *) DetectorDriver::get()->
                GetProcessor("DoubleBetaProcessor"))->GetLowResBars();
    }
    static const vector<ChanEvent *> &labr3Evts =
            event.GetSummary("labr3")->GetList();
    static const vector<ChanEvent *> &naiEvts =
            event.GetSummary("nai")->GetList();
    static const vector<ChanEvent *> &geEvts =
            event.GetSummary("ge")->GetList();


    /// PLOT ANALYSIS HISTOGRAMS-------------------------------------------------------------------------------------------------------------------------------------

    rootGstrutInit(sing); // initalize the root structures
    rootGstrutInit2(Gpro);
    rootGstrutInit2(Lpro);
    rootGstrutInit2(Npro);

    //Seting vars for addback
    double LrefTime = -2.0 * LsubEventWindow_;
    double NrefTime = -2.0 * NsubEventWindow_;
    double GrefTime = -2.0 * GsubEventWindow_;




    //Cycle timing
    static double cycleLast = 2;
    static int cycleNum = 0;
    if (TreeCorrelator::get()->place("Cycle")->status()) {
        double cycleTime = TreeCorrelator::get()->place("Cycle")->last().time;
        cycleTime *= (Globals::get()->clockInSeconds() * 1.e9);
        if (cycleTime != cycleLast) {
            double tdiff = (cycleTime - cycleLast) / 1e6; //Outputs cycle length in msecs.
            if (cycleNum == 0) {
                cout
                        << " #  There are some events at the beginning of the first segment missing from Histograms that use cycleNum."
                        << endl << " #  This is a product of not starting the cycle After the LDF." << endl
                        << " #  This First TDIFF is most likly nonsense" << endl;
            }
            cycleLast = cycleTime;
            cycleNum++;
            cout << "Cycle Change " << endl << "Tdiff (Cycle start and Now) (ms)= " << tdiff << endl
                 << "Starting on Cycle #" << cycleNum << endl;
        }
    }
    sing.cycle = cycleNum;

    //set multiplicys for sing branch based on the size of the detector maps for the event. limitation: sub event is smaller than full event this will end up being too large
    sing.gMulti = geEvts.size();
    sing.nMulti = naiEvts.size();
    sing.lMulti = labr3Evts.size();
    sing.bMulti = lrtBetas.size();


    //Betas------------------------------------------------------------------------------------------------------------
    for (map<unsigned int, pair<double, double> >::iterator bIt = lrtBetas.begin();
         bIt != lrtBetas.end(); bIt++) {
        plot(D_BETASCALARRATE, cycleNum);//PLOTTING BETA SCALAR SUM per CYCLE (LIKE 759 but per cycle vs per second
        plot(D_BETAENERGY, bIt->second.second);
        sing.beta = (bIt->second.second);
    }

    //NaI ----------------------------------------------------------------------------------------------------------------------------------------------
    for (vector<ChanEvent *>::const_iterator itNai = naiEvts.begin();
         itNai != naiEvts.end(); itNai++) {
        int naiNum = (*itNai)->GetChanID().GetLocation();
        sing.NaI[naiNum] = (*itNai)->GetCalEnergy();
        plot(D_NAISUM, (*itNai)->GetCalEnergy()); //plot totals

        if (hasBeta) {  //Beta Gate
            plot(D_NAIBETA, (*itNai)->GetCalEnergy()); //plot beta-gated totals

            //begin addback calulations for NaI
            double energy = (*itNai)->GetCalEnergy();
            double time = (*itNai)->GetCorrectedTime();

            if (energy < NgammaThreshold_) {
                continue;
            }//end energy comp if statment
            double t1 = Globals::get()->clockInSeconds();
            double dtime = abs(time - NrefTime) * t1;

            if (dtime >
                NsubEventWindow_) { //if event time is outside sub event window start new addback after filling tree
                Npro.AbEvtNum = evtNum;
                Npro.AbE = NaddBack_.back().energy;
                Npro.Multi= NaddBack_.back().multiplicity;
                Taux->Fill();
                NaddBack_.push_back(ScintAddBack());
            }//end subEvent IF
            NaddBack_.back().energy += energy; // if still inside sub window: incrament
            NaddBack_.back().time = time;
            NaddBack_.back().multiplicity += 1;
            NrefTime = time;

            // Begin Symplot inner loop
            for (vector<ChanEvent *>::const_iterator itNai2 = itNai+1;
                 itNai2 != geEvts.end(); itNai2++) {
                double energy2=(*itNai2)->GetCalEnergy();
                int naiNum2=(*itNai2)->GetChanID().GetLocation();
                //double time2=(*itGe2)->GetCorrectedTime();
                if (energy2 < NgammaThreshold_ ) {
                    continue;
                }//end energy comp if statement
                if (naiNum2 != naiNum) {
                    Npro.SymX=energy;
                    Npro.SymY=energy2;
                    Taux->Fill();
                }


            } //end symplot inner loop
        }//end beta gate
    } //NaI loop End

    //HPGe CLOVER--------------------------------------------------------------------------------------------------------------------------------------
    for (vector<ChanEvent *>::const_iterator itGe = geEvts.begin();
         itGe != geEvts.end(); itGe++) {
        int geNum = (*itGe)->GetChanID().GetLocation();
        sing.Ge[geNum] = (*itGe)->GetCalEnergy();

        if (hasBeta) { //beta-gated Processing to cut LaBr contamination out

            //begin addback calulations for clover
            double energy = (*itGe)->GetCalEnergy();
            double time = (*itGe)->GetCorrectedTime();
            if (energy < GgammaThreshold_) {
                continue;
            }//end energy comp if statment
            double t1 = Globals::get()->clockInSeconds();
            double dtime = abs(time - GrefTime) * t1;

            if (dtime >
                GsubEventWindow_) { //if event time is outside sub event window start new addback after filling tree
                Gpro.AbEvtNum= evtNum;
                Gpro.Multi= GaddBack_.back().multiplicity;
                Gpro.AbE= GaddBack_.back().energy;
                Taux->Fill();
                GaddBack_.push_back(ScintAddBack());
            } //end subEvent IF

            GaddBack_.back().energy += energy;
            GaddBack_.back().time = time;
            GaddBack_.back().multiplicity += 1;
            GrefTime = time;

            // Begin Symplot inner loop
            for (vector<ChanEvent *>::const_iterator itGe2 = itGe+1;
                 itGe2 != geEvts.end(); itGe2++) {
                double energy2=(*itGe2)->GetCalEnergy();
                int geNum2=(*itGe2)->GetChanID().GetLocation();
                //double time2=(*itGe2)->GetCorrectedTime();
                if (energy2 < GgammaThreshold_ ) {
                    continue;
                }//end energy comp if statement
                if (geNum2 != geNum) {
                    Gpro.SymX=energy;
                    Gpro.SymY=energy2;
                    Taux->Fill();
                }


            } //end symplot inner loop




        } //end BetaGate
    } //GE loop end

    //Hagrid  ----------------------------------------------------------------------------------------------------------------------------------------------
    for (vector<ChanEvent *>::const_iterator itLabr = labr3Evts.begin();
         itLabr != labr3Evts.end(); itLabr++) {
        int labrNum = (*itLabr)->GetChanID().GetLocation();
        plot(D_LABR3SUM, (*itLabr)->GetCalEnergy()); //plot non-gated totals

        if (hasBeta) {

            plot(D_LABR3BETA, (*itLabr)->GetCalEnergy()); //plot beta-gated totals
            //begin addback calculations for LaBr | Beta Gated to Remove La Contamination

            double energy = (*itLabr)->GetCalEnergy();
            double time = (*itLabr)->GetCorrectedTime();

            if (energy < LgammaThreshold_) {
                continue;
            }//end energy comp if statment

            double t1 = Globals::get()->clockInSeconds();
            double dtime = abs(time - LrefTime) * t1;

            if (dtime >
                LsubEventWindow_) { //if event time is outside sub event window start new addback after filling tree

                Lpro.AbEvtNum= evtNum;
                Lpro.AbE= LaddBack_.back().energy;
                Lpro.Multi= LaddBack_.back().multiplicity;
                Taux->Fill();
                LaddBack_.push_back(ScintAddBack());
            }// end if for new entry in vector

            LaddBack_.back().energy += energy;
            LaddBack_.back().time = time;
            LaddBack_.back().multiplicity += 1;
            LrefTime = time;

            for (vector<ChanEvent *>::const_iterator itLabr2 = itLabr+1;
                 itLabr2 != geEvts.end(); itLabr2++) {
                double energy2=(*itLabr2)->GetCalEnergy();
                int labrNum2=(*itLabr2)->GetChanID().GetLocation();
                //double time2=(*itGe2)->GetCorrectedTime();
                if (energy2 < LgammaThreshold_ ) {
                    continue;
                }//end energy comp if statement
                if (labrNum2 != labrNum) {
                    Lpro.SymX=energy;
                    Lpro.SymY=energy2;
                    Taux->Fill();
                }


            } //end symplot inner loop

        }//end beta gate

        sing.LaBr[labrNum] = (*itLabr)->GetCalEnergy();
    } //Hagrid loop end


    sing.eventNum = evtNum;
    Taux->Fill();

    evtNum++;
    EndProcess();
    return (true);
}
