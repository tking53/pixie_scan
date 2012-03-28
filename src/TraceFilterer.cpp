/** \file  TraceFilterer.cpp
 *  \brief Implements the analysis of traces using trapezoidal filters
 *
 *  This trace plots the filtered traces as well as energy and timing info
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>

#include "damm_plotids.h"
#include "RandomPool.h"
#include "Trace.h"
#include "TraceFilterer.h"

using namespace std;

const string TraceFilterer::defaultFilterFile = "filter.txt";
const int TraceFilterer::energyBins = SC;
const double TraceFilterer::energyScaleFactor = 2.198; //< TO BE USED WITH MAGIC +40 ENERGY SAMPLE LOCATION
// const double TraceFilterer::energyScaleFactor = 2.547; //< multiply the energy filter sums by this to gain match to raw spectra

extern RandomPool randoms;

/** 
 *  A do nothing constructor
 */
TraceFilterer::PulseInfo::PulseInfo()
{
    isFound = false;
}

/**
 *  Constructor so we can build the pulseinfo in place
 */
TraceFilterer::PulseInfo::PulseInfo(Trace::size_type theTime, double theEnergy) :
    time(theTime), energy(theEnergy)
{
    isFound = true;
}

TraceFilterer::TraceFilterer() : 
    TracePlotter(), fastParms(5,5), 
    energyParms(100,100), thirdParms(20,10)
{
    //? this uses some legacy values for third parms, are they appropriate
    name = "Filterer";
    
    fastThreshold = fastParms.GetRiseSamples()  * 3;
    slowThreshold = thirdParms.GetRiseSamples() * 2;
    useThirdFilter = false;
}

TraceFilterer::~TraceFilterer()

{
    // do nothing
}

bool TraceFilterer::Init(const string &filterFile)
{
    const int maxTraceLength = 6400;

    TracePlotter::Init();

    fastFilter.reserve(maxTraceLength);
    energyFilter.reserve(maxTraceLength);
    thirdFilter.reserve(maxTraceLength);

    // read in the filter parameters
    ifstream in(filterFile.c_str());
    if (!in) {
	cout << "Failed to open the filter parameter file" << endl;
	cout << "  Using default values instead" << endl;
	return true;
    }

    while (!in.eof()) {
	if ( isdigit(in.peek()) ) {
	    Trace::size_type rise, gap;
	    
	    cout << "Trace analysis parameters are: " << endl;
	    in >> rise >> gap >> fastThreshold;
	    cout << "  Fast filter: " << rise << " " << gap 
		 << " " << fastThreshold << endl;
	    fastParms = TrapezoidalFilterParameters(gap, rise);

	    in >> rise >> gap;
	    cout << "  Energy filter: " << rise << " " << gap << endl;
	    energyParms = TrapezoidalFilterParameters(gap, rise);

	    in >> rise >> gap >> slowThreshold;	    
	    cout << "  Third filter: " << rise << " " << gap 
		 << " " << slowThreshold;
	    if (!useThirdFilter) {
		cout << ", not used";
	    }
	    cout << endl;
	    thirdParms = TrapezoidalFilterParameters(gap, rise);
	    
	    // scale thresholds by the length of integration (i.e. rise time)
	    fastThreshold *= fastParms.GetRiseSamples();
	    slowThreshold *= thirdParms.GetRiseSamples();
	} else {
	    // assume this is a comment
	    in.ignore(1000, '\n');
	}
    }

    if (in.fail()) {
	cerr << "Problem reading filter parameters file" << endl;
	return true;
    }

    return true;
}


void TraceFilterer::DeclarePlots(void) const
{
    using namespace dammIds::trace;

    TracePlotter::DeclarePlots();

    DeclareHistogram2D(DD_FILTER1, traceBins, numTraces, "fast filter");
    DeclareHistogram2D(DD_FILTER2, traceBins, numTraces, "energy filter");
    if (useThirdFilter) {
	DeclareHistogram2D(DD_FILTER3, traceBins, numTraces, "3rd filter");
    }
    DeclareHistogram2D(DD_REJECTED_TRACE, traceBins, numTraces, "rejected traces");

    DeclareHistogram1D(D_ENERGY1, energyBins, "E1 from trace"); 
}

void TraceFilterer::Analyze(Trace &trace,
			    const string &type, const string &subtype)
{
    using namespace dammIds::trace;
    TracePlotter::Analyze(trace, type, subtype);

    if (level >= 5) {
	trace.DoBaseline(0, 20); //? make the number of bins a parameter for the filter
	if (trace.GetValue("sigmaBaseline") * fastParms.GetRiseSamples() > fastThreshold / 2) {
	    static int rejectedTraces = 0;
	    cout << "Rejected trace with bad baseline = " << trace.GetValue("sigmaBaseline") << endl;
	    trace.Plot(DD_REJECTED_TRACE, rejectedTraces++);
	    EndAnalyze(); // update timing
	    return;
	}

	fastFilter.clear();
	energyFilter.clear();

	// determine trace filters, these are trapezoidal filters characterized
	//   by a risetime and a gaptime and a range of the filter 
	trace.TrapezoidalFilter(fastFilter, fastParms);
	trace.TrapezoidalFilter(energyFilter, energyParms);

	if (useThirdFilter) {
	    thirdFilter.clear();
	    trace.TrapezoidalFilter(thirdFilter, thirdParms);
	}
	FindPulse(fastFilter.begin(), fastFilter.end());

	if (pulse.isFound) {
	    trace.SetValue("filterTime", (int)pulse.time);
	    trace.SetValue("filterEnergy", pulse.energy);
	}

	// now plot some stuff
	/* quick plot function still needs scaling and handling of negative numbers
	fastFilter.Plot(DD_FILTER1, numTracesAnalyzed);
	energyFilter.Plot(DD_FILTER2, numTracesAnalyzed);
	thirdFilter.Plot(DD_FILTER3, numTracesAnalyzed);
	*/
	for (Trace::size_type i = 0; i < trace.size(); i++) {
	    if (i < fastFilter.size())
		plot(DD_FILTER1, i, numTracesAnalyzed, 
		     abs(fastFilter[i]) / fastParms.GetRiseSamples() );
	    if (i < energyFilter.size())
		plot(DD_FILTER2, i, numTracesAnalyzed, 
		     abs(energyFilter[i]) / energyParms.GetRiseSamples() );
	    if (i < thirdFilter.size() && useThirdFilter)
		plot(DD_FILTER3, i, numTracesAnalyzed, 
		     abs(thirdFilter[i]) / thirdParms.GetRiseSamples() );

	}
	// calculated values at end of traces
	/*
	plot(DD_TRACE, trace.size() + 10, numTracesAnalyzed, energy);
	plot(DD_TRACE, trace.size() + 11, numTracesAnalyzed, time);
	*/

	plot(D_ENERGY1, pulse.energy);
    } // sufficient analysis level

    EndAnalyze(trace);
}

const TraceFilterer::PulseInfo& TraceFilterer::FindPulse(Trace::iterator begin, Trace::iterator end)
{
    // class to see if fast filter is above threshold
    static binder2nd< greater<Trace::value_type> > crossesThreshold
	(greater<Trace::value_type>(), fastThreshold);

    Trace::size_type sample;
    pulse.isFound = false;

    while (begin < end) {
	begin = find_if(begin, end, crossesThreshold);
	if (begin == end) {
	    break;
	}

	pulse.time    = begin - fastFilter.begin();
	pulse.isFound = true;

	// sample the slow filter in the middle of its size
	if (useThirdFilter) {
	    sample = pulse.time + (thirdParms.GetSize() - fastParms.GetSize()) / 2;
	    if (sample >= thirdFilter.size() ||
		thirdFilter[sample] < slowThreshold) {
		begin++; 
		continue;
	    }
	}
	//? some investigation needed here for good resolution
	// add a presample location
	// sample = pulse.time + (energyParms.GetSize() - fastParms.GetSize()) / 2;
	sample = pulse.time + 40;
	
	if (sample < energyFilter.size()) {
	    pulse.energy = energyFilter[sample] + randoms.Get();	    
	    // scale to the integration time
	    pulse.energy /= energyParms.GetRiseSamples();
	    pulse.energy *= energyScaleFactor;
	} else pulse.energy = NAN;
	break;
    }

    return pulse;
}
