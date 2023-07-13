// This program runs the charm hadrons decay simulation
// keywords: hadronic decay
// feeding a single particle gun at given energy
// and lifetime to get decay products
//
//

#include "Pythia8/Pythia.h"
//#include "Pythia8Plugins/HepMC3.h"
#include <iostream>
#include <fstream>
#include <string>
#include<cmath> 

//for generator
#include<vector>
#include<algorithm>
#include<iterator>
using namespace Pythia8;

//==========================================================================
//==========================================================================

// Single-particle gun. The particle must be a colour singlet.
// Input: flavour, energy, direction (theta, phi).
// If theta < 0 then random choice over solid angle.
// Optional final argument to put particle at rest => E = m.

void fillParticle(int id, double ee, double thetaIn, double phiIn,
  Event& event, ParticleData& pdt, Rndm& rndm, bool atRest = false,
  bool hasLifetime = true) {

  // Reset event record to allow for new event.
  event.reset();

  // Select particle mass; where relevant according to Breit-Wigner.
  double mm = pdt.mSel(id);
  double pp = sqrtpos(ee*ee - mm*mm);
  // Special case when particle is supposed to be at rest.
  if (atRest) {
    ee = mm;
    pp = 0.;
  }

  // Angles as input or uniform in solid angle.
  double cThe, sThe, phi;
  if (thetaIn >= 0.) {
    cThe = cos(thetaIn);
    sThe = sin(thetaIn);
    phi  = phiIn;
  } else {
    cThe = 2. * rndm.flat() - 1.;
    sThe = sqrtpos(1. - cThe * cThe);
    phi = 2. * M_PI * rndm.flat();
  }

  // Store the particle in the event record.
  int iNew = event.append( id, 1, 0, 0, pp * sThe * cos(phi),
    pp * sThe * sin(phi), pp * cThe, ee, mm);

  // Generate lifetime, to give decay away from primary vertex.
  if (hasLifetime) event[iNew].tau( event[iNew].tau0() * rndm.exp() );

};

// following function generates pythia events (decays from particle gun)
// and extracts the muon information (multiplicity, energy fraction and angle theta)

void run_get_events(int pdgid, double energy, int nevents, string filename) {
//I/O stuff
//initialize simple text file to store event info
  std::ofstream myfile;
  myfile.open(filename);
  myfile << "#Columns: Muon multiplicity, energy fraction, theta, Emin: 10GeV, Emax=100PeV, gridpoints = 350" << endl;

	//pdg id
	int idGun        = pdgid;
	//energy
	double eeGun     = energy;
	// if the particle is at rest
	bool atRest      = false;
	// if the particle has non-zero lifetime
	bool hasLifetime = true;
	// number of events to generate
	int nEvent       = nevents;
	// number of events to list
	int nList        = 3;

	//Generator; shorthand for event and particleData.
	Pythia pythia;
	Event& event      = pythia.event;
	ParticleData& pdt = pythia.particleData;

	pythia.readString("Print:quiet = on");
	// Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
	pythia.readString("ProcessLevel:all = off");
	pythia.readString("PartonLevel:all = off");
	//pythia.readString("HadronLevel:all = off");
	pythia.readString("HadronLevel:Hadronize = off");
	pythia.readString("HadronLevel:Decay = on");

	bool showScaleAndVertex = true;

	//turn ON all the charm hadron decays explicitely
	pythia.readString("411:mayDecay = true");
	pythia.readString("421:mayDecay = true");
	pythia.readString("431:mayDecay = true");
	pythia.readString("4122:mayDecay = true");
	pythia.readString("4132:mayDecay = true");
	pythia.readString("4232:mayDecay = true");
	pythia.readString("4332:mayDecay = true");

	pythia.init();

// Book histogram (for testing purpose ONLY)
//  Hist muMult("Muon multiplicity", 10, -0.5, 9.5);
//  Hist muFrac("Muon energy fraction", 100, 0.0, 1.0);
//  Hist muTheta("Muon theta", 100, 0.0, M_PI);

  // Begin of event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    // Set up single particle, along the z direction (theta=0, phi=0)
	fillParticle( idGun, eeGun, 0., 0.,
      event, pdt, pythia.rndm, atRest, hasLifetime);

    // Generate events. Quit if failure.
    if (!pythia.next()) {
    //if (!pythia.forceHadronLevel()) {
    //cout << " Event generation aborted prematurely, owing to error!\n";
    //break;
      continue;
    }

    // List first few events. (for testing purpose only)
//    if (iEvent < nList) {
//      event.list(showScaleAndVertex);
      // Also list junctions.
//      event.listJunctions();
//    }

	// get statistical variables
	int mult = 0;// muon multiplicity per event
	double frac = 0.0;// muon energy fraction
	double theta = 0.0;// muon opening angle

	for (int i=0; i < event.size(); i++) {
	// get the particle adsolute pdg id
	  int pid = event[i].idAbs();
	  // check for muon in the decay products
	  if (pid == 13) {
	  	++mult;
		// get the variables temporarily (in case there are >1 muon)
		double tfrac = event[i].e()/eeGun;
		double ttheta = event[i].theta();
		// select the variables for the highest energy muons
		if (tfrac > frac) {
			frac = tfrac;
			theta = ttheta;
		}
	  }
	}
	//filling histogram (testing only)
//	if (mult>0) {
//		muFrac.fill(frac);
//		muTheta.fill(theta);
//	}
//	muMult.fill(mult);

	// filling the text file
	myfile << mult << ',' << frac << ',' << theta << endl;


  //end of event loop
  }

  pythia.stat();

  //printing the histograms (testing only)
//    cout << muMult << muFrac << muTheta;
  //Done.
  myfile.close();
}


//following class creats the numpy logspace like arrays
template<typename T = double>
class Logspace {
private:
    T curValue, base, step;

public:
    Logspace(T first, T last, int num, T base = 10.0) : curValue(first), base(base){
       step = (last - first)/(num-1);
    }

    T operator()() {
        T retval = pow(base, curValue);
        curValue += step;
        return retval;
    }
};

int main(int argc, char* argv[]) {
	//Get the command line variables
	//particle pdg code
	int pid = std::stoi(argv[1]);
	//outputfilename
	string fname = argv[2];
	// number of event sto generate per particle energies
	int neve = std::stoi(argv[3]);
	// number of grid points to loop over
	int ngrid = std::stoi(argv[4]);
	// setting the minimum and maximum energy power index
	double emin_pow = 1.0;//10^1GeV
	double emax_pow = 8.0;//10^8GeV
	//
	//for DEBUG ONLY ---->
	//double emin_pow = 8.0;//10^1GeV
	//double emax_pow = 8.0;//10^8GeV
	//------------------->
	// get the energy grid point array
	std::vector<double> enarr;
	std::generate_n(std::back_inserter(enarr), ngrid, Logspace<>(emin_pow,emax_pow,ngrid,10));

	cout << "Starting event generation for particle decay of: "
		<< pid << "  ;Minimum energy: " << pow(10, emin_pow)
		<< "   ;Maximum energy: " << pow(10, emax_pow)
		<< "   ;NGridPoints: " << ngrid <<
		"   ; NEvents: " << neve << endl;
	// loop over the energy points and store each run into iterative files
 	for (int gridpoint=0; gridpoint<ngrid; gridpoint++) {
		// get the new filename
		string newf = fname;
		int rpos = newf.find(".txt");
		string tag = "_" + std::to_string(gridpoint) + ".txt";
		newf.replace(rpos, 4, tag);
		cout << "Generating decay for: "<< pid << "\t" 
			<< enarr[gridpoint]<< "GeV" << "\t" << newf << endl;
		run_get_events(pid, enarr[gridpoint], neve, newf);
	}

      
return 0;
}
