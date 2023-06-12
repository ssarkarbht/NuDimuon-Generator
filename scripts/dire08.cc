
// DIRE includes.
#include "Dire/Dire.h"

// Pythia includes.
#include "Pythia8/Pythia.h"

// Standard Libraries
#include <iostream>
#include <fstream>
#include <string>
#include<cmath>
#include <sstream>
#include <vector>
#include <iomanip>

using namespace Pythia8;

//==========================================================================
//following function checks if the given particle is a charm hadron
bool isCharmHadron(Particle& part){
	//if ((part.status()>80 && part.status()<92) &&
	if (part.status()>0 &&
	   ((part.idAbs()>400 && part.idAbs()<500) ||
	    (part.idAbs()>4000 && part.idAbs()<5000)))
	  return true;
     return false;
}

//==========================================================================
// following function runs pythia+dire for every neutrino
// energy instance
void run_pythia(double nuen, int nutype, int target_type,
	       	string settings, std::ofstream& outfile){
//void run_pythia(std::ofstream& outfile){
  Pythia pythia;

  // Create and initialize DIRE shower plugin.
  Dire dire;
  dire.init(pythia, settings.c_str());

  // Flag to check for charm event
  bool charmFound = false;

  //overwrite the neutrino energy and target type with
  //values from the input event configuration file
  pythia.settings.parm("Beams:eA", nuen);
  pythia.settings.mode("Beams:idA", nutype);
  //target 1:proton, 2:neutron
  if (target_type==1) pythia.settings.mode("Beams:idB", 2212);
  if (target_type==2) pythia.settings.mode("Beams:idB", 2112);

  //use quiet mode (don't print anything unless debugging or logging)
  //pythia.readString("Print:quiet = on");

  //re-intialize the pythia settings (after changing few settigns)
  pythia.init();

  // Start generation loop
  int nEvent = pythia.settings.mode("Main:numberOfEvents");
  for( int iEvent=0; iEvent<nEvent; ++iEvent ){

    // Generate next event
    if( !pythia.next() ) {
      if( pythia.info.atEndOfFile() )
        break;
      else continue;
    }

    // Get event weight(s).
    double evtweight         = pythia.info.weight();

    // Do not print zero-weight events.
    if ( evtweight == 0. ) continue;

    // Retrieve the shower weight.
    dire.weightsPtr->calcWeight(0.);
    dire.weightsPtr->reset();
    double wt = dire.weightsPtr->getShowerWeight();

    if (abs(wt) > 1e4) {
        cout << "Warning in DIRE main program dire06.cc: Shower weight larger"
        << " than 10000. Discarding event with rare shower weight fluctuation."
        << endl;
        evtweight = 0.;
      }
    // Do not print zero-weight events.
    if ( evtweight == 0. ) continue;

    //Event Analysis (check for charm quark production)
    for (int i=0; i < pythia.event.size(); i++){
	if (pythia.event[i].status()==-23 &&
	   pythia.event[i].idAbs()==4) {
             charmFound=true;
	     //For testing only -->
	     //pythia.event.list();
             break;
        }
    } // end of particle loop

    //list the event with charm production
    if (charmFound){
	//Testing only (commented)--------->
	//pythia.event.list();
	//keep track of the number of particles to store
	int npart = 0;
	//loop over all the particles in the event
	for (int i=0; i < pythia.event.size(); i++){
	  Particle part = pythia.event[i];
	  //get the outgoing lepton
	  if (part.status()==23 &&
		part.idAbs()>10 &&
		part.idAbs()<17){
		//outgoing lepton found, update the particle number
		npart++;
		//write the output to file
		outfile << "P" << "\t" << npart << "\t" <<
		part.id() << "\t"<<std::fixed << setprecision(8)<<
		part.e() << "\t" << part.theta() << "\t"
		<< part.phi() << endl;
          }
	  //get the outgoing quark info
	  if (part.status()==-23 &&
		part.idAbs()==4){
		//charm quark found, uodate the particle number
		npart++;
		//write the output to file
		outfile << "P" << "\t" << npart << "\t" <<
		part.id() << "\t"<< std::fixed << setprecision(8)<<
		part.e() << "\t" << part.theta() << "\t"
		<< part.phi() << endl;
          }
	  //get the outgoing charm hadrons
	  if (isCharmHadron(part)){
		//charm hadron found, uodate the particle number
		npart++;
		//write the output to file
		outfile << "P" << "\t" << npart << "\t" <<
		part.id() << "\t"<< std::fixed << setprecision(8)<<
		part.e() << "\t" << part.theta() << "\t"
		<< part.phi() << endl;
	  }
	}
	
        break;//end the pythia run once th3e event is extracted
    }//end of charm event loop

  } // end loop over events to generate

  // print cross section, errors
  //pythia.stat();

}

//==========================================================================

int main( int argc, char* argv[]  ){

  // Check that correct number of command-line arguments
  if (argc < 4) {
    cerr << " Unexpected number of command-line arguments ("<<argc-1<<"). \n"
         << " You are expected to provide the arguments" << endl
         << " 1. Input file for settings" << endl
         << " 2. Input file for Neutrino energies and target (p/n)" << endl
         << " 2. Output event text file" << endl
         << argc-1 << " arguments provided:";
         for ( int i=1; i<argc; ++i) cerr << " " << argv[i];
         cerr << "\n Program stopped. " << endl;
    return 1;
  }
  string f_settings = argv[1];
  string f_points = argv[2];
  string f_out = argv[3];

  // I/O  settings
  
  // Set the output file settings
  std::ofstream outfile;
  outfile.open(f_out.c_str());

  //Read and loop over event points file
  std::ifstream energy_input;
  energy_input.open(f_points.c_str());

  vector<string> temp;
  string line;//collector of each row
  string val;//collector of each field in a row

  while(getline(energy_input, line, '\n')){
    //get the data from each line
    temp.clear();
    stringstream s(line);
    while(getline(s, val, ',')){\
	temp.push_back(val);
    }

    stringstream evenum;
    evenum << temp[0];
    int evtNum;
    evenum >> evtNum;

    stringstream ntype;
    ntype << temp[1];
    int nuType;
    ntype >> nuType;

    stringstream nuen;
    nuen << temp[2];
    double nuEnergy;
    nuen >> nuEnergy;

    stringstream ttype;
    ttype << temp[3];
    int targetType;
    ttype >> targetType;

    outfile << "E" << "\t" << nuType << "\t" << evtNum << "\t" <<
    std::fixed << setprecision(8)<< nuEnergy << "\t" << targetType << endl;

	run_pythia(nuEnergy, nuType, targetType, f_settings, outfile);
	//run_pythia(outfile);

  }//end of event points file loop
  // Done
  outfile.close();
  energy_input.close();
  return 0;

}
