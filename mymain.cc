// main02.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the pT_Z spectrum at the Tevatron.

#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "Pythia8Plugins/FastJet3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <set>
#include <stdlib.h> 
using namespace std;
using namespace Pythia8;
double deltaPHI(double phi1,double phi2){
	double R;
	R=phi1-phi2;
	if (R>M_PI){R=2*M_PI-R;
	}
	if (R<-M_PI){R=R+2*M_PI;
	}
	return R;
}
double deltaR(double phi1,double phi2,double eta1,double eta2){
	double R;
	R=sqrt(pow(deltaPHI(phi1,phi2),2)+pow(eta1-eta2,2));
	return R;
}
void print(std::vector<int> const &input)
{
	for (auto const& i: input) {
		std::cout << i << " ";
	}
}

int main() {
  // Number of events, generated and listed ones (for jets).
  int  nEvent = 30;
  int nListJets = nEvent; //this value should <= nEvent
  
  // Select common parameters for SlowJet and FastJet analyses.
  int    power   = -1;     // -1 = anti-kT; 0 = C/A; 1 = kT.
  double R       = 0.4;    // Jet size.
  double pTMin   = 500.0;    // Min jet pT.
  double pTMax   = 550.0;    // Max jet pT.(not applied)
  double etaMax  = 2.0;    // Pseudorapidity range of detector.
  double d_R     = 0.2;    //matching radius for tagging jet 
  int gg         = 21;
  set<int> qq;
  qq.insert(qq.begin(),1);
  qq.insert(qq.begin(),-1);
  qq.insert(qq.begin(),2);
  qq.insert(qq.begin(),-2);
  qq.insert(qq.begin(),3);
  qq.insert(qq.begin(),-3);
//  int    xjet    = 21;      //finding x jet (x=gluon or quark)
//  int    xjetp              //out put x jet nth vector component
  int    select  = 2;      // Which particles are included?
  int    massSet = 2;      // Which mass are they assumed to have?
  
  // Generator. Process selection. Tevatron initialization. Histogram.
  Pythia pythia;
  Event& event = pythia.event;
  
  pythia.readString("Beams:eCM = 14000.");
  pythia.readString("WeakBosonAndParton:qg2gmZq = on");  //"WeakBosonAndParton:qqbar2gmZg = on" for z+g,WeakBosonAndParton:qg2gmZq for z+q
  pythia.readString("PhaseSpace:pTHatMin = 500.");  //PhaseSpace:pTHatMax   (default = -1.)
  pythia.readString("PhaseSpace:pTHatMax = 550.");
 
  

  pythia.init();
//========================================================================================
  //open file
ofstream myfile;
  myfile.open ("myevents.txt");//qq for quark ,gg for gluon
ofstream myout;
  myout.open ("myout.txt");//qq for quark ,gg for gluon
//========================================================================================


  // Set up FastJet jet finder.
  //   one can use either explicitly use antikt, cambridge, etc., or
  //   just use genkt_algorithm with specification of power
  //fastjet::JetAlgorithm algorithm;
  //if (power == -1)      algorithm = fastjet::antikt_algorithm;
  //if (power ==  0)      algorithm = fastjet::cambridge_algorithm;
  //if (power ==  1)      algorithm = fastjet::kt_algorithm;
  //fastjet::JetDefinition jetDef(algorithm, R);
  // there's no need for a pointer to the jetDef (it's a fairly small object)
  fastjet::JetDefinition jetDef(fastjet::genkt_algorithm, R, power);
  std::vector <fastjet::PseudoJet> fjInputs;

 std::cout<<"===================================================================111111111111111111========================================================================================================"<<
		endl;

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    clock_t befGen = clock();
    if (!pythia.next()) continue;
    clock_t aftGen = clock();
    

   
   
 std::cout<<"===================================================================2222222222222222222========================================================================================================"<<
		endl;

    // Begin FastJet analysis: extract particles from event record.
    clock_t befFast = clock();
    fjInputs.resize(0);
    Vec4   pTemp;
    vector<double> phi0, eta0;
    double mTemp,ptmp,etamp,phimp,idmp;
    std::set<int>::iterator qq_index;
    int nAnalyze = 0;
     //Show the haedest event
     std::cout<<"============================hardest-event================================"<<endl;
     std::cout<<"no.\teventid\tstatus\tmothers\t\tdaughters\tpT\teta\t\phi\n"<<endl;
     for (int i = 0; i < event.size(); ++i){
	if (abs(event[i].status())<=30){	
		std::cout<<i<<"\t"<<event[i].name()<<"\t"<<event[i].status()<<"\t"<<event[i].mother1()<<"   "<<event[i].mother2()<<"\t\t"<<event[i].daughter1()<<"   "<<event[i].daughter2()<<"\t\t"<<event[i].pT()<<"\t"<<event[i].eta()<<"\t"<<event[i].phi()<<endl;
	}
	

    }
    
    std::cout<<"=========================event-END==============================="<<endl;
    myout<<"=========================event-END==============================="<<endl;

    //to find g for matching jets
    std::cout<<"================================status-23==================================="<<endl;
    myout<<"================================status-23==================================="<<endl;
    std::cout<<"no.\teventid\tstatus\tmothers\t\tdaughters\tpT\teta\t\phi\n"<<endl;
    myout<<"no.\teventid\tstatus\tmothers\t\tdaughters\tpT\teta\t\phi\n"<<endl;
    for (int i = 0; i < event.size(); ++i){
	//qq_index=qq.find(event[i].id());
 	//std::cout<<"=====================qq_index==========================="<<*qq_index<<endl;
	
	if (event[i].status()==-23){
		
		std::cout<<i<<"\t"<<event[i].name()<<"\t"<<event[i].status()<<"\t"<<event[i].mother1()<<"   "<<event[i].mother2()<<"\t\t"<<event[i].daughter1()<<"   "<<event[i].daughter2()<<"\t\t"<<event[i].pT()<<"\t"<<event[i].eta()<<"\t"<<event[i].phi()<<endl;

		 myout<<i<<"\t"<<event[i].name()<<"\t"<<event[i].status()<<"\t"<<event[i].mother1()<<"   "<<event[i].mother2()<<"\t\t"<<event[i].daughter1()<<"   "<<event[i].daughter2()<<"\t\t"<<event[i].pT()<<"\t"<<event[i].eta()<<"\t"<<event[i].phi()<<endl;
		if(qq.find(event[i].id())!=qq.end()){    //for qq (qq.find(event[i].id())!=qq.end() ) , for gg event[i].id()==21
			std::cout<<"=====================ture==========================="<<endl;
			myout<<"=====================ture==========================="<<endl;
			eta0.insert(eta0.begin(), event[i].eta());
			phi0.insert(phi0.begin(), event[i].phi());
		}
	
	}
	

    }
    std::cout<<"==============================status-23-END=============================="<<endl;
    myout<<"==============================status-23-END=============================="<<endl;
    std::cout<<"===============================================================================phi0size==============================================================================================="<<
			phi0.size()<<endl;
    myout<<"===============================================================================phi0size==============================================================================================="<<
			phi0.size()<<endl;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {

      // Require visible/charged particles inside detector.
      if      (select > 2 &&  event[i].isNeutral() ) continue;
      else if (select == 2 && !event[i].isVisible() ) continue;
      if (etaMax < 20. && abs(event[i].eta()) > etaMax) continue;

      // Create a PseudoJet from the complete Pythia particle.
      fastjet::PseudoJet particleTemp = event[i];

      // Optionally modify mass and energy.
      pTemp = event[i].p();
      
      mTemp = event[i].m();
     
      if (massSet < 2) {
        
        mTemp = (massSet == 0 || event[i].id() == 22) ? 0. : 0.13957;
        pTemp.e( sqrt(pTemp.pAbs2() + mTemp*mTemp) );
        particleTemp.reset_momentum( pTemp.px(), pTemp.py(),
           pTemp.pz(), pTemp.e() );
      }

      // Store acceptable particles as input to Fastjet.
      // Conversion to PseudoJet is performed automatically
      // with the help of the code in FastJet3.h.
      fjInputs.push_back( particleTemp);

      ++nAnalyze;
    }

    // Run Fastjet algorithm and sort jets in pT order.
    vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    inclusiveJets = clustSeq.inclusive_jets(pTMin);
    sortedJets    = sorted_by_pt(inclusiveJets);
    clock_t aftFast = clock();
    
    // List first few FastJet jets and some info about them.
    // Note: the final few columns are illustrative of what information
    // can be extracted, but does not exhaust the possibilities.
    if (iEvent < nListJets) {
      cout << "\n --------  FastJet jets, p = " << setw(2) << power
           << "  --------------------------------------------------\n\n "
           << "  i         pT        y      phi  mult chgmult photons"
           << "      hardest  pT in neutral " << endl
           << "                                                       "
           << "  constituent        hadrons " << endl;
   
      for (int i = 0; i < int(sortedJets.size()); ++i) {
        
        vector<fastjet::PseudoJet> constituents
          = sortedJets[i].constituents();
	//print jet kenetic
        std::cout << "===========================================================================\nsortedJets\tpt\teta\tphi\tconstituents_size=\t"<<sortedJets[i].pt()<<"\t"<<sortedJets[i].eta()<<"\t"<<sortedJets[i].phi()
			<<"\t"<<constituents.size()<<endl;
        myout<< "===========================================================================\nsortedJets\tpt\teta\tphi\tconstituents_size=\t"<<sortedJets[i].pt()<<"\t"<<sortedJets[i].eta()<<"\t"<<sortedJets[i].phi()
			<<"\t"<<constituents.size()<<endl;
	//match jet with gluon/quark eta and phi.
	for (int j = 0;j < int(phi0.size()); ++j)
	if(deltaR(phi0[j],sortedJets[i].phi(),eta0[j],sortedJets[i].eta())<d_R){
		std::cout<<"gluonjet"<<endl;  //print ''gluon jet''
		myfile<<"\n"<<sortedJets[i].e()<<"\t"<<sortedJets[i].pt()<<"\t"<<sortedJets[i].eta()<<"\t"<<sortedJets[i].phi()<<"\t"<<constituents.size()<<"\n"<<endl; //save jet pt eta phi in myfile
		myout<<"gluonjet"<<endl;  //print ''gluon jet''
		for (int j=0;j<int(constituents.size()); ++j){
			
			std::cout<<"\nJet_constituent\tpt\teta\tphi\tpID=\t"<<constituents[j].pt()<<"\t"<<constituents[j].eta()<<
			"\t"<<constituents[j].phi()<<"\t"<<constituents[j].user_info<Particle>().name()<<endl;   
			myfile<<constituents[j].e()<<"\t"<<constituents[j].pt()<<"\t"<<constituents[j].eta()<<"\t"<<constituents[j].phi()<<endl; //save jet constituents pt eta phi in myfile

			myout<<"\nJet_constituent\tpt\teta\tphi\tpID=\t"<<constituents[j].pt()<<"\t"<<constituents[j].eta()<<
			"\t"<<constituents[j].phi()<<"\t"<<constituents[j].user_info<Particle>().name()<<endl;   
			
			}
	}
        fastjet::PseudoJet hardest
          = fastjet::SelectorNHardest(1)(constituents)[0];
        vector<fastjet::PseudoJet> neutral_hadrons
          = ( fastjet::SelectorIsHadron()
           && fastjet::SelectorIsNeutral())(constituents);
        double neutral_hadrons_pt = join(neutral_hadrons).perp();
	
//        cout << setw(4) << i << fixed << setprecision(3) << setw(11)
//             << sortedJets[i].perp() << setw(9)  << sortedJets[i].rap()
//             << setw(9) << sortedJets[i].phi_std()
//             << setw(6) << constituents.size()
//             << setw(8) << fastjet::SelectorIsCharged().count(constituents)
//             << setw(8) << fastjet::SelectorId(22).count(constituents)
//             << setw(13) << hardest.user_info<Particle>().name()
//             << "     " << setw(10) << neutral_hadrons_pt << endl;
      }
      myfile << "\t";
      cout << "\n --------  End FastJet Listing  ------------------"
           << "---------------------------------" << endl;
    }
  
   

 
  // Done.
} 
 //=================================================================
   //close file
  myfile.close();
  myout.close();
  return 0;
}
