// -*- C++ -*-
//#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Analysis.hh"
//#include "Rivet/Cuts.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...


namespace Rivet {


  //using namespace Cuts;
  

  class CMS_2014_I9999999 : public Analysis {
  public:

    /// Constructor
    CMS_2014_I9999999()
      : Analysis("CMS_2014_I9999999")
    {    }


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// @todo Initialise and register projections here

      /// @todo Book histograms here, e.g.:
      // _h_XXXX = bookProfile1D(1, 1, 1);
      // _h_YYYY = bookHisto1D(2, 1, 1);

      //Cut cuts = etaIn(-10,10) & (pT >= 0.0*GeV);
      //WFinder wfinder_bare_el(fs, cuts, PID::ELECTRON, 0*GeV, 1000*GeV, 25*GeV, 0.0, WFinder::NOCLUSTER);
      //addProjection(wfinder_el, "WFinder_el");

      IdentifiedFinalState ifs1(-2.1, 2.1, 6.0*GeV);
      ifs1.acceptIdPair(PID::TAU);
      //ifs1.acceptIdPair(PID::ELECTRON);
      addProjection(ifs1, "IFS1");

      const FastJets jets(FinalState(-10, 10, 0.0*GeV), FastJets::ANTIKT, 0.5);

      addProjection(jets, "Jets");

      _h_fiducial = bookHisto1D(1,1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      //const double weight = event.weight();

      //Particles taus = applyProjection<IdentifiedFinalState>(event, "IFS1").particlesByPt();

      //if (taus.size() > 0)
      //	std::cout << "taus.size() = " << taus.size() << std::endl;

      vector<const GenParticle*> quarks;
      vector<const GenParticle*> leptons;

      std::cout << "andrew debug -1" << std::endl;

      foreach (const GenParticle* p, particles(event.genEvent())) {
	//n_particles++;

	std::cout << "andrew debug -2" << std::endl;

	int n_2212_ancestors = 0;
	int n_quark_ancestors = 0;

	const vector<GenParticle*> ancestors = particles_in(p, HepMC::ancestors);

	if (p->status() == 3){

	  std::cout << "andrew debug -3" << std::endl;

	  if (ancestors.size() == 1){
	    
	    foreach (const GenParticle* pa, ancestors) {
	      
	      if(pa->pdg_id() == 2212)
		n_2212_ancestors++;
	    }
	    
	    assert(n_2212_ancestors == 1);
	    
	    continue;
	    
	  }

	  std::cout << "andrew debug -3.1" << std::endl;

	  if (ancestors.size() == 2){
	    
	    foreach (const GenParticle* pa, ancestors) {
	      
	      if(pa->pdg_id() == 2212)
		n_2212_ancestors++;
	      
	      if(pa->pdg_id() <= 6 && pa->pdg_id() >= -6 && pa->pdg_id() != 0)
		n_quark_ancestors++;
	      
	    }
	    
	    assert(n_2212_ancestors == 1 && n_quark_ancestors == 1);
	    
	    continue;
	    
	  }

	  std::cout << "andrew debug -3.2" << std::endl;

	  if (abs(p->pdg_id()) == 1 || abs(p->pdg_id()) == 2 || abs(p->pdg_id()) == 3 || abs(p->pdg_id()) == 4 || abs(p->pdg_id()) == 5 || abs(p->pdg_id()) == 6){
	    quarks.push_back(p); 
	  }
	  
	  if (abs(p->pdg_id()) == 11 || abs(p->pdg_id()) == 13 || abs(p->pdg_id()) == 15){
	    leptons.push_back(p); 
	  }
	  
	}
      }

      std::cout << "andrew debug -3.5" << std::endl;

      assert(leptons.size() == 2);
      assert(quarks.size() == 2);

      if (leptons[0]->momentum().perp() < 10)
	vetoEvent;
      
      std::cout << "andrew debug -4" << std::endl;

      if (leptons[1]->momentum().perp() < 10)
	vetoEvent;

      std::cout << "andrew debug -4.1" << std::endl;

      if (quarks[0]->momentum().perp() < 20)
	vetoEvent;

      std::cout << "andrew debug -4.2" << std::endl;

      if (quarks[1]->momentum().perp() < 20)
	vetoEvent;

      std::cout << "andrew debug -4.3" << std::endl;

      if (fabs(quarks[0]->momentum().eta()) > 5.0)
	vetoEvent;
      
      std::cout << "andrew debug -4.4" << std::endl;

      if (fabs(quarks[1]->momentum().eta()) > 5.0)
	vetoEvent;

            std::cout << "andrew debug -4.5" << std::endl;

      if (fabs(leptons[0]->momentum().eta()) > 2.5)
	vetoEvent;

            std::cout << "andrew debug -4.6" << std::endl;

      if (fabs(leptons[1]->momentum().eta()) > 2.5)
	vetoEvent;

      std::cout << "fabs(leptons[1]->momentum().eta()) = " << fabs(leptons[1]->momentum().eta()) << std::endl;

      std::cout << "andrew debug -5.1" << std::endl;


      Rivet::FourMomentum quark0_four_momentum(quarks[0]->momentum());
      Rivet::FourMomentum quark1_four_momentum(quarks[1]->momentum());

      if ((quark0_four_momentum + quark1_four_momentum).mass() < 300)
	vetoEvent;

      std::cout << "(quark0_four_momentum + quark0_four_momentum).mass() = " << (quark0_four_momentum + quark0_four_momentum).mass() << std::endl;

      std::cout << "andrew debug -5.2" << std::endl;

      if ( fabs(quarks[0]->momentum().eta() - quarks[1]->momentum().eta()) < 2.5 )
	vetoEvent;

      std::cout << "fabs(quarks[0]->momentum().eta() - quarks[1]->momentum().eta()) = " << fabs(quarks[0]->momentum().eta() - quarks[1]->momentum().eta()) << std::endl;

      std::cout << "andrew debug -5.3" << std::endl;

      std::cout << "event.weight() = " << event.weight() << std::endl;

      _h_fiducial->fill(7000, event.weight());

      /// @todo Do the event by event analysis here

    }    


    /// Normalise histograms etc., after the run
    void finalize() {

      /// @todo Normalise, scale and otherwise manipulate histograms here

      // scale(_h_YYYY, crossSection()/sumOfWeights()); // norm to cross section
      // normalize(_h_YYYY); // normalize to unity

    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_fiducial;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2014_I9999999);

}
