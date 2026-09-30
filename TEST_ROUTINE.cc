//Use Rivet as analyze-framework for already MC-generated (e.g. Pythia, Herwig) events.
//One event from MC-Simulation (based on theory) contains particle-list, 4-moments, PID...
//PID: particle id, one for every particle from final state, + (particle) and - (anti particle) 
//Classic Rivet analysis separated in init(), analyze(), finalize().
//init(): define projections (particle level with jets), define histograms, set cuts
//analyze(): use projections, calculate observables, fill histograms (loop over every event)
//finalize(): normalize histograms
//particle level: stable particles; final state: sum of all stable particles per event from MC,...
//...Rivet reconstructs Jets, Leptons, etc. from final state particle list 


// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/VisibleFinalState.hh"




//ATTENTION!!!
//SOME OF THE COMMENTS MADE IN THE FOLLOWING
//MIGHT NOT BE QUITE RIGHT!!!
//isabel.molinos-pineiro@tu-dortmund.de




namespace Rivet {

  //@brief ttbar + gamma at 13 TeV; TEST_ROUTINE name of analysis.
  class TEST_ROUTINE : public Analysis {
  public:

	  
    //standard constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TEST_ROUTINE);


    //Book histograms and initialise projections before the run.
    void init() {

      //Example for a global variable, here: counter.
      eventCounter = 0;
      
      //First cut to veto leptons (extremly soft).
      Cut eta_full = Cuts::abseta < 5.0 && Cuts::pT > 1.0*MeV;

      // All final state particles, necessary to get neutrinos.
      FinalState fs(eta_full);
      
      
      //Get all photons of final state to dress leptons. Create projection-object with type
      //FinalState and save as variable "photons". Filter condition inside brackets: pid contains 
      //identif. number of particle/photon. Take all particles, surviving selection an save in list
      //(done per event)
      FinalState photons(Cuts::abspid == PID::PHOTON);

      
      //Projection to find electrons (in two steps, finding and then reconstructing).
      //prompt_elec finds all elecs not from hadron decay in MC (e.g. W, Z, tau; b/c-elecs removed)...
      //...(MC-Simulations record decay chains) --> leptonic spin analyzer
      PromptFinalState prompt_el(Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);

      //LeptonFinder reconstructs obj. (add photons on cone: dressed leptons) + use analysis-cuts...
      //...from PromptFinalState above  
      LeptonFinder elecs(prompt_el, photons, 0.1, (Cuts::abseta < 2.5 && Cuts::pT > 27*GeV));

      //Before jets are clustered, all prompt elecs are removed. First LF has hard cuts for...
      //...further analysis, second LF contains all soft ones to cut them completly from jet. 
      LeptonFinder veto_elecs(prompt_el, photons, 0.1, eta_full);
      declare(elecs, "elecs");
      

      //Projection to find the muons. See explanations above.
      PromptFinalState prompt_mu(Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT);
      LeptonFinder muons(prompt_mu, photons, 0.1, (Cuts::abseta < 2.5 && Cuts::pT > 27*GeV));
      LeptonFinder veto_muons(prompt_mu, photons, 0.1, eta_full);
      declare(muons, "muons");
      
      
      //Projection for PROMPT neutrinos (invisible), should not be part of jet.
      //From jet input, veto leptons are removed with implemented lepton projections
      const InvisibleFinalState invis(OnlyPrompt::YES, TauDecaysAs::PROMPT);
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(veto_elecs);
      vfs.addVetoOnThisFinalState(veto_muons);
      vfs.addVetoOnThisFinalState(invis);

      //smallR jet projection: every final state particle with deltaR less x clustered later...
      //...on inside one jet (minus prompt leptons with very soft cuts). It is important...
      //...to remove prompt leptons because statist. they can be inside jet-cone and be...
      //...and may be incorrectly counted twice, which distorts jet-kinematics.
      FastJets smallRjets(vfs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::ALL);
      declare(smallRjets, "smallRjets");
 
      
      //Projection of largeR/reclustered (!) jets. Construction of boosted phase space.
      FastJets largeRjets(vfs, JetAlg::ANTIKT, 1.0);
      declare(largeRjets, "largeRjets");

      //Get all neutrinos/invisbles, not only prompt for MET-construction (do not use...
      //...neutrino-projection above, since only prompt ones are used there).
      VetoedFinalState inv_fs(fs);
      inv_fs.addVetoOnThisFinalState(VisibleFinalState(fs));
      declare(inv_fs, "InvisibleFS");


      //Here, define the histograms you want to make, first number is number of bins...
      //...and other two numbers are upper and lower range of plot.
      book(_h["smallRjets_n"],  "smallRjets_n",  11, 26.5, 10.5);
      book(_h["smallRjets_pT_all"], "smallRjets_pT_all", 30, 26.5, 389.5);
      book(_h["smallRjets_pT_lead"], "smallRjets_pT_lead", 30, 26.5, 519.5);
      book(_h["lep_pT"], "lep_pT", 30, -0.5, 430.5);
      book(_h["elec_pT"], "elec_pT", 30, -0.5, 399.5);
      book(_h["muon_pT"], "muon_pT", 30, -0.5, 399.5);
      book(_h["bT_smallRjets_n_precut"], "bT_smallRjets_n_precut", 7, -0.5, 6.5);
      book(_h["bT_smallRjets_n"], "bT_smallRjets_n", 6, -0.5, 5.5);
      book(_h["bT_smallRjets_pT_all"], "bT_smallRjets_pT_all", 30, -0.5, 299.5);
      book(_h["largeRjets_n"], "largeRjets_n", 7, -0.5, 6.5); 
      book(_h["largeRjets_mass_all"], "largeRjets_mass_all", 50, -0.5, 269.5);
      book(_h["hadrTop_smallRjets_n"], "hadrTop_smallRjets_n", 4, 0.5, 4.5);
      book(_h["hadrTopJet_mass"], "hadrTopJet_mass", 25, 99.5, 239.5);
      book(_h["hadrTopJet_pT"], "hadrTopJet_pT", 30, 219.5, 829.5);
      book(_h["sumLep_bTjet_mass"], "sumLep_bTjet_mass", 30, -0.5, 199.5);
      book(_h["sumLep_bTjet_pT"], "sumLep_bTjet_pT", 30, -0.5, 699.5);
      book(_h["lepTop_mass"], "lepTop_mass", 30, 79.5, 339.5);
      book(_h["lepTop_pT"], "lepTop_pT", 30, -0.5, 889.5);    
      book(_h["ttbar_mass"], "ttbar_mass", 40, 699.5, 1849.5);
      book(_h["ttbar_pT"], "ttbar_pT", 30, -0.5, 740.5);
      book(_h["ttbar_CMttbar_pT"], "ttbar_CMttbar_pT", 10, 0.0, 5.0);
      book(_h["ttbar_CMttbar_pz"], "ttbar_CMttbar_pz", 10, -5.0, 5.0);
      book(_h["ttbar_CMttbar_mass"], "ttbar_CMttbar_mass", 45, 500, 1400.0);

      book(_h["W_hadrAnalyzer_top"], "W_hadrAnalyzer_top", 4, 0, 500);
      book(_h["W_hadrAnalyzer_antitop"], "W_hadrAnalyzer_antitop", 4, 0, 500);
      
      //Polarisations.
      book(_h["cosTheta_k"], "cosTheta_k", 4, -1.0, 1.0);
      book(_h["cosThetaBar_k"], "cosThetaBar_k", 4, -1.0, 1.0);
      book(_h["cosTheta_r"], "cosTheta_r", 4, -1.0, 1.0);
      book(_h["cosThetaBar_r"], "cosThetaBar_r", 4, -1.0, 1.0);
      book(_h["cosTheta_n"], "cosTheta_n", 4, -1.0, 1.0);
      book(_h["cosThetaBar_n"], "cosThetaBar_n", 4, -1.0, 1.0);

      //Sanity plots for spin analyzer.
      book(_h["cosTheta_k_lep"], "cosTheta_k_lep", 4, -1.0, 1.0);
      book(_h["cosTheta_r_lep"], "cosTheta_r_lep", 4, -1.0, 1.0);
      book(_h["cosTheta_n_lep"], "cosTheta_n_lep", 4, -1.0, 1.0);
      book(_h["cosThetaBar_k_hadr"], "cosThetaBar_k_hadr", 4, -1.0, 1.0);
      book(_h["cosThetaBar_r_hadr"], "cosThetaBar_r_hadr", 4, -1.0, 1.0);
      book(_h["cosThetaBar_n_hadr"], "cosThetaBar_n_hadr", 4, -1.0, 1.0);

      book(_h["cosTheta_k_hadr"], "cosTheta_k_hadr", 4, -1.0, 1.0);
      book(_h["cosTheta_r_hadr"], "cosTheta_r_hadr", 4, -1.0, 1.0);
      book(_h["cosTheta_n_hadr"], "cosTheta_n_hadr", 4, -1.0, 1.0);
      book(_h["cosThetaBar_k_lep"], "cosThetaBar_k_lep", 4, -1.0, 1.0);
      book(_h["cosThetaBar_r_lep"], "cosThetaBar_r_lep", 4, -1.0, 1.0);
      book(_h["cosThetaBar_n_lep"], "cosThetaBar_n_lep", 4, -1.0, 1.0);


      //Correlations.
      book(_h["C_kk"], "C_kk", {-1.0, -0.3, 0, 0.3, 1.0});
      book(_h["C_rr"], "C_rr", {-1.0, -0.5, -0.2, 0, 0.2, 0.5, 1.0});
      book(_h["C_nn"], "C_nn", {-1.0, -0.4, 0, 0.4, 1.0});

      book(_h["C_nk"], "C_nk", 4, -1.0, 1);
      book(_h["C_kn"], "C_kn", 4, -1.0, 1);
      book(_h["C_nr"], "C_nr", 4, -1.0, 1);
      book(_h["C_rn"], "C_rn", 4, -1.0, 1);
      book(_h["C_rk"], "C_rk", 4, -1.0, 1);
      book(_h["C_kr"], "C_kr", 4, -1.0, 1);


    }


    void analyze(const Event& event) { 
     

      eventCounter++;

      if (eventCounter % 1000 == 0) {
        std::cout << "Processing event: " << eventCounter << std::endl;
      }
	    

      //smallR jet-projection and analysis-cuts. Jets of each event sorted by pT.
      DressedLeptons elecs = apply<LeptonFinder>(event, "elecs").dressedLeptons();
      DressedLeptons muons = apply<LeptonFinder>(event, "muons").dressedLeptons();
      Jets smallRjets = apply<FastJets>(event, "smallRjets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      
      //Projection for MET (all invisible particles, not only prompt).
      const FinalState& ifs = apply<FinalState>(event, "InvisibleFS");

      //Calculate MET. Actually, this gives px, py and pz on truth-level; detectors are...
      //...not able to meassure pz due to beam geometry. Get pz from W-Massconstrain.
      FourMomentum met;
      for (const Particle& p : ifs.particles())  met += p.momentum();
      

      //largeR jet-projection with recommended analysis-cuts. Jets of each event sorted by pT.
      Jets largeRjets = apply<FastJets>(event, "largeRjets").jetsByPt(Cuts::pT > 150*GeV && Cuts::abseta < 2.0);
      

      //PRE-SELECTION of smallRjets (non-boosted)
      //JET-SELECTION: at least 3 smallR jets AND at least 2 bTagged smallR jets per event
      if (smallRjets.size() < 3) vetoEvent; 

      //Filling some sanity jet-histograms.
      _h["smallRjets_n"] -> fill(smallRjets.size());
      _h["smallRjets_pT_lead"] -> fill(smallRjets[0].pT()/GeV);
      

      //Mechanism: "For jet in container smallRjets fill...".
      for (const Jet& jet : smallRjets) {
        _h["smallRjets_pT_all"] -> fill(jet.pT()/GeV);
      }
      
      
      //b-tagging: loop over jet list of each event and analyse, if b-hadrons are contained (hint...
      //...that jet originates from b-quark). Two are requested for event selection.
      Jets bT_smallRjets;

      for (const Jet& jet : smallRjets) {
        if (jet.bTagged()) bT_smallRjets += jet;
      }

      _h["bT_smallRjets_n_precut"] -> fill(bT_smallRjets.size());

      if (bT_smallRjets.size() < 2) vetoEvent;

      
      //Fill bTag-histograms after cut.
      _h["bT_smallRjets_n"] -> fill(bT_smallRjets.size());
      
      for (const Jet& jet : bT_smallRjets) {
        _h["bT_smallRjets_pT_all"] -> fill(jet.pT()/GeV);
      }


      //BOOSTED SELECTION: hadronic decay
      //First cut: at least one largeRjet with specif. pT and mass in range. Loop over largeRjets list...
      //...and find candidate for hadronic top decay (hadrTopCandid collects all rcJet selected by cuts)
      Jets hadrTopCandid;
 
      for (const Jet& jet : largeRjets) {
        if (jet.pT() > 300*GeV && jet.mass() > 120*GeV && jet.mass() < 220*GeV) {
          hadrTopCandid += jet;
        }
      }

      if (hadrTopCandid.empty()) vetoEvent;


      //largeRjet-histograms: multiplicity per event and mass of all largeRjets (after first cuts).
      _h["largeRjets_n"] -> fill(largeRjets.size());

      for (const Jet& jet : largeRjets) {
        _h["largeRjets_mass_all"] -> fill(jet.mass()/GeV);
      }
       
     
      //Second cut: select the rcJet, with mass nearest to top-mass an call it hadrTopJet.
      //Declare single jet object and loop over every hadrTopCandid to find the rcJet with lowest mass difference...
      //...to initialise it with hadrTopJet (plot some properties of this jet).
      Jet hadrTopJet;

      //Select large starting value for mass difference and loop. If jet inside lR_hadrTopCandid has smaller difference,...
      //...overwrite hadrTopJet.
      double minMassDiff = 1e9;
      const double mTop = 172.5*GeV;

      for (const Jet& jet : hadrTopCandid) {
        double massDiff = fabs(jet.mass() - mTop);
        
	if (massDiff < minMassDiff) {
          minMassDiff = massDiff;
	  hadrTopJet = jet;
	}	
      }



      // Check overlap of jets/leptons (jets already reconstructed). Prompt leptons already...              
      //...removed from jet, but NOW, also remove prompt leptons that are to near to jet.
      for (const Jet& jet : smallRjets) {
	idiscard(elecs, deltaRLess(jet, 0.4));
	idiscard(muons, deltaRLess(jet, 0.4));
      }

     
      //LEPTON SELECTION: If the event does not contain a lepton, do not fill the histograms.
      Particle lepPart;
      FourMomentum lep;

      if (elecs.size() == 1) {
        lepPart = elecs[0];
	lep = lepPart.momentum();
	_h["elec_pT"] -> fill(lep.pT()/GeV);
      }
      else if (muons.size() == 1) {
        lepPart = muons[0];
	lep = lepPart.momentum();
	_h["muon_pT"] -> fill(lep.pT()/GeV);
      }
      else {
        vetoEvent;
      } 

      _h["lep_pT"] -> fill(lep.pT()/GeV);

      
      //Continue BOOSTED SELECTION: third Cut, hadronic decay.
      if (deltaR(hadrTopJet, lep) <= 1.0) vetoEvent;
      

      //Fourth Cut: Check if hadrTopJet contains 2-3 smallR jets (called hadrTop_smallRjets).
      //Loop over all smallR jets and see if they lay around hadrTopJet.    
      Jets hadrTop_smallRjets;

      for (const Jet& jet : smallRjets) {
        if (deltaR(jet, hadrTopJet) < 1.0) {
          hadrTop_smallRjets += jet;  
	}
      }

      if (hadrTop_smallRjets.size() < 2 || hadrTop_smallRjets.size() > 3) vetoEvent;

      _h["hadrTop_smallRjets_n"] -> fill(hadrTop_smallRjets.size());

     

      //Fifth Cut: Check if one constituent inside hadrTopJet is bTagged.
      int nr_bT_smallRjets = 0;

      for (const Jet& jet : hadrTop_smallRjets) {
        if (jet.bTagged()) {
          nr_bT_smallRjets++;
        }
      }
    
      if (nr_bT_smallRjets != 1) vetoEvent;


      FourMomentum hadrTop = hadrTopJet.momentum();

      _h["hadrTopJet_mass"] -> fill(hadrTopJet.mass()/GeV);
      _h["hadrTopJet_pT"] -> fill(hadrTopJet.pT()/GeV);



      //BOOSTED SELECTION: leptonic decay.
      //Sixth cut. Finding bT_smallRjet on leptonic side.
      Jet lept_bT_smallRjet;

      double minDeltaR = 1e9;

      for (const Jet& jet : bT_smallRjets) {
         
	 double dR = deltaR(jet, lep);

	 if (dR <= minDeltaR) {
           minDeltaR = dR;
           lept_bT_smallRjet = jet;
	}
      }
    
      if (minDeltaR > 2) vetoEvent;
    


      //Last Cut: See if invariant mass of lept_bT_smallRjet and lep less than 180 GeV.
      //lep already defined als 4-moment, jet has more information. Use jet.momentum()...
      //...to fulfill proper vector addition. ".mass()" calculates inv. mass as seen in theory.
      FourMomentum sumLep_bTjet = lep + lept_bT_smallRjet.momentum();
      
      if (sumLep_bTjet.mass() > 180*GeV) vetoEvent;

      _h["sumLep_bTjet_mass"] -> fill(sumLep_bTjet.mass()/GeV);
      _h["sumLep_bTjet_pT"] -> fill(sumLep_bTjet.pT()/GeV);

    
    //MET-cuts (new cuts not given in ATLAS draft).
    double mTW = sqrt(2 * lep.pT() * met.pT() * (1 - cos(deltaPhi(lep, met))));

    if (met.pT() < 40*GeV) vetoEvent;
    if (mTW < 60*GeV) vetoEvent;

    //Neutrino-reconstruction.
    double pz_nu = _computeneutrinoz(lep, met);
    double E_nu = sqrt(sqr(met.px()) + sqr(met.py()) + sqr(pz_nu));

    FourMomentum neutrino(E_nu, met.px(), met.py(), pz_nu);


    //Full 4-moment reconstruction of leptonic top.
    
    FourMomentum lepTop = lep + neutrino + lept_bT_smallRjet.momentum();


    _h["lepTop_mass"] -> fill(lepTop.mass()/GeV);
    _h["lepTop_pT"] -> fill(lepTop.pT()/GeV);


    //ttbar 4-momentum!
    //UND NEUER CUT AUF MTTBAR (ATLAS Draft).
    FourMomentum ttbar = hadrTop + lepTop;
    if (ttbar.mass() < 800*GeV) vetoEvent;

    _h["ttbar_mass"] -> fill(ttbar.mass()/GeV);   
    _h["ttbar_pT"] -> fill(ttbar.pT()/GeV);

    //---------------------------------------------------------------------
    //HELICITY BASIS-------------------------------------------------------
    //---------------------------------------------------------------------

    //Start with k (direction of top and antitop in CM). First move lepTop and hadrTop into CM-ttbar frame.
    LorentzTransform boost_CM_ttbar = LorentzTransform::mkFrameTransformFromBeta(ttbar.betaVec());

    //To get angular distributions of t and tbar, differentiation between top and antitop is necessary...
    //...to reconstruct helicity basis correct, since both can decay hadronically and leptonically...
    //...(decide by charge of lepton)
    FourMomentum top;
    FourMomentum antitop;

    if (lepPart.charge() > 0) {
      top = lepTop;
      antitop = hadrTop;
    }
    else {
      top = hadrTop;
      antitop = lepTop;
    }

    
    //Boosting top and antitop in CMttbar.
    FourMomentum top_CM_ttbar = boost_CM_ttbar.transform(top);
    FourMomentum antitop_CM_ttbar = boost_CM_ttbar.transform(antitop);


    //Construct k and kbar along direction of each top.
    //Everything implemented as in 1508.05271v2!
    Vector3 k = top_CM_ttbar.p3().unit();
    Vector3 kbar = -k;


    //Now, r and rbar:
    Vector3 p(0.0, 0.0, 1.0); //Beam axis.
    double y_p = p.dot(k); //cos(angle): p and k.
    double r_p = sqrt(1.0 - y_p * y_p); //Normalisation.  

    double sign_yp;

    if (y_p >= 0.0) {
      sign_yp = 1.0;
    }
    else {
      sign_yp = -1.0;
    }
    
    
    Vector3 r = sign_yp * (p - y_p * k) / r_p;
    Vector3 rbar = -r;


    //Lastly, n and nbar:
    Vector3 n = sign_yp * (p.cross(k)) / r_p;
    Vector3 nbar = -n;


    //Now, get boost of top and antitop into their rest frame since angles meassured here. k and kbar do NOT...
    //...need to be boosted, since the boost into rest frame is along k or kbar, and therefore, not necessary.
    LorentzTransform boost_rF_top = LorentzTransform::mkFrameTransformFromBeta(top_CM_ttbar.betaVec());
    LorentzTransform boost_rF_antitop = LorentzTransform::mkFrameTransformFromBeta(antitop_CM_ttbar.betaVec());


    //Sanity Output. Extremly small values or 0 for dot-products and amounts equal to 1.
    //std::cout << "k·r = " << k.dot(r)
    //          << "  k·n = " << k.dot(n)
    //          << "  r·n = " << r.dot(n) << std::endl;

    //std::cout << "|k| = " << k.mod()
    //           << "  |r| = " << r.mod() << "  |n| = " << n.mod() <<
    //           std::endl;
    


    //Theta-Cut (ATLAS Draft):
    if (fabs(y_p) >= 0.4) vetoEvent;
    
    

    //Construction of the hadronic spin analyzer: W_bos (4-Moment) and q_optHadrPol (vector).
    FourMomentum hadrAnalyzer_rF;
    Vector3 hadrAnalyzer;

    if (hadrTop_smallRjets.size() == 2) {

      //Get W-boson as spin analyzer. hadrTop_smallRjets contain either 2 or 3 smallRjets (existing in largeRjet).
      for (const Jet& jet : hadrTop_smallRjets) {
        if (jet.bTagged() == false) {
          FourMomentum W = jet.momentum();

	  //Sanity Plot to see, how many W bosons are taken into account for hadr spin analyzer (in .yoda available).


          if (lepPart.charge() > 0) {
	    FourMomentum hadrAnalyzer_CM_ttbar = boost_CM_ttbar.transform(W);
	    hadrAnalyzer_rF = boost_rF_antitop.transform(hadrAnalyzer_CM_ttbar);
	    hadrAnalyzer = hadrAnalyzer_rF.p3().unit();
	     _h["W_hadrAnalyzer_antitop"]->fill(hadrAnalyzer_rF.mass());

	  } else {
            FourMomentum hadrAnalyzer_CM_ttbar = boost_CM_ttbar.transform(W);
	    hadrAnalyzer_rF = boost_rF_top.transform(hadrAnalyzer_CM_ttbar);
	    hadrAnalyzer = hadrAnalyzer_rF.p3().unit();
	     _h["W_hadrAnalyzer_top"]->fill(hadrAnalyzer_rF.mass());
          }
	  
        }
      }

    } else {

      //For 3 constituents, use the hadronic polarimeter.
      Jet bjet, wjet1, wjet2;
      Jets wjets;

      for (const Jet& jet : hadrTop_smallRjets) {
        if (jet.bTagged() == true) {
          bjet = jet;
        } else {
          wjets += jet;
	}
      }

      wjet1 = wjets[0];
      wjet2 = wjets[1];

      //Boost all constituent jets into top/antitop RF to see which wjet soft and which one hard.
      FourMomentum bjet_rF, wjet1_rF, wjet2_rF; 
      FourMomentum q_hard_rF, q_soft_rF;
      
     
      if (lepPart.charge() > 0) {
	FourMomentum bjet_CM_ttbar = boost_CM_ttbar.transform(bjet.momentum());
	FourMomentum wjet1_CM_ttbar = boost_CM_ttbar.transform(wjet1.momentum());
	FourMomentum wjet2_CM_ttbar = boost_CM_ttbar.transform(wjet2.momentum());

        bjet_rF = boost_rF_antitop.transform(bjet_CM_ttbar);
        wjet1_rF = boost_rF_antitop.transform(wjet1_CM_ttbar);
        wjet2_rF = boost_rF_antitop.transform(wjet2_CM_ttbar);  
        
	if (wjet1_rF.E() > wjet2_rF.E()) {
          q_hard_rF = wjet1_rF;
	  q_soft_rF = wjet2_rF;
        } else {
          q_soft_rF = wjet1_rF;
          q_hard_rF = wjet2_rF;	  
	}
        	
      } else {
        FourMomentum bjet_CM_ttbar = boost_CM_ttbar.transform(bjet.momentum());
        FourMomentum wjet1_CM_ttbar = boost_CM_ttbar.transform(wjet1.momentum());
        FourMomentum wjet2_CM_ttbar = boost_CM_ttbar.transform(wjet2.momentum());

        bjet_rF = boost_rF_top.transform(bjet_CM_ttbar);
        wjet1_rF = boost_rF_top.transform(wjet1_CM_ttbar);
        wjet2_rF = boost_rF_top.transform(wjet2_CM_ttbar);

        if (wjet1_rF.E() > wjet2_rF.E()) { 
          q_hard_rF = wjet1_rF;
          q_soft_rF = wjet2_rF;
        } else {
          q_soft_rF = wjet1_rF;
          q_hard_rF = wjet2_rF;
        }
      }

      //Reconstruct W boson in t or tbar RF and boost bjet and soft jet (d-type, 61%) into W RF.
      //Then calculate the helicity angle.
      FourMomentum W_rF = q_hard_rF + q_soft_rF;  
      LorentzTransform boost_rF_W = LorentzTransform::mkFrameTransformFromBeta(W_rF.betaVec());
      
      FourMomentum bjet_rF_W = boost_rF_W.transform(bjet_rF);
      FourMomentum q_soft_rF_W = boost_rF_W.transform(q_soft_rF);

      double c_Whel_abs = fabs(q_soft_rF_W.p3().unit().dot(-bjet_rF_W.p3().unit())); 


      //Calculate right weights with it. Function outside analyze().

      //Weights, propabilities: 
      double rho_pos = rho(c_Whel_abs);
      double rho_neg = rho(-c_Whel_abs); 

      double prop_d_soft = rho_neg / (rho_pos + rho_neg);
      double prop_d_hard = rho_pos / (rho_pos + rho_neg);     


      //Hadronic Polarimeter (3-vector) inside t or tbar rF.
      hadrAnalyzer = (prop_d_soft * q_soft_rF.p3().unit()) + (prop_d_hard * q_hard_rF.p3().unit());    
      hadrAnalyzer = hadrAnalyzer.unit();
    
    }   
    
    
    //LEP into top-/ antitop-RF. Boost lepton in CM-ttbar.
    FourMomentum lep_CM_ttbar = boost_CM_ttbar.transform(lep);
    FourMomentum lep_rF;
    
    if (lepPart.charge() > 0) {
      lep_rF = boost_rF_top.transform(lep_CM_ttbar);
    } else {
      lep_rF = boost_rF_antitop.transform(lep_CM_ttbar);
    } 
    

    //Polarisation and correlation plots. 

    double cosTheta_k_lep, cosTheta_r_lep, cosTheta_n_lep;
    double cosThetaBar_k_lep, cosThetaBar_r_lep, cosThetaBar_n_lep;
    double cosTheta_k_hadr, cosTheta_r_hadr, cosTheta_n_hadr;
    double cosThetaBar_k_hadr, cosThetaBar_r_hadr, cosThetaBar_n_hadr;


    if (lepPart.charge() > 0) {

      cosTheta_k_lep = lep_rF.p3().unit().dot(k);
      cosTheta_k_lep = clamp(cosTheta_k_lep, -1.0, 1.0);
      _h["cosTheta_k"]->fill(cosTheta_k_lep);

      cosTheta_r_lep = lep_rF.p3().unit().dot(r);
      cosTheta_r_lep = clamp(cosTheta_r_lep, -1.0, 1.0);
      _h["cosTheta_r"]->fill(cosTheta_r_lep);

      cosTheta_n_lep = lep_rF.p3().unit().dot(n);
      cosTheta_n_lep = clamp(cosTheta_n_lep, -1.0, 1.0);
      _h["cosTheta_n"]->fill(cosTheta_n_lep); 

      cosThetaBar_k_hadr = hadrAnalyzer.dot(kbar);
      cosThetaBar_k_hadr = clamp(cosThetaBar_k_hadr, -1.0, 1.0);
      _h["cosThetaBar_k"]->fill(cosThetaBar_k_hadr);

      cosThetaBar_r_hadr = hadrAnalyzer.dot(rbar);
      cosThetaBar_r_hadr = clamp(cosThetaBar_r_hadr, -1.0, 1.0);
      _h["cosThetaBar_r"]->fill(cosThetaBar_r_hadr);

      cosThetaBar_n_hadr = hadrAnalyzer.dot(nbar);
      cosThetaBar_n_hadr = clamp(cosThetaBar_n_hadr, -1.0, 1.0);
      _h["cosThetaBar_n"]->fill(cosThetaBar_n_hadr);

      //Separated in hadr and lept decay.
      _h["cosTheta_k_lep"]->fill(cosTheta_k_lep);
      _h["cosTheta_r_lep"]->fill(cosTheta_r_lep);
      _h["cosTheta_n_lep"]->fill(cosTheta_n_lep);
      _h["cosThetaBar_k_hadr"]->fill(cosThetaBar_k_hadr);
      _h["cosThetaBar_r_hadr"]->fill(cosThetaBar_r_hadr);
      _h["cosThetaBar_n_hadr"]->fill(cosThetaBar_n_hadr);
 

      //Correlation, diagonals.
      double C_kk = cosTheta_k_lep * cosThetaBar_k_hadr;
      _h["C_kk"]->fill(C_kk);

      double C_rr = cosTheta_r_lep * cosThetaBar_r_hadr;
      _h["C_rr"]->fill(C_rr);

      double C_nn = cosTheta_n_lep * cosThetaBar_n_hadr;
      _h["C_nn"]->fill(C_nn);


      //Non-diagonals
      double C_nk = cosTheta_n_lep * cosThetaBar_k_hadr;
      _h["C_nk"]->fill(C_nk);

       double C_kn = cosTheta_k_lep * cosThetaBar_n_hadr;
      _h["C_kn"]->fill(C_kn);

      double C_nr = cosTheta_n_lep * cosThetaBar_r_hadr;
      _h["C_nr"]->fill(C_nr);

      double C_rn = cosTheta_r_lep * cosThetaBar_n_hadr;
      _h["C_rn"]->fill(C_rn);

      double C_rk = cosTheta_r_lep * cosThetaBar_k_hadr;
      _h["C_rk"]->fill(C_rk);

      double C_kr = cosTheta_k_lep * cosThetaBar_r_hadr;
      _h["C_kr"]->fill(C_kr);


    } else {


      cosThetaBar_k_lep = lep_rF.p3().unit().dot(kbar);
      cosThetaBar_k_lep = clamp(cosThetaBar_k_lep, -1.0, 1.0);
      _h["cosThetaBar_k"]->fill(cosThetaBar_k_lep);

      cosThetaBar_r_lep = lep_rF.p3().unit().dot(rbar);
      cosThetaBar_r_lep = clamp(cosThetaBar_r_lep, -1.0, 1.0);
      _h["cosThetaBar_r"]->fill(cosThetaBar_r_lep);

      cosThetaBar_n_lep = lep_rF.p3().unit().dot(nbar);
      cosThetaBar_n_lep = clamp(cosThetaBar_n_lep, -1.0, 1.0);
      _h["cosThetaBar_n"]->fill(cosThetaBar_n_lep);

      cosTheta_k_hadr = hadrAnalyzer.dot(k);
      cosTheta_k_hadr = clamp(cosTheta_k_hadr, -1.0, 1.0);
      _h["cosTheta_k"]->fill(cosTheta_k_hadr);

      cosTheta_r_hadr = hadrAnalyzer.dot(r);
      cosTheta_r_hadr = clamp(cosTheta_r_hadr, -1.0, 1.0);
      _h["cosTheta_r"]->fill(cosTheta_r_hadr);

      cosTheta_n_hadr = hadrAnalyzer.dot(n);
      cosTheta_n_hadr = clamp(cosTheta_n_hadr, -1.0, 1.0);
      _h["cosTheta_n"]->fill(cosTheta_n_hadr);


      _h["cosThetaBar_k_lep"]->fill(cosThetaBar_k_lep);
      _h["cosThetaBar_r_lep"]->fill(cosThetaBar_r_lep);
      _h["cosThetaBar_n_lep"]->fill(cosThetaBar_n_lep);
      _h["cosTheta_k_hadr"]->fill(cosTheta_k_hadr);
      _h["cosTheta_r_hadr"]->fill(cosTheta_r_hadr);
      _h["cosTheta_n_hadr"]->fill(cosTheta_n_hadr);


      //Correlation.
      double C_kk = cosThetaBar_k_lep * cosTheta_k_hadr;
      _h["C_kk"]->fill(C_kk);

      double C_rr = cosThetaBar_r_lep * cosTheta_r_hadr;
      _h["C_rr"]->fill(C_rr);

      double C_nn = cosThetaBar_n_lep * cosTheta_n_hadr;
      _h["C_nn"]->fill(C_nn);


      //Non-diagonals.
      double C_nk = cosTheta_n_hadr * cosThetaBar_k_lep;
      _h["C_nk"]->fill(C_nk);

      double C_kn = cosTheta_k_hadr * cosThetaBar_n_lep;
      _h["C_kn"]->fill(C_kn);

      double C_nr = cosTheta_n_hadr * cosThetaBar_r_lep;
      _h["C_nr"]->fill(C_nr);

      double C_rn = cosTheta_r_hadr * cosThetaBar_n_lep;
      _h["C_rn"]->fill(C_rn);

      double C_rk = cosTheta_r_hadr * cosThetaBar_k_lep;
      _h["C_rk"]->fill(C_rk);      

      double C_kr = cosTheta_k_hadr * cosThetaBar_r_lep;
      _h["C_kr"]->fill(C_kr);

     }


   }



    void finalize() {

   
      // Normalise histogram to absolute fiducial cross section or to unity
      scale(_h, crossSection() / femtobarn / sumOfWeights());
      normalize(_n, 1.0, true);

    }

    

    void dualbook(const string& name, unsigned int id) {
     // book(_h["abs_"+name],    id, 1, 1);
      book(_n["norm_"+name], 1+id, 1, 1);
    }

    void dualfill(const string& name, const double value) {
      _h["abs_"+name]->fill(value);
      _n["norm_"+name]->fill(value);
    }


    ConstGenParticlePtr findDecay(ConstGenParticlePtr particle, int ID){

      //      const int particle_ID = particle->pdg_id();
      
      ConstGenVertexPtr dv = particle->end_vertex();
      
      if (dv) {
	
	for(ConstGenParticlePtr pp: HepMCUtils::particles(dv, Relatives::CHILDREN)){

	  if(pp->pdg_id() == ID){

	    return findDecay(pp, ID);

	  }
	  else{
	    
	    return particle;
	      
	  }

	  
	}
	
      }      
      
    }
    
    double findParent(const Particle particle) {
      
      const int particle_ID = particle.pid();
      
      Particle current = particle;

      int counter = 0;
      
      while(current.pid() == particle_ID){
	
	for (const Particle& parent : current.parents()) {
	  
	  current = parent;
	  
	}
	counter += 1;

	if(counter == 15)
	  return current.pid();
	
      }
      return current.pid();
    }

    bool fromTop(const Particle particle) {

      //      const int particle_ID = particle.pid();

      Particle current = particle;

      //      std::cout << "check from top = " << particle_ID << std::endl;
      
      for (const Particle& parent : current.parents()) {
	
	current = parent;
	
	if(fabs(parent.pid()) == 6)
	    return true;
	else
	  fromTop(parent);
	  
	
      }

    return false;
  }
    
    
    
  private:

    // Compute z component of neutrino momentum given lepton and met
    double _computeneutrinoz(const FourMomentum& lep, const FourMomentum& met) const {
      double m_W = 80.399; // in GeV, given in the paper
      double k = (( sqr( m_W ) - sqr( lep.mass() ) ) / 2 ) + (lep.px() * met.px() + lep.py() * met.py());
      double a = sqr ( lep.E() ) - sqr ( lep.pz() );
      double b = -2*k*lep.pz();
      double c = sqr( lep.E() ) * sqr( met.pT() ) - sqr( k );
      double discriminant = sqr(b) - 4 * a * c;
      double quad[2] = { (- b - sqrt(discriminant)) / (2 * a), (- b + sqrt(discriminant)) / (2 * a) }; //two possible quadratic solns

      double pzneutrino;
      if (discriminant < 0) { // if the discriminant is negative:
        pzneutrino = - b / (2 * a);
      } else { // if the discriminant is positive, take the soln with smallest absolute value
        pzneutrino = (fabs(quad[0]) < fabs(quad[1])) ? quad[0] : quad[1];
      }
      return pzneutrino;
    }


    //Weights hadronic polarimeter:
      double rho(double c) const {
        double f_R = 0.0;
        double f_0 = 0.7;
        double f_L = 0.3;

        return (3.0/8.0) * f_R * (1+c)*(1+c) + (3.0/4.0) * f_0 * (1-c*c) + (3.0/8.0) * f_L * (1-c)*(1-c);
      }


    /// @name Histograms
    map<string, Histo1DPtr> _h, _n;

    size_t _mode;

    int eventCounter;
    
  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TEST_ROUTINE);
}
