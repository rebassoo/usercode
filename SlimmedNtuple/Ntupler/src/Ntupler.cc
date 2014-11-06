#include "SlimmedNtuple/Ntupler/interface/Ntupler.h"

//
// constructors and destructor
//
Ntupler::Ntupler(const edm::ParameterSet& iConfig)

{

  isMC = iConfig.getParameter<bool>("ismc");
  isoValInputTags = iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags");

  muon_id_ = new std::vector<bool>;
  muon_pt_ = new std::vector<float>;
  muon_px_ = new std::vector<float>;
  muon_py_ = new std::vector<float>;
  muon_pz_ = new std::vector<float>;
  muon_et_ = new std::vector<float>;
  muon_energy_ = new std::vector<float>;
  muon_charge_ = new std::vector<float>;
  muon_eta_ = new std::vector<float>;
  muon_phi_ = new std::vector<float>;
  muon_track_pt_ = new std::vector<float>;
  muon_track_eta_ = new std::vector<float>;
  muon_track_phi_ = new std::vector<float>;


  gen_pt_ = new std::vector<float>;
  gen_px_ = new std::vector<float>;
  gen_py_ = new std::vector<float>;
  gen_pz_ = new std::vector<float>;
  gen_et_ = new std::vector<float>;
  gen_energy_ = new std::vector<float>;
  gen_charge_ = new std::vector<float>;
  gen_eta_ = new std::vector<float>;
  gen_phi_ = new std::vector<float>;
  gen_pdgId_ = new std::vector<int>;


  electron_id_ = new std::vector<bool>;
  electron_pt_ = new std::vector<float>;
  electron_px_ = new std::vector<float>;
  electron_py_ = new std::vector<float>;
  electron_pz_ = new std::vector<float>;
  electron_et_ = new std::vector<float>;
  electron_energy_ = new std::vector<float>;
  electron_charge_ = new std::vector<float>;
  electron_eta_ = new std::vector<float>;
  electron_phi_ = new std::vector<float>;
  electron_track_pt_ = new std::vector<float>;
  electron_track_eta_ = new std::vector<float>;
  electron_track_phi_ = new std::vector<float>;

  trigger_prescaleValue_ = new std::vector<int>;
  trigger_name_ = new std::vector<std::string>;
  trigger_decision_ = new std::vector<int>;
  vertex_nextra_tracks_ = new std::vector<int>;
  vertexCandMuIndex_ = new std::vector<int>;
  vertexCandEIndex_ = new std::vector<int>;

  vertex_nextra_tracks_ee_ = new std::vector<int>;
  vertexCandE1Index_ee_ = new std::vector<int>;
  vertexCandE2Index_ee_ = new std::vector<int>;

  vertex_nextra_tracks_mumu_ = new std::vector<int>;
  vertexCandMu1Index_mumu_ = new std::vector<int>;
  vertexCandMu2Index_mumu_ = new std::vector<int>;

  met_et_ = new std::vector<float>;
  met_px_ = new std::vector<float>;
  met_py_ = new std::vector<float>;
  met_sumEt_ = new std::vector<float>;

  vertex_ntracks_ = new std::vector<int>;
  vertex_x_ = new std::vector<float>;
  vertex_y_ = new std::vector<float>;
  vertex_z_ = new std::vector<float>;
  vertexerr_x_ = new std::vector<float>;
  vertexerr_y_ = new std::vector<float>;
  vertexerr_z_ = new std::vector<float>;
  vertex_chi2_ = new std::vector<float>;
  vertex_ndof_ = new std::vector<float>;
  vertex_isFake_ = new std::vector<int>;
  vertex_isValid_ = new std::vector<int>;
  vertex_tracks_pt_ = new std::vector<std::vector<float>>;
  vertex_tracks_eta_ = new std::vector<std::vector<float>>;
  vertex_tracks_phi_ = new std::vector<std::vector<float>>;
  vertex_tracks_d0_ = new std::vector<std::vector<float>>;
  vertex_tracks_d0Err_ = new std::vector<std::vector<float>>;
  vertex_tracks_dz_ = new std::vector<std::vector<float>>;
  vertex_tracks_dzErr_ = new std::vector<std::vector<float>>;
  vertex_tracks_ptErr_ = new std::vector<std::vector<float>>;
  vertex_tracks_highPurity_ = new std::vector<std::vector<int>>;
  vertex_tracks_weight_ = new std::vector<std::vector<float>>;
  vertex_tracks_dxy_ = new std::vector<std::vector<float>>;
  vertex_tracks_dxyErr_ = new std::vector<std::vector<float>>;
  vertex_tracks_normChi2_ = new std::vector<std::vector<float>>;
  vertex_tracks_pxlhits_ = new std::vector<std::vector<float>>;
  vertex_tracks_silcnhits_ = new std::vector<std::vector<float>>;


  tracks_chi2_ = new std::vector<float>;
  tracks_ndof_ = new std::vector<float>;
  tracks_charge_ = new std::vector<float>;
  tracks_pt_ = new std::vector<float>;
  tracks_px_ = new std::vector<float>;
  tracks_py_ = new std::vector<float>;
  tracks_pz_ = new std::vector<float>;
  tracks_eta_ = new std::vector<float>;
  tracks_phi_ = new std::vector<float>;
  tracks_d0_ = new std::vector<float>;
  tracks_d0Err_ = new std::vector<float>;
  tracks_pxlhits_ = new std::vector<float>;
  tracks_silcnhits_ = new std::vector<float>;
  tracks_dz_ = new std::vector<float>;
  tracks_vx_ = new std::vector<float>;
  tracks_vy_ = new std::vector<float>;
  tracks_vz_ = new std::vector<float>;
  tracks_sim_pt_ = new std::vector<float>;
  tracks_sim_eta_ = new std::vector<float>;
  tracks_sim_phi_ = new std::vector<float>;
  tracks_sim_id_ = new std::vector<float>;
  tracks_sim_motherid_ = new std::vector<float>;
  tracks_gen_id_ = new std::vector<float>;
  tracks_gen_motherid_ = new std::vector<float>;



  ev_ = new uint;
  run_ = new uint;
  lumiblock_ = new uint;
  experimentType_ = new uint;
  bunchCrossing_ = new uint;
  orbitNumber_ = new uint;
  Tnpv_ = new float;
  pileupWeight_ = new float;

  edm::Service<TFileService> fs;

  tree_=fs->make<TTree>("SlimmedNtuple","SlimmedNtuple");

  h_trueNumInteractions = fs->make<TH1F>("h_trueNumInteractions" , "PU" , 50 , -0.5 , 49.5 );
  h_trueNumInteractions0 = fs->make<TH1F>("h_trueNumInteractions0" , "PU0" , 50 , -0.5 , 49.5 );
  h_pt_extra_track = fs->make<TH1F>("h_pt_extra_track" , "ExtraTracks" , 5000 , 0 , 500 );

  if(isMC){
    LumiWeights = new edm::LumiReWeighting(
					   //  LumiWeights = new edm::LumiReWeighting(
					   //"/uscms/home/rebassoo/2013_04_12_MakeMCPileupDist/RunningCrab/PileupMC_total.root",
					   "PileupMC_total.root",
					   //"/uscms/home/rebassoo/2013_04_12_MakeMCPileupDist/RunningCrab/Data_8TeV_PileupDistribution.root",
					   //"/uscms/home/rebassoo/2013_04_12_MakeMCPileupDist/RunningCrab/Data_8TeV_PileupDistribution.root",
					   "Data_8TeV_PileupDistribution.root",
					   "demo/h_trueNumInteractions0",
					   //"pileup",
					   "pileup"
					   );
    //  LumiWeights->weight3D_init(1.0);
  }

   //now do what ever initialization is needed


  tree_->Branch("run",run_,"run/i");
  tree_->Branch("event",ev_,"event/i");
  tree_->Branch("lumiblock",lumiblock_,"lumiblock/i");
  tree_->Branch("experimentType",experimentType_,"experimentType/i");
  tree_->Branch("bunchCrossing",bunchCrossing_,"bunchCrossing/i");
  tree_->Branch("orbitNumber",orbitNumber_,"orbitNumber/i");
  tree_->Branch("Tnpv",Tnpv_,"Tnpv/f");
  tree_->Branch("pileupWeight",pileupWeight_,"pileupWeight/f");
  tree_->Branch("muon_id",&muon_id_);
  tree_->Branch("muon_pt",&muon_pt_);
  tree_->Branch("muon_px",&muon_px_);
  tree_->Branch("muon_py",&muon_py_);
  tree_->Branch("muon_pz",&muon_pz_);
  tree_->Branch("muon_et",&muon_et_);
  tree_->Branch("muon_energy",&muon_energy_);
  tree_->Branch("muon_charge",&muon_charge_);
  tree_->Branch("muon_eta",&muon_eta_);
  tree_->Branch("muon_phi",&muon_phi_);
  tree_->Branch("muon_track_pt",&muon_track_pt_);
  tree_->Branch("muon_track_eta",&muon_track_eta_);
  tree_->Branch("muon_track_phi",&muon_track_phi_);

  tree_->Branch("gen_pt",&gen_pt_);
  tree_->Branch("gen_px",&gen_px_);
  tree_->Branch("gen_py",&gen_py_);
  tree_->Branch("gen_pz",&gen_pz_);
  tree_->Branch("gen_et",&gen_et_);
  tree_->Branch("gen_energy",&gen_energy_);
  tree_->Branch("gen_charge",&gen_charge_);
  tree_->Branch("gen_eta",&gen_eta_);
  tree_->Branch("gen_phi",&gen_phi_);
  tree_->Branch("gen_pdgId",&gen_pdgId_);


  tree_->Branch("electron_id",&electron_id_);
  tree_->Branch("electron_pt",&electron_pt_);
  tree_->Branch("electron_px",&electron_px_);
  tree_->Branch("electron_py",&electron_py_);
  tree_->Branch("electron_pz",&electron_pz_);
  tree_->Branch("electron_et",&electron_et_);
  tree_->Branch("electron_energy",&electron_energy_);
  tree_->Branch("electron_charge",&electron_charge_);
  tree_->Branch("electron_eta",&electron_eta_);
  tree_->Branch("electron_phi",&electron_phi_);
  tree_->Branch("electron_track_pt",&electron_track_pt_);
  tree_->Branch("electron_track_eta",&electron_track_eta_);
  tree_->Branch("electron_track_phi",&electron_track_phi_);

  tree_->Branch("trigger_prescaleValue",&trigger_prescaleValue_);
  tree_->Branch("trigger_name",&trigger_name_);
  tree_->Branch("trigger_decision",&trigger_decision_);
  tree_->Branch("vertex_nextra_tracks",&vertex_nextra_tracks_);
  tree_->Branch("vertexCandMuIndex",&vertexCandMuIndex_);
  tree_->Branch("vertexCandEIndex",&vertexCandEIndex_);

  tree_->Branch("vertex_nextra_tracks_ee",&vertex_nextra_tracks_ee_);
  tree_->Branch("vertexCandE1Index_ee",&vertexCandE1Index_ee_);
  tree_->Branch("vertexCandE2Index_ee",&vertexCandE2Index_ee_);

  tree_->Branch("vertex_nextra_tracks_mumu",&vertex_nextra_tracks_mumu_);
  tree_->Branch("vertexCandMu1Index_mumu",&vertexCandMu1Index_mumu_);
  tree_->Branch("vertexCandMu2Index_mumu",&vertexCandMu2Index_mumu_);

  tree_->Branch("met_et",&met_et_);
  tree_->Branch("met_px",&met_px_);
  tree_->Branch("met_py",&met_py_);
  tree_->Branch("met_sumEt",&met_sumEt_);

  tree_->Branch("vertex_ntracks",&vertex_ntracks_);
  tree_->Branch("vertex_x",&vertex_x_);
  tree_->Branch("vertex_y",&vertex_y_);
  tree_->Branch("vertex_z",&vertex_z_);
  tree_->Branch("vertexerr_x",&vertexerr_x_);
  tree_->Branch("vertexerr_y",&vertexerr_y_);
  tree_->Branch("vertexerr_z",&vertexerr_z_);
  tree_->Branch("vertex_chi2",&vertex_chi2_);
  tree_->Branch("vertex_ndof",&vertex_ndof_);
  tree_->Branch("vertex_isFake",&vertex_isFake_);
  tree_->Branch("vertex_isValid",&vertex_isValid_);
  tree_->Branch("vertex_tracks_pt",&vertex_tracks_pt_);
  tree_->Branch("vertex_tracks_eta",&vertex_tracks_eta_);
  tree_->Branch("vertex_tracks_phi",&vertex_tracks_phi_);
  tree_->Branch("vertex_tracks_d0",&vertex_tracks_d0_);
  tree_->Branch("vertex_tracks_d0Err",&vertex_tracks_d0Err_);
  tree_->Branch("vertex_tracks_dz",&vertex_tracks_dz_);
  tree_->Branch("vertex_tracks_dzErr",&vertex_tracks_dzErr_);
  tree_->Branch("vertex_tracks_ptErr",&vertex_tracks_ptErr_);
  tree_->Branch("vertex_tracks_highPurity",&vertex_tracks_highPurity_);
  tree_->Branch("vertex_tracks_weight",&vertex_tracks_weight_);
  tree_->Branch("vertex_tracks_dxy_",&vertex_tracks_dxy_);
  tree_->Branch("vertex_tracks_dxyErr_",&vertex_tracks_dxyErr_);
  tree_->Branch("vertex_tracks_normChi2_",&vertex_tracks_normChi2_);
  tree_->Branch("vertex_tracks_pxlhits_",&vertex_tracks_pxlhits_);
  tree_->Branch("vertex_tracks_silcnhits_",&vertex_tracks_silcnhits_);
     

  tree_->Branch("tracks_chi2",&tracks_chi2_);
  tree_->Branch("tracks_ndof",&tracks_ndof_);
  tree_->Branch("tracks_charge",&tracks_charge_);
  tree_->Branch("tracks_pt",&tracks_pt_);
  tree_->Branch("tracks_px",&tracks_px_);
  tree_->Branch("tracks_py",&tracks_py_);
  tree_->Branch("tracks_pz",&tracks_pz_);
  tree_->Branch("tracks_eta",&tracks_eta_);
  tree_->Branch("tracks_phi",&tracks_phi_);
  tree_->Branch("tracks_d0",&tracks_d0_);
  tree_->Branch("tracks_d0Err",&tracks_d0Err_);
  tree_->Branch("tracks_pxlhits_",&tracks_pxlhits_);
  tree_->Branch("tracks_silcnhits_",&tracks_silcnhits_);
  tree_->Branch("tracks_dz",&tracks_dz_);
  tree_->Branch("tracks_vx",&tracks_vx_);
  tree_->Branch("tracks_vy",&tracks_vy_);
  tree_->Branch("tracks_vz",&tracks_vz_);
  tree_->Branch("tracks_sim_pt",&tracks_sim_pt_);
  tree_->Branch("tracks_sim_eta",&tracks_sim_eta_);
  tree_->Branch("tracks_sim_phi",&tracks_sim_phi_);
  tree_->Branch("tracks_sim_id",&tracks_sim_id_);
  tree_->Branch("tracks_sim_motherid",&tracks_sim_motherid_);
  tree_->Branch("tracks_gen_id",&tracks_gen_id_);
  tree_->Branch("tracks_gen_motherid",&tracks_gen_motherid_);


}


Ntupler::~Ntupler()
{
  delete ev_;
  delete run_;
  delete lumiblock_;
  delete experimentType_;
  delete bunchCrossing_;
  delete orbitNumber_;
  delete Tnpv_;
  delete pileupWeight_;
  delete muon_id_;
  delete muon_pt_;
  delete muon_px_;
  delete muon_py_;
  delete muon_pz_;
  delete muon_et_;
  delete muon_energy_;
  delete muon_charge_;
  delete muon_eta_;
  delete muon_phi_;
  delete muon_track_pt_;
  delete muon_track_eta_;
  delete muon_track_phi_;

  delete gen_pt_;
  delete gen_px_;
  delete gen_py_;
  delete gen_pz_;
  delete gen_et_;
  delete gen_energy_;
  delete gen_charge_;
  delete gen_eta_;
  delete gen_phi_;
  delete gen_pdgId_;


  delete electron_id_;
  delete electron_pt_;
  delete electron_px_;
  delete electron_py_;
  delete electron_pz_;
  delete electron_et_;
  delete electron_energy_;
  delete electron_charge_;
  delete electron_eta_;
  delete electron_phi_;
  delete electron_track_pt_;
  delete electron_track_eta_;
  delete electron_track_phi_;


  delete trigger_prescaleValue_;
  delete trigger_name_;
  delete trigger_decision_;
  delete vertex_nextra_tracks_;
  delete vertexCandMuIndex_;
  delete vertexCandEIndex_;

  delete vertex_nextra_tracks_ee_;
  delete vertexCandE1Index_ee_;
  delete vertexCandE2Index_ee_;

  delete vertex_nextra_tracks_mumu_;
  delete vertexCandMu1Index_mumu_;
  delete vertexCandMu2Index_mumu_;

  delete met_et_;
  delete met_px_;
  delete met_py_;
  delete met_sumEt_;

  delete vertex_ntracks_;
  delete vertex_x_;
  delete vertex_y_;
  delete vertex_z_;
  delete vertexerr_x_;
  delete vertexerr_y_;
  delete vertexerr_z_;
  delete vertex_chi2_;
  delete vertex_ndof_;
  delete vertex_isFake_;
  delete vertex_isValid_;
  delete vertex_tracks_pt_;
  delete vertex_tracks_eta_;
  delete vertex_tracks_phi_;
  delete vertex_tracks_highPurity_;
  delete vertex_tracks_weight_;
  delete vertex_tracks_d0_;
  delete vertex_tracks_d0Err_;
  delete vertex_tracks_dz_;
  delete vertex_tracks_dzErr_;
  delete vertex_tracks_ptErr_;
  delete vertex_tracks_dxy_;
  delete vertex_tracks_dxyErr_;
  delete vertex_tracks_normChi2_;
  delete vertex_tracks_pxlhits_;
  delete vertex_tracks_silcnhits_;

  delete tracks_chi2_;
  delete tracks_ndof_;
  delete tracks_charge_;
  delete tracks_pt_;
  delete tracks_px_;
  delete tracks_py_;
  delete tracks_pz_;
  delete tracks_eta_;
  delete tracks_phi_;
  delete tracks_d0_;
  delete tracks_d0Err_;
  delete tracks_pxlhits_;
  delete tracks_silcnhits_;
  delete tracks_dz_;
  delete tracks_vx_;
  delete tracks_vy_;
  delete tracks_vz_;
  delete tracks_sim_pt_;
  delete tracks_sim_eta_;
  delete tracks_sim_phi_;
  delete tracks_sim_id_;
  delete tracks_sim_motherid_;
  delete tracks_gen_id_;
  delete tracks_gen_motherid_;




//    // do anything here that needs to be done at desctruction time
//    // (e.g. close files, deallocate resources etc.)

 }


//
// member functions
//

// ------------ method called for each event  ------------
void
Ntupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;


   ////////////////////////////////////////////////////////////////////////////////
   //Fill basic event variables
   ////////////////////////////////////////////////////////////////////////////////
   //Event information
   *run_ = iEvent.id().run();
   *ev_ = iEvent.id().event();
   *lumiblock_ = iEvent.luminosityBlock();
   *experimentType_ = iEvent.experimentType();
   *bunchCrossing_ = iEvent.bunchCrossing();
   *orbitNumber_ = iEvent.orbitNumber();
   *Tnpv_=-999.;
   *pileupWeight_=1.;
   int nextra_tracks = 999;
   int nextra_tracks_mumu = 999;
   int nextra_tracks_ee = 999;


   ////////////////////////////////////////////////////////////////////////////////
   //Look at basic pileup information and plot true num interactions, only for MC
   ////////////////////////////////////////////////////////////////////////////////
   if(isMC){
     
     edm::InputTag PileupSrc_("addPileupInfo");
     Handle<std::vector< PileupSummaryInfo > >  PupInfo;
     iEvent.getByLabel(PileupSrc_, PupInfo);
     std::vector<PileupSummaryInfo>::const_iterator PVI;
     for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
       //     std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
       //     std::cout << " True Num Interactions: " << PVI->getTrueNumInteractions() << endl;
       h_trueNumInteractions->Fill(PVI->getTrueNumInteractions());
       if(PVI->getBunchCrossing()==0){     h_trueNumInteractions0->Fill(PVI->getTrueNumInteractions());
	 *Tnpv_=PVI->getTrueNumInteractions();
	 *pileupWeight_ = LumiWeights->weight( PVI->getTrueNumInteractions() );
       }
     }
   }





   ////////////////////////////////////////////////////////////////////////////////
   //Fill MET information into ntuple
   ////////////////////////////////////////////////////////////////////////////////
   Handle< edm::View<reco::PFMET> > metHandle;
   iEvent.getByLabel("pfType1CorrectedMet",metHandle);
   edm::View<reco::PFMET> mets = *metHandle;
   edm::View<reco::PFMET>::const_iterator iMet;

   for (iMet = mets.begin(); iMet != mets.end(); ++iMet) {
     (*met_et_).push_back(iMet->et());
     (*met_px_).push_back(iMet->px());
     (*met_py_).push_back(iMet->py());
     (*met_sumEt_).push_back(iMet->sumEt());
   }


   ////////////////////////////////////////////////////////////////////////////////
   //Fill Trigger information into ntuple
   ////////////////////////////////////////////////////////////////////////////////
   Handle<TriggerResults> hltResults;
   iEvent.getByLabel(InputTag("TriggerResults","","HLT"),hltResults);
   const TriggerNames & trigNames = iEvent.triggerNames(*hltResults);

   for(unsigned int i=0; i<trigNames.size();i++){
     (*trigger_prescaleValue_).push_back(hltConfig_.prescaleValue(iEvent,iSetup,trigNames.triggerName(i)));
     (*trigger_name_).push_back(trigNames.triggerName(i));
     (*trigger_decision_).push_back(hltResults->accept(i));
   }

   if(isMC){

     ////////////////////////////////////////////////////////////////////////////////
     //If isMC fill Gen particle variables in ntuple 
     ////////////////////////////////////////////////////////////////////////////////
     
     Handle<reco::GenParticleCollection> genP;
     iEvent.getByLabel("genParticles",genP);
     //     cout<<" "<<endl;
     //     cout<<"Just before gen particles"<<endl;
     //Loop over gen particles and store particles that come from a W in the documentation line with pt and eta cuts
     /*
     double countMus=0;
     bool dimuon=false;
     for (reco::GenParticleCollection::const_iterator mcIter=genP->begin(); mcIter != genP->end(); mcIter++ ) {
       if(mcIter->pdgId()==13&& mcIter->status() == 3){countMus++;}
       if(mcIter->pdgId()==-13&& mcIter->status() == 3){countMus++;}
     }
     if(countMus==2){dimuon=true;}
     */
     for (reco::GenParticleCollection::const_iterator mcIter=genP->begin(); mcIter != genP->end(); mcIter++ ) {
       /*
       cout<<"MC id is: "<<mcIter->pdgId()<<endl;
       cout<<"Pt, eta, phi is: "<<mcIter->pt()<<", "<<mcIter->eta()<<", "<<mcIter->phi()<<endl;
       cout<<"E, px,py, pz is: "<<mcIter->energy()<<", "<<mcIter->px()<<", "<<mcIter->py()<<", "<<mcIter->pz()<<endl;
       cout<<"Status is: "<<mcIter->status()<<endl;
       */
       /*
       if(dimuon){

       */       
       //       if (mcIter->pt()>20 && fabs(mcIter->eta())<2.5){ 
       if ( (fabs(mcIter->pdgId())==11|| fabs(mcIter->pdgId())==12 || fabs(mcIter->pdgId())==13 || fabs(mcIter->pdgId())==14 || fabs(mcIter->pdgId())==15 || fabs(mcIter->pdgId())==16 ) && mcIter->status() == 3 ){ 

       //       if ( (fabs(mcIter->pdgId())==11|| fabs(mcIter->pdgId())==12 || fabs(mcIter->pdgId())==13 || fabs(mcIter->pdgId())==14 || fabs(mcIter->pdgId())==15 || fabs(mcIter->pdgId())==16 ) ){ 
       //       if(mcIter->pdgId()==15){
	 //cout<<"Status is: "<<mcIter->status();
	 //cout<<", MC id is: "<<mcIter->pdgId()<<endl;
	 //cout<<"Pt, eta, phi is: "<<mcIter->pt()<<", "<<mcIter->eta()<<", "<<mcIter->phi()<<endl;
	 //const Candidate * mom = p.mother();
	 //	 cout<<", MC id is: "<<mcIter->pdgId()<<endl;

	 //cout<<"MC id is: "<<mcIter->pdgId()<<endl;
	 //	 cout<<"Pt, eta, phi is: "<<mcIter->pt()<<", "<<mcIter->eta()<<", "<<mcIter->phi()<<endl;
	 //cout<<"Status is: "<<mcIter->status()<<endl;

	 /*
	 int n = mcIter->numberOfDaughters();
	 for(int j = 0; j < n; ++ j) {
	   const reco::Candidate * d = mcIter->daughter( j );
	   int dauId = d->pdgId();
	   cout<<"Daughter pdg Id: "<<dauId<<endl;
	 }
	 */
	 (*gen_pt_).push_back(mcIter->pt());
	 (*gen_px_).push_back(mcIter->px());
	 (*gen_py_).push_back(mcIter->py());
	 (*gen_pz_).push_back(mcIter->pz());
	 (*gen_et_).push_back(mcIter->et());
	 (*gen_energy_).push_back(mcIter->energy());
	 (*gen_charge_).push_back(mcIter->charge());
	 (*gen_eta_).push_back(mcIter->eta());
	 (*gen_phi_).push_back(mcIter->phi());
	 (*gen_pdgId_).push_back(mcIter->pdgId());
	 //       }// end of looking at documentation line	   
       }
       //       }//end of if statement to make sure particle has a mother
       //}//end of dimuon      
     }//end of looking at gen particles
     
     //     cout<<"End of looking at gen particles"<<endl;
     //     cout<<" "<<endl;

   }//end of specifying that it needs to be MC
   
   // beam spot
   edm::Handle<reco::BeamSpot> beamspot_h;
   //   iEvent.getByLabel(beamSpotInputTag_, beamspot_h);
   iEvent.getByLabel("offlineBeamSpot", beamspot_h);
   const reco::BeamSpot &beamSpot = *(beamspot_h.product());

   // vertices
   edm::Handle<reco::VertexCollection> vtxs;
   iEvent.getByLabel("offlinePrimaryVertices", vtxs);
   const reco::VertexCollection* vertexs = vtxs.product();
   reco::VertexCollection::const_iterator vertex_i;


   ////////////////////////////////////////////////////////////////////////////////
   //Muon Selection and filling of ntuple
   ////////////////////////////////////////////////////////////////////////////////

   edm::Handle< edm::View<reco::Muon> > muonHandle;
   iEvent.getByLabel("muons", muonHandle);
   edm::View<reco::Muon> muons = *muonHandle;
   edm::View<reco::Muon>::const_iterator iMuon;

   //loop over muons and store muons that pass id cuts, including pt>20 GeV and |eta|<2.4
   int n_MuonsPassingCuts = 0;
   vector<int> muonCharge;
   vector<double> muonTrackPt;
   vector<double> muonTrackEta;
   vector<double> muonTrackPhi;
   for (iMuon = muons.begin(); iMuon != muons.end(); ++iMuon) {
     bool globalMuon = false;
     if(iMuon->isGlobalMuon()){globalMuon=true;}

     bool pfMuon = false;
     if(iMuon->isPFMuon()){pfMuon=true;}

     bool chi2 = false;
     bool MuonChamberHits = false;
     if(iMuon->globalTrack().isNonnull()){
       if(iMuon->globalTrack()->normalizedChi2()< 10){chi2=true;}
       if(iMuon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0){MuonChamberHits=true;}
     }

     bool MuonStations = false;
     if(iMuon->numberOfMatchedStations() > 1){MuonStations=true;}

     bool NPxlHits = false;
     bool NtrackerLayers = false;     
     if(iMuon->innerTrack().isNonnull()){
       if(iMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0){NPxlHits = true;}
       if(iMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5){NtrackerLayers = true;}
       //       if( fabs(iMuon->innerTrack()->dxy(beamSpot.position())) < 0.2){d0cut = true;}
       
     }

     bool d0cut = false;
     bool dzcut = false;
     double d0vtx = 999.;
     double dzvtx = 999.;
     if(iMuon->muonBestTrack().isNonnull()){
       if (vtxs->size() > 0) {
	 reco::VertexRef vtx(vtxs, 0);    
	 d0vtx = iMuon->muonBestTrack()->dxy(vtx->position());
	 dzvtx = iMuon->muonBestTrack()->dz(vtx->position());
       } 
       else {
	 d0vtx = iMuon->muonBestTrack()->dxy();
	 dzvtx = iMuon->muonBestTrack()->dz();
       }
       if( fabs(d0vtx) < 0.2){d0cut = true;}
       if( fabs(dzvtx) < 0.5){dzcut = true;}
     }

     bool pass_muonId =false;

     if(iMuon->pt()>20&&fabs(iMuon->eta())<2.4&&d0cut&&dzcut&&NtrackerLayers&&NPxlHits&&MuonStations&&MuonChamberHits&&chi2&&pfMuon&&globalMuon){
       pass_muonId=true;
     }
     /*
     cout<<"Muon pt, eta"<<iMuon->pt()<<", "<<iMuon->eta()<<endl;

     if(d0cut){
       cout<<"pass d0 cut"<<endl;
     }

     if(dzcut){
       cout<<"pass dz cut"<<endl;
     }

     if(NtrackerLayers){
       cout<<"pass NTrackerLayers cut"<<endl;
     }

     if(NPxlHits){
       cout<<"pass NPxlHits cut"<<endl;
     }

     if(chi2){
       cout<<"pass chi2 cut"<<endl;
     }

     if(MuonStations){
       cout<<"pass Muon Stations cut"<<endl;}

     if(MuonChamberHits){
       cout<<"pass Muon ChamberHits cut"<<endl;}

     if(pfMuon){
       cout<<"pass pfMuon"<<endl;}

     if(globalMuon){
       cout<<"pass globalMuon"<<endl;}
     */
     if(iMuon->pt()>20&&fabs(iMuon->eta())<2.4){
       
       (*muon_id_).push_back(pass_muonId);
       (*muon_pt_).push_back(iMuon->pt());
       (*muon_px_).push_back(iMuon->px());
       (*muon_py_).push_back(iMuon->py());
       (*muon_pz_).push_back(iMuon->pz());
       (*muon_et_).push_back(iMuon->et());
       (*muon_energy_).push_back(iMuon->energy());
       (*muon_charge_).push_back(iMuon->charge());
       (*muon_eta_).push_back(iMuon->eta());
       (*muon_phi_).push_back(iMuon->phi());
       if(iMuon->innerTrack().isNonnull()){
	 muonTrackPt.push_back(iMuon->innerTrack()->pt());
	 muonTrackEta.push_back(iMuon->innerTrack()->eta());
	 muonTrackPhi.push_back(iMuon->innerTrack()->phi());
	 (*muon_track_pt_).push_back(iMuon->innerTrack()->pt());
	 (*muon_track_eta_).push_back(iMuon->innerTrack()->eta());
	 (*muon_track_phi_).push_back(iMuon->innerTrack()->phi());
       }
       else{
	 muonTrackPt.push_back(-9999.);
	 muonTrackEta.push_back(-9999.);
	 muonTrackPhi.push_back(-9999.);	 
	 (*muon_track_pt_).push_back(-9999.);
	 (*muon_track_eta_).push_back(-9999.);
	 (*muon_track_phi_).push_back(-9999.);
	 cout<<"There is something wrong with the Muon track"<<endl;}
       muonCharge.push_back(iMuon->charge());
       n_MuonsPassingCuts++;
     }
     
   }// end of loop over muons
        


   ////////////////////////////////////////////////////////////////////////////////
   //Electron Selection and filling of ntuple
   ////////////////////////////////////////////////////////////////////////////////

   // New 2012 electron ID variables
   // conversions 
   edm::Handle<reco::ConversionCollection> conversions_h; 
   iEvent.getByLabel("allConversions", conversions_h); 

   // iso deposits 
   IsoDepositVals isoVals(isoValInputTags.size()); 
   for (size_t j = 0; j < isoValInputTags.size(); ++j) { 
     iEvent.getByLabel(isoValInputTags[j], isoVals[j]); 
   } 

   // rho for isolation 
   edm::Handle<double> rhoIso_h; 
   iEvent.getByLabel("kt6PFJetsForIsolation","rho", rhoIso_h); 
   double rhoIso = *(rhoIso_h.product()); 

   // electrons
   edm::Handle<reco::GsfElectronCollection> els_h;
   iEvent.getByLabel("gsfElectrons", els_h);
   
   double checkduplicates[4] = {0,0,0,0};


   unsigned int n = els_h->size();
   int n_ElsPassingCuts = 0;
   vector<int> electronCharge;
   vector<double> electronTrackPt;
   vector<double> electronTrackEta;
   vector<double> electronTrackPhi;

   for(unsigned int i = 0; i < n; ++i) {

     // get reference to electron
     reco::GsfElectronRef ele(els_h, i);

     float pt            = ele->pt();
     float eta           = ele->superCluster()->eta();

     double iso_ch = (*(isoVals)[0])[ele];
     double iso_em = (*(isoVals)[1])[ele];
     double iso_nh = (*(isoVals)[2])[ele];
     bool passmediumEid=false;
     if(PassTriggerCuts(EgammaCutBasedEleId::TRIGGERTIGHT, ele ) == true){

       passmediumEid = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, ele, conversions_h, beamSpot, vtxs, iso_ch, iso_em, iso_nh, rhoIso);
     }

     //     if(passmediumEid&&pt>20&&fabs(eta)<2.4){
     bool electronId = false;
     if(passmediumEid){electronId = true;}
     if(pt>20&&fabs(eta)<2.4){
       
       bool duplicate = false;
       if(pt==checkduplicates[0] && 
	  checkduplicates[1]==eta && 
	  checkduplicates[2]== ele->phi()&&
	  checkduplicates[3]==ele->charge()){
	 cout<<"Duplicate Electron"<<endl;
	 duplicate = true;
       }
       else{
	 checkduplicates[0] = pt;checkduplicates[1]=eta;
	 checkduplicates[2] = ele->phi(); checkduplicates[3]=ele->charge();
       }

       
       if(!duplicate){
	 
	 (*electron_id_).push_back(electronId);
	 (*electron_pt_).push_back(ele->pt());
	 (*electron_px_).push_back(ele->px());
	 (*electron_py_).push_back(ele->py());
	 (*electron_pz_).push_back(ele->pz());
	 (*electron_et_).push_back(ele->et());
	 (*electron_energy_).push_back(ele->energy());
	 (*electron_charge_).push_back(ele->charge());
	 (*electron_eta_).push_back(ele->eta());
	 (*electron_phi_).push_back(ele->phi());
	 n_ElsPassingCuts++;
	 electronCharge.push_back(ele->charge());
	 if(ele->closestCtfTrackRef().isNonnull()){
	   electronTrackPt.push_back(ele->closestCtfTrackRef()->pt());
	   electronTrackEta.push_back(ele->closestCtfTrackRef()->eta());
	   electronTrackPhi.push_back(ele->closestCtfTrackRef()->phi());
	   (*electron_track_pt_).push_back(ele->closestCtfTrackRef()->pt());
	   (*electron_track_eta_).push_back(ele->closestCtfTrackRef()->eta());
	   (*electron_track_phi_).push_back(ele->closestCtfTrackRef()->phi());
	 }
	 else{
	   electronTrackPt.push_back(-9999.);
	   electronTrackEta.push_back(-9999.);
	   electronTrackPhi.push_back(-9999.);
	   (*electron_track_pt_).push_back(-9999.);
	   (*electron_track_eta_).push_back(-9999.);
	   (*electron_track_phi_).push_back(-9999.);
	 }

       }//end fo requiring not a duplicate (this is for FAST SIM)
     }//end of requiring Electron Id
     
     //     float trackIso      = ele.dr03TkSumPt();
     //     float ecalIso       = ele.dr03EcalRecHitSumEt();
     //     float hcalIso       = ele.dr03HcalTowerSumEt();


   }//end of loop over electrons



   //get the builder:
   edm::ESHandle<TransientTrackBuilder> theBuilder;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theBuilder);
   //   cout <<" Got a "<<typeid(*theBuilder).name()<<endl;
   //   cout << "Field at origin (in Testla): "<< (*theBuilder).field()->inTesla(GlobalPoint(0.,0.,0.))<<endl;


   ////////////////////////////////////////////////////////////////////////////////
   //Store track information 
   ////////////////////////////////////////////////////////////////////////////////
   
   edm::Handle<reco::TrackCollection> trackCollection;
   iEvent.getByLabel("generalTracks", trackCollection);
   const reco::TrackCollection *tracks = trackCollection.product();
   reco::TrackCollection::const_iterator tciter;

   if ( tracks->size() > 0 )
     {
       // Loop on tracks
       int countTracks=0;
       //       cout<<"RIGHT BEFORE TRACKS"<<endl;
       
       for ( tciter=tracks->begin(); tciter!=tracks->end(); ++tciter)
	 {
	   

	   const reco::Track tr = (*tciter);
	   
	   reco::TransientTrack transientTrack = theBuilder->build(tr);
	   transientTrack.setBeamSpot(beamSpot);
	   TrajectoryStateClosestToBeamLine traj = transientTrack.stateAtBeamLine();
	   Measurement1D meas = traj.transverseImpactParameter();
	   //       Measurement1D measSign = traj.transverseImpactParameter().significance();
	   //double d0 = meas.value();
	   //double d0_error = meas.error();
	   //	   d0_vector.push_back(d0);
	   //d0Err_vector.push_back(d0_error);

	   
	   //	   dxy_vector.push_back(tciter->dxy(beamSpot.position()));
	   //dxyErr_vector.push_back(tciter->dxyError());
	   //normalizedChi2_vector.push_back(tciter->normalizedChi2());
	   //pxlhits_vector.push_back(tciter->hitPattern().pixelLayersWithMeasurement());
	   //silcnhits_vector.push_back(tciter->hitPattern().trackerLayersWithMeasurement());
	   
	   //	   if(tciter->pt()>10000){
	   //	   cout<<"vertex_track pt: "<<tciter->pt()<<", trackNum: "<<countTracks<<endl;
	   //	   cout<<"d0, error, signif: "<<d0<<", "<<d0_error<<", "<<meas.significance()<<endl;
	   //	   double dxy_sign=-999;
	   //	   if(tciter->dxyError()!=0){dxy_sign = tciter->dxy(beamSpot.position())/tciter->dxyError();}
	   //	   cout<<"dxy beamspot, error, sig: "<<tciter->dxy(beamSpot.position())<<", "<<tciter->dxyError()<<", "<<dxy_sign<<endl;
	     //	   }
	   (*tracks_chi2_).push_back(tciter->chi2());
	   (*tracks_ndof_).push_back(tciter->ndof());
	   (*tracks_charge_).push_back(tciter->charge());
	   (*tracks_pt_).push_back(tciter->pt());
	   (*tracks_px_).push_back(tciter->px());
	   (*tracks_py_).push_back(tciter->py());
	   (*tracks_pz_).push_back(tciter->pz());
	   (*tracks_eta_).push_back(tciter->eta());
	   (*tracks_phi_).push_back(tciter->phi());
	   (*tracks_d0_).push_back(tciter->dxy(beamSpot.position()));
	   (*tracks_d0Err_).push_back(tciter->dxyError());
	   (*tracks_pxlhits_).push_back(tciter->hitPattern().pixelLayersWithMeasurement());
	   (*tracks_silcnhits_).push_back(tciter->hitPattern().trackerLayersWithMeasurement());
	   (*tracks_dz_).push_back(tciter->dz());
	   (*tracks_vx_).push_back(tciter->vx());
	   (*tracks_vy_).push_back(tciter->vy());
	   (*tracks_vz_).push_back(tciter->vz());
	   countTracks++;
	 }

     }
   
   //   cout<<"RIGHT AFTER TRACKS"<<endl;
   //   cout<<"Get right before looping over vertices:"<<endl;
   
   ////////////////////////////////////////////////////////////////////////////////
   //Store vertex information 
   ////////////////////////////////////////////////////////////////////////////////
   for(vertex_i = vertexs->begin(); vertex_i<vertexs->end(); vertex_i++){
     (*vertex_ntracks_).push_back(vertex_i->tracksSize());
     (*vertex_x_).push_back(vertex_i->x());
     (*vertex_y_).push_back(vertex_i->y());
     (*vertex_z_).push_back(vertex_i->z());
     (*vertexerr_x_).push_back(vertex_i->xError());
     (*vertexerr_y_).push_back(vertex_i->yError());
     (*vertexerr_z_).push_back(vertex_i->zError());
     (*vertex_chi2_).push_back(vertex_i->chi2());
     (*vertex_ndof_).push_back(vertex_i->ndof());
     (*vertex_isFake_).push_back(vertex_i->isFake());
     (*vertex_isValid_).push_back(vertex_i->isValid());

     //     GlobalPoint vert(vertex_i->x(), vertex_i->y(), vertex_i->z());

     std::vector<float> pt_vector;
     std::vector<float> eta_vector;
     std::vector<float> phi_vector;
     std::vector<float> d0_vector;
     std::vector<float> d0Err_vector;
     std::vector<float> dxy_vector;
     std::vector<float> dxyErr_vector;
     std::vector<float> normalizedChi2_vector;
     std::vector<float> pxlhits_vector;
     std::vector<float> silcnhits_vector;
     std::vector<float> dz_vector;
     std::vector<float> dzErr_vector;
     std::vector<float> ptErr_vector;
     std::vector<int> high_purity_vector;
     std::vector<float> weight_vector;

     int trackNum = 0;
     //     cout<<"Get right before looping over tracks at vertex:"<<endl;

     //Should I store all tracks for each vertex? So vector of tracks for each vertex?
     for(reco::Vertex::trackRef_iterator vertex_Tracks = vertex_i->tracks_begin();vertex_Tracks<vertex_i->tracks_end(); vertex_Tracks++){

       //       const reco::Track tr = (*vertex_Tracks[trackNum]);
       //       const reco::Track tr = ();

       reco::TransientTrack transientTrack =  theBuilder->build(vertex_Tracks->castTo<reco::TrackRef>());

       //       const reco::TransientTrack transientTrack = theBuilder->build(tr,beamSpot);
       transientTrack.setBeamSpot(beamSpot);
       TrajectoryStateClosestToBeamLine traj = transientTrack.stateAtBeamLine();
       Measurement1D meas = traj.transverseImpactParameter();
       //       Measurement1D measSign = traj.transverseImpactParameter().significance();
       double d0 = meas.value();
       double d0_error = meas.error();
       d0_vector.push_back(d0);
       d0Err_vector.push_back(d0_error);




       dxy_vector.push_back((*vertex_Tracks)->dxy(beamSpot.position()));
       //dxy_vector.push_back((*vertex_Tracks)->dxy(vert));
       dxyErr_vector.push_back((*vertex_Tracks)->dxyError());
       normalizedChi2_vector.push_back((*vertex_Tracks)->normalizedChi2());
       pxlhits_vector.push_back((*vertex_Tracks)->hitPattern().pixelLayersWithMeasurement());
       silcnhits_vector.push_back((*vertex_Tracks)->hitPattern().trackerLayersWithMeasurement());

       pt_vector.push_back((*vertex_Tracks)->pt());

       //       cout<<"Is stateAtBeamLine valid: "<<transientTrack.stateAtBeamLine().isValid()<<endl;
       //       if(vertex_Tracks[0]->pt()>10000){
       //	 cout<<"vertex_track pt: "<<vertex_Tracks[0]->pt()<<", trackNum: "<<trackNum<<endl;
       //	 cout<<"d0, error, signif: "<<d0<<", "<<d0_error<<", "<<meas.significance()<<endl;
       // cout<<"dxy, error: "<<(*vertex_Tracks)->dxy(beamSpot.position())<<", "<<(*vertex_Tracks)->dxyError()<<endl;;
       // cout<<"Pt, eta, phi: "<<(*vertex_Tracks)->pt()<<", "<<(*vertex_Tracks)->eta()<<", "<<(*vertex_Tracks)->phi()<<endl;;
       // cout<<"Normalized chi2: "<<(*vertex_Tracks)->normalizedChi2()<<endl;
       // cout<<"Pxl hits: "<<(*vertex_Tracks)->hitPattern().pixelLayersWithMeasurement()<<endl;
       // cout<<"Silicon hits: "<<(*vertex_Tracks)->hitPattern().trackerLayersWithMeasurement()<<endl;
       eta_vector.push_back((*vertex_Tracks)->eta());
       phi_vector.push_back((*vertex_Tracks)->phi());
       //       }
       dz_vector.push_back((*vertex_Tracks)->dz());
       dzErr_vector.push_back((*vertex_Tracks)->dzError());
       ptErr_vector.push_back((*vertex_Tracks)->ptError());
       high_purity_vector.push_back((*vertex_Tracks)->quality(reco::TrackBase::highPurity));
       weight_vector.push_back(vertex_i->trackWeight(*vertex_Tracks));
       
       //       std::cout<<"The track quality is: "<<(*vertex_Tracks)->quality(reco::TrackBase::highPurity)<<endl;
       trackNum++;
     }

     (*vertex_tracks_pt_).push_back(pt_vector);
     (*vertex_tracks_eta_).push_back(eta_vector);
     (*vertex_tracks_phi_).push_back(phi_vector);
     (*vertex_tracks_d0_).push_back(d0_vector);
     (*vertex_tracks_d0Err_).push_back(d0Err_vector);
     (*vertex_tracks_dz_).push_back(dz_vector);
     (*vertex_tracks_dzErr_).push_back(dzErr_vector);
     (*vertex_tracks_ptErr_).push_back(ptErr_vector);
     (*vertex_tracks_highPurity_).push_back(high_purity_vector);       
     (*vertex_tracks_weight_).push_back(weight_vector);       
     (*vertex_tracks_dxy_).push_back(dxy_vector);
     (*vertex_tracks_dxyErr_).push_back(dxyErr_vector);
     (*vertex_tracks_normChi2_).push_back(normalizedChi2_vector);
     (*vertex_tracks_pxlhits_).push_back(pxlhits_vector);
     (*vertex_tracks_silcnhits_).push_back(silcnhits_vector);
     pt_vector.clear();
     eta_vector.clear();
     phi_vector.clear();
     d0_vector.clear();
     d0Err_vector.clear();
     dz_vector.clear();
     dzErr_vector.clear();
     ptErr_vector.clear();
     high_purity_vector.clear();
     weight_vector.clear();
     dxy_vector.clear();
     dxyErr_vector.clear();
     normalizedChi2_vector.clear();
     pxlhits_vector.clear();
     silcnhits_vector.clear();

   }

   //   cout<<"Get right after looping over vertices:"<<endl;

   ////////////////////////////////////////////////////////////////////////////////
   //Store vertex information for e's that pass e e criteria
   ////////////////////////////////////////////////////////////////////////////////
   int count_Nvertices_match_ee = 0;
   for(vertex_i = vertexs->begin(); vertex_i<vertexs->end(); vertex_i++){

     double vertex_ntracks = vertex_i->tracksSize();
     double vertex_z = vertex_i->z();
     if(fabs(vertex_z)<24){

       for(int l = 0;l<n_ElsPassingCuts;l++){
	 for(int k = l+1;k<n_ElsPassingCuts;k++){
	   
	   if((electronCharge[l]*electronCharge[k])<0){
	     if(vertex_ntracks<=100){
	       
	       bool pass_electron_assoc = false;
	       bool pass_electron2_assoc = false;
	       for(reco::Vertex::trackRef_iterator vertex_Tracks = vertex_i->tracks_begin();vertex_Tracks<vertex_i->tracks_end(); vertex_Tracks++){
		 if( fabs((*vertex_Tracks)->pt()-electronTrackPt[l])<0.001 && 
		     fabs((*vertex_Tracks)->eta()-electronTrackEta[l])<0.001 && 
		     fabs((*vertex_Tracks)->phi()-electronTrackPhi[l])<0.001){
		   pass_electron_assoc = true;
		 }
		 
		 if( fabs((*vertex_Tracks)->pt()-electronTrackPt[k])<0.001 && 
		     fabs((*vertex_Tracks)->eta()-electronTrackEta[k])<0.001 && 
		     fabs((*vertex_Tracks)->phi()-electronTrackPhi[k])<0.001){
		   pass_electron2_assoc = true;
		 }
	       }//Loop over tracks to see if the leptons are there
	       
	       if(pass_electron_assoc&&pass_electron2_assoc){
		 count_Nvertices_match_ee++;		 
		 nextra_tracks_ee = vertex_ntracks-2;
		 (*vertex_nextra_tracks_ee_).push_back(nextra_tracks_ee);
		 (*vertexCandE1Index_ee_).push_back(l);
		 (*vertexCandE2Index_ee_).push_back(k);
		 //		 cout<<"This event passes mu+ mu- criteria and should be vetoed"<<endl;
	       }
	       
	     }//end of vertex requirement
	   }//end of if statement for charge
	 }
       }
     }}// end of loop over vertices
   if(count_Nvertices_match_ee>1){
     cout<<"More than one matching e e pair"<<endl;
   }
   



   ////////////////////////////////////////////////////////////////////////////////
   //Veto events that pass mu mu criteria
   ////////////////////////////////////////////////////////////////////////////////
   bool veto_event = false;
   int count_Nvertices_match_mumu = 0;
   for(vertex_i = vertexs->begin(); vertex_i<vertexs->end(); vertex_i++){
     double vertex_ntracks = vertex_i->tracksSize();
     double vertex_z = vertex_i->z();
     if(fabs(vertex_z)<24){

       for(int l = 0;l<n_MuonsPassingCuts;l++){
	 for(int k = l+1;k<n_MuonsPassingCuts;k++){
	   
	   if((muonCharge[l]*muonCharge[k])<0){
	     //	     if(vertex_ntracks<=17){
	     //	     cout<<"vertex extra tracks "<<vertex_ntracks-2<<endl;
	     if(vertex_ntracks<=100){
	       
	       bool pass_muon_assoc = false;
	       bool pass_muon2_assoc = false;
	       for(reco::Vertex::trackRef_iterator vertex_Tracks = vertex_i->tracks_begin();vertex_Tracks<vertex_i->tracks_end(); vertex_Tracks++){
		 if( fabs((*vertex_Tracks)->pt()-muonTrackPt[l])<0.001 && 
		     fabs((*vertex_Tracks)->eta()-muonTrackEta[l])<0.001 && 
		     fabs((*vertex_Tracks)->phi()-muonTrackPhi[l])<0.001){
		   pass_muon_assoc = true;
		 }
		 
		 if( fabs((*vertex_Tracks)->pt()-muonTrackPt[k])<0.001 && 
		     fabs((*vertex_Tracks)->eta()-muonTrackEta[k])<0.001 && 
		     fabs((*vertex_Tracks)->phi()-muonTrackPhi[k])<0.001){
		   pass_muon2_assoc = true;
		 }
	       }//Loop over tracks to see if the leptons are there
	       
	       if(pass_muon_assoc&&pass_muon2_assoc){

		 /*
		 if(nextra_tracks_mumu<15){
		   for(reco::Vertex::trackRef_iterator vertex_Tracks = vertex_i->tracks_begin();vertex_Tracks<vertex_i->tracks_end(); vertex_Tracks++){
		     (*vertex_Tracks)->pt();
		     (*vertex_Tracks)->eta();
		     (*vertex_Tracks)->eta();
		     cout<<"Pt, eta, phi: "<<(*vertex_Tracks)->pt()<<", "<<(*vertex_Tracks)->eta()<<", "<<(*vertex_Tracks)->phi()<<endl;;
		   }
		 }
		 */

		 count_Nvertices_match_mumu++;
		 veto_event=true;		 
		 nextra_tracks_mumu = vertex_ntracks-2;
		 //		 cout<<"pass vertexing, extra tracks "<<vertex_ntracks-2<<endl;
		 (*vertex_nextra_tracks_mumu_).push_back(nextra_tracks_mumu);
		 (*vertexCandMu1Index_mumu_).push_back(l);
		 (*vertexCandMu2Index_mumu_).push_back(k);
		 //		 cout<<"This event passes mu+ mu- criteria and should be vetoed"<<endl;
	       }
	       
	     }//end of vertex requirement
	   }//end of if statement for charge
	 }
       }
     }}// end of loop over vertices
   if(count_Nvertices_match_mumu>1){
     cout<<"More than one matching mu mu pair"<<endl;
   }


   //////////////////////////////////////////////////////////////////////////////////////////////
   //Loop to see how many extra tracks are on the emu vertex (and if they are on the same vertex)
   //////////////////////////////////////////////////////////////////////////////////////////////
   int count_Nvertices_match_emus = 0;
   for(vertex_i = vertexs->begin(); vertex_i<vertexs->end(); vertex_i++){
     double vertex_ntracks = vertex_i->tracksSize();
     double vertex_z = vertex_i->z();
     if(fabs(vertex_z)<24){

     for(int k = 0;k<n_ElsPassingCuts;k++){
       
       for(int l = 0;l<n_MuonsPassingCuts;l++){
	 
	 //	 if(n_MuonsPassingCuts==1&&n_ElsPassingCuts==1&&(muonCharge*electronCharge)<0){
	 if((muonCharge[l]*electronCharge[k])<0){
	   
	   //Loop over tracks from vertex and look for matches

	   if(vertex_ntracks<=100){
	     bool pass_muon_assoc = false;
	     bool pass_electron_assoc = false;
	     double pt_track=0;
	     for(reco::Vertex::trackRef_iterator vertex_Tracks = vertex_i->tracks_begin();vertex_Tracks<vertex_i->tracks_end(); vertex_Tracks++){

	       if( fabs((*vertex_Tracks)->pt()-muonTrackPt[l])<0.001 && 
		   fabs((*vertex_Tracks)->eta()-muonTrackEta[l])<0.001 && 
		   fabs((*vertex_Tracks)->phi()-muonTrackPhi[l])<0.001){
		 pass_muon_assoc = true;
	       }
	       else if (fabs((*vertex_Tracks)->pt()-electronTrackPt[k])<0.001 && 
		   fabs((*vertex_Tracks)->eta()-electronTrackEta[k])<0.001 && 
		   fabs((*vertex_Tracks)->phi()-electronTrackPhi[k])<0.001){
		 pass_electron_assoc = true;
	       }
	       else{pt_track = (*vertex_Tracks)->pt();
	       }

	     }//Loop over tracks to see if the leptons are there
	     
	     if(pass_muon_assoc&&pass_electron_assoc&&!veto_event){
	       count_Nvertices_match_emus++;
	       nextra_tracks = vertex_ntracks-2;
	       if(nextra_tracks==1){h_pt_extra_track->Fill(pt_track);}
	       (*vertex_nextra_tracks_).push_back(nextra_tracks);
	       (*vertexCandMuIndex_).push_back(l);
	       (*vertexCandEIndex_).push_back(k);
	     }

	   }//require less than 16 tracks matched to the vertex
	 
	 }//end of loop over requiring at least one electron and muon in the event
       }}//end of for loops over leptons
     
     }//end of requiring vertex have |z|<24 cm

   }//end of loop over vertices
   if(count_Nvertices_match_emus>1){
     cout<<"More than one matching e mu pair"<<endl;
   }







   ////////////////////////////////////////////////////////////////////////////////
   //Reco to Sim vertex association
   ////////////////////////////////////////////////////////////////////////////////
   bool isRECO = false;
   if(isRECO){

     //   edm::Handle<reco::TrackCollection> trackCollection;
     //   iEvent.getByLabel("generalTracks", trackCollection);
     //   const reco::TrackCollection *tracks = trackCollection.product();
     //   reco::TrackCollection::const_iterator tciter;
     
     edm::Handle<edm::View<reco::Track> > trackCollectionH;
     iEvent.getByLabel("generalTracks",trackCollectionH);
     //   const reco::TrackCollection *tracks = trackCollectionH.product();

     //open a collection of Tracking Particles
     edm::Handle<TrackingParticleCollection> TPtracks;
     iEvent.getByLabel("mergedtruth","MergedTrackTruth",TPtracks);
     
     //get the associator
     iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theAssociator);
     
     //associate
     //TrackAssociatorBase* assoc = ((TrackAssociatorBase*)theAssociator.product());
     reco::RecoToSimCollection recSimColl = theAssociator->associateRecoToSim(trackCollectionH,TPtracks, &iEvent);
     //   cout<<"I have found: "<<recSimColl.size()<<" associations in total."<<endl;
     
     //   int iter = 0;
     //   reco::TrackCollection::size_type iter=0;
     edm::View<reco::Track>::size_type iter=0;
     //look at each track associated to a particular sim info
     //begining of loop over map
     for (reco::RecoToSimCollection::const_iterator RtSit=recSimColl.begin();RtSit!=recSimColl.end();++RtSit){
     const std::vector<std::pair<TrackingParticleRef,double> > & tp = RtSit->val;
     
     //     reco::TrackRef tracks(trackCollectionH, iter);
     RefToBase<reco::Track> tracks(trackCollectionH, iter);
     //need to open collection of tracks before doing this
     //     double muon_pt_reco_assoc = (*tracks)[iter].pt();
     //double muon_pt_reco_assoc = tracks->pt();
     //double muon_eta_reco_assoc = tracks->eta();
     //double muon_phi_reco_assoc = tracks->phi();
     
     (*tracks_pt_).push_back(tracks->pt());
     (*tracks_eta_).push_back(tracks->eta());
     (*tracks_phi_).push_back(tracks->phi());

     //cout<<"Reco track with pt, eta, phi: "<<muon_pt_reco_assoc<<", "<<muon_eta_reco_assoc<<", "<<muon_phi_reco_assoc<<endl;
     
     //     cout<<"I have found: "<<tp.size()<<" tracking particle associated with this reco::Track."<<endl;
     if(tp.size()!=0)
       {
	 //take the match with best quality
	 std::vector<std::pair<TrackingParticleRef,double> >::const_iterator vector_iterator = tp.begin();
	 const std::pair<TrackingParticleRef,double> & pair_in_vector = *vector_iterator;
	 const TrackingParticleRef & trp = pair_in_vector.first;


	 //const double & quality = pair_in_vector.second;
	 //int particle_ID = trp->pdgId();
	 //	 cout<<"The quality of the closest for tight association is:"<<quality<<"The particle ID is: "<<particle_ID<<endl;
	 //	 quality_tight->Fill(quality);
	 //	 h_ID_vs_quality_tight->Fill(particle_ID,quality);
	 //cout<<"Particle id is: "<<particle_ID<<endl;

	 
	 //double sim_pt_reco_assoc = trp->pt();
	 //double sim_eta_reco_assoc = trp->eta();
	 //double sim_phi_reco_assoc = trp->phi();
	 //double sim_id = trp->pdgId();
	 //	 cout<<"Reco track with pt, eta, phi: "<<sim_pt_reco_assoc<<", "<<sim_eta_reco_assoc<<", "<<sim_phi_reco_assoc<<endl;
	 (*tracks_sim_pt_).push_back(trp->pt());
	 (*tracks_sim_eta_).push_back(trp->eta());
	 (*tracks_sim_phi_).push_back(trp->phi());
	 (*tracks_sim_id_).push_back(trp->pdgId());


	 int parent_ID = 0;
	 const TrackingParticle::TrackingVertexRef tvr = trp->parentVertex();
	 const TrackingParticleRefVector tvfv = tvr->sourceTracks();
	 

	 //cout<<"The size of the TrackingParticleRefVector is: "<<tvfv.size()<<endl;
	 
	 if(tvfv.size()==0)
	   {//cout<<"There is no parent particle for this particle"<<endl;
	   }
	 else
	   {
	     TrackingVertex::tp_iterator TVit = tvr->sourceTracks_begin();
	     const TrackingParticleRef & trp_parent = *TVit;
	     parent_ID = trp_parent->pdgId();
	     //cout<<"The parent ID number is: "<<parent_ID<<endl;
	   }
	 (*tracks_sim_motherid_).push_back(parent_ID);	 
	 //return parent_ID;
	 
	 
	 int motherID=0;
	 int num_gen_part=0;
	 int ipdg = 0;
	 for(TrackingParticle::genp_iterator tpgi = trp->genParticle_begin();tpgi!=trp->genParticle_end();tpgi++)
	   {
	     
	     const HepMC::GenParticle *part  = &(*(*tpgi));
	     
	     //	     cout<<"The part pointer for GenParticle is: "<<part<<endl;
	     ipdg = part->pdg_id();
	     //cout<<"The Id for the GenParticle associated with track is: "<<ipdg<<endl;

	     if(!part->production_vertex())continue;
	     if(part->production_vertex()->particles_in_size()==0)continue;
	     const HepMC::GenParticle *mother = *(part->production_vertex()->particles_in_const_begin());
	     //	     cout<<"the mother Id is: "<<mother->pdg_id()<<endl;
	     motherID = mother->pdg_id();
	     num_gen_part++;
	     
	   }

	 (*tracks_gen_id_).push_back(ipdg);	 
	 (*tracks_gen_motherid_).push_back(motherID);

       }//end of if statement making sure there is at least one associated sim particle.     
     
     iter++;
     
   }//end of looping over map of recoToSim collection

}

   ////////////////////////////////////////////////////////////////////////////////
   //End of Reco to Sim vertex association
   ////////////////////////////////////////////////////////////////////////////////

   /*
   //   reco::SimToRecoCollection s2rTracks = associatorByHits->associateSimToReco (trackCollectionH,TPCollectionH,&event );
   reco::RecoToSimCollection r2sTracks = theAssociator->associateRecoToSim(trackCollectionH,TPtracks,&iEvent );
   reco::SimToRecoCollection s2rTracks = theAssociator->associateSimToReco(trackCollectionH,TPtracks,&iEvent );

   //Sim vertices
   edm::Handle<TrackingVertexCollection>  TVCollection;
   iEvent.getByLabel("mergedtruth","MergedTrackTruth",TVCollection);
   const TrackingVertexCollection tVC   = *(TVCollection.product());
 
   //Reco Vertex Collection
   edm::Handle<edm::View<reco::Vertex> > vertexCollection;
   iEvent.getByLabel("offlinePrimaryVertices", vertexCollection);
   const edm::View<reco::Vertex> vC = *(vertexCollection.product());

   edm::ESHandle<VertexAssociatorBase> theTracksAssociator;
   iSetup.get<VertexAssociatorRecord>().get("VertexAssociatorByTracks",theTracksAssociator);
   VertexAssociatorBase* associatorByTracks;
   associatorByTracks = (VertexAssociatorBase *) theTracksAssociator.product();

   //reco::VertexRecoToSimCollection r2sVertices = associatorByTracks->associateRecoToSim(vertexCollection,TVCollection,iEvent,r2sTracks);
   //   reco::VertexSimToRecoCollection s2rVertices = associatorByTracks->associateSimToReco(vertexCollection,TVCollection,iEvent,s2rTracks);
   

   //Looking at sim vertices
   int iTV = 0;
   cout<<"Num of sim vertices: "<<tVC.size()<<endl;
   for (TrackingVertexCollection::const_iterator tV = tVC.begin();tV != tVC.end(); ++tV,++iTV) {
     TrackingVertexRef tVertexR = TrackingVertexRef(TVCollection,iTV);
     double nSimTracks = (tV->daughterTracks()).size();
       // Loop over daughter tracks of TrackingVertex
     int numMu= 0;
     for (TrackingVertex::tp_iterator simDaughter = tV->daughterTracks_begin();simDaughter != tV->daughterTracks_end(); ++simDaughter) {
       TrackingParticleRef tp = *simDaughter;
       // cout<<"Particle id at this vertex is: "<<tp->pdgId()<<endl;
       if(fabs(tp->pdgId())==13&&tp->pt()>5){numMu++;}
     }
     if(numMu>1){
       cout<<"Vertex num is: "<<iTV<<endl;
       cout<<"Size of this vertex is: "<<nSimTracks<<endl;
       
       for (TrackingVertex::tp_iterator simDaughter = tV->daughterTracks_begin();simDaughter != tV->daughterTracks_end(); ++simDaughter) {
	 TrackingParticleRef tp = *simDaughter;
	 cout<<"Particle id,pt at this vertex is: "<<tp->pdgId()<<",pt: "<<tp->pt()<<", eta: "<<tp->eta()<<endl;
       }
     }//end of requiring more than one muon

   }
   */

   /*
   
   int sim_vertex_count = 0;
   for (reco::VertexSimToRecoCollection::const_iterator iS2R = s2rVertices.begin();iS2R != s2rVertices.end(); ++iS2R) {
     cout<<"Sim vertex count: "<<sim_vertex_count<<endl;
     TrackingVertexRef simVertex = (iS2R -> key);
     //     HepLorentzVector simVec = simVertex->position();
     //     math::XYZPoint   simPos = math::XYZPoint(simVec.x(),simVec.y(),simVec.z());
     double ntrue = simVertex->daughterTracks().size();
     //     std::vector<std::pair<reco::VertexRef, double> > recoVertices = iS2R->val;
     int nummus=0;
     for(TrackingVertex::tp_iterator TVit = simVertex->sourceTracks_begin();TVit!=simVertex->sourceTracks_end();TVit++){

       const TrackingParticleRef & trp = *TVit;
       //cout<<"Sim particle id: "<<trp->pdgId()<<endl;
       if(fabs(trp->pdgId())==13&&trp->pt()>10){
	 nummus++;
       }
       
     }
     if(nummus>1){
       cout<<"Num particles at simvertex "<<simVertex->daughterTracks().size()<<endl;
       for(TrackingVertex::tp_iterator TVit = simVertex->sourceTracks_begin();TVit!=simVertex->sourceTracks_end();TVit++){
	 const TrackingParticleRef & trp = *TVit;
	 cout<<"Sim particle id: "<<trp->pdgId()<<endl;
       }
     }


          for (std::vector<std::pair<VertexRef, double> >::const_iterator iMatch = recoVertices.begin();iMatch != recoVertices.end(); ++iMatch) {
       VertexRef recoV = iMatch->first;
       double qual  = iMatch->second;
       math::XYZPoint recoPos = (iMatch -> first) -> position();
       double nreco = (iMatch->first)->tracksSize();
       
       double xmiss = simPos.X() - recoPos.X();
       double ymiss = simPos.Y() - recoPos.Y();
       double zmiss = simPos.Z() - recoPos.Z();
       double rmiss = sqrt(xmiss*xmiss+ymiss*ymiss+zmiss*zmiss);
       
       sr_xMiss->Fill(simPos.X() - recoPos.X());
       sr_yMiss->Fill(simPos.Y() - recoPos.Y());
       sr_zMiss->Fill(simPos.Z() - recoPos.Z());
       sr_rMiss->Fill(rmiss);

       sr_zVert->Fill(simPos.Z());
       sr_zTrue->Fill(recoPos.Z());

       sr_nTrue->Fill(ntrue);
       sr_nReco->Fill(nreco);
       sr_qual->Fill(qual);
     }
     
   }//end of loop over sim vertices
   */
   tree_->Fill();
   (*muon_id_).clear();
   (*muon_pt_).clear();
   (*muon_px_).clear();
   (*muon_py_).clear();
   (*muon_pz_).clear();
   (*muon_et_).clear();
   (*muon_energy_).clear();
   (*muon_charge_).clear();
   (*muon_eta_).clear();
   (*muon_phi_).clear();
   (*muon_track_pt_).clear();
   (*muon_track_eta_).clear();
   (*muon_track_phi_).clear();

   (*gen_pt_).clear();
   (*gen_px_).clear();
   (*gen_py_).clear();
   (*gen_pz_).clear();
   (*gen_et_).clear();
   (*gen_energy_).clear();
   (*gen_charge_).clear();
   (*gen_eta_).clear();
   (*gen_phi_).clear();
   (*gen_pdgId_).clear();
   (*electron_id_).clear();
   (*electron_pt_).clear();
   (*electron_px_).clear();
   (*electron_py_).clear();
   (*electron_pz_).clear();
   (*electron_et_).clear();
   (*electron_energy_).clear();
   (*electron_charge_).clear();
   (*electron_eta_).clear();
   (*electron_phi_).clear();
   (*electron_track_pt_).clear();
   (*electron_track_eta_).clear();
   (*electron_track_phi_).clear();

   (*trigger_prescaleValue_).clear();
   (*trigger_name_).clear();
   (*trigger_decision_).clear();
   (*vertex_nextra_tracks_).clear();
   (*vertexCandMuIndex_).clear();
   (*vertexCandEIndex_).clear();

   (*vertex_nextra_tracks_ee_).clear();
   (*vertexCandE1Index_ee_).clear();
   (*vertexCandE2Index_ee_).clear();

   (*vertex_nextra_tracks_mumu_).clear();
   (*vertexCandMu1Index_mumu_).clear();
   (*vertexCandMu2Index_mumu_).clear();

   (*met_et_).clear();
   (*met_px_).clear();
   (*met_py_).clear();
   (*met_sumEt_).clear();
   
   (*vertex_ntracks_).clear();
   (*vertex_x_).clear();
   (*vertex_y_).clear();
   (*vertex_z_).clear();
   (*vertexerr_x_).clear();
   (*vertexerr_y_).clear();
   (*vertexerr_z_).clear();
   (*vertex_chi2_).clear();
   (*vertex_ndof_).clear();
   (*vertex_isFake_).clear();
   (*vertex_isValid_).clear();
   (*vertex_tracks_pt_).clear();
   (*vertex_tracks_eta_).clear();
   (*vertex_tracks_phi_).clear();
   (*vertex_tracks_highPurity_).clear();
   (*vertex_tracks_d0_).clear();
   (*vertex_tracks_d0Err_).clear();
   (*vertex_tracks_dz_).clear();
   (*vertex_tracks_dzErr_).clear();
   (*vertex_tracks_ptErr_).clear();
   (*vertex_tracks_weight_).clear();
   (*vertex_tracks_dxy_).clear();
   (*vertex_tracks_dxyErr_).clear();
   (*vertex_tracks_normChi2_).clear();
   (*vertex_tracks_pxlhits_).clear();
   (*vertex_tracks_silcnhits_).clear();

   (*tracks_chi2_).clear();
   (*tracks_ndof_).clear();
   (*tracks_charge_).clear();
   (*tracks_pt_).clear();
   (*tracks_px_).clear();
   (*tracks_py_).clear();
   (*tracks_pz_).clear();
   (*tracks_eta_).clear();
   (*tracks_phi_).clear();
   (*tracks_d0_).clear();
   (*tracks_d0Err_).clear();
   (*tracks_pxlhits_).clear();
   (*tracks_silcnhits_).clear();
   (*tracks_dz_).clear();
   (*tracks_vx_).clear();
   (*tracks_vy_).clear();
   (*tracks_vz_).clear();
   (*tracks_sim_pt_).clear();
   (*tracks_sim_eta_).clear();
   (*tracks_sim_phi_).clear();
   (*tracks_sim_id_).clear();
   (*tracks_sim_motherid_).clear();
   (*tracks_gen_id_).clear();
   (*tracks_gen_motherid_).clear();



}


// ------------ method called once each job just before starting event loop  ------------
void 
Ntupler::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Ntupler::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
Ntupler::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  //  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
  if (hltConfig_.init(iRun,iSetup,"HLT",changed)) {
    // if init returns TRUE, initialisation has succeeded!
    if (changed) {
      // The HLT config has actually changed wrt the previous Run, hence rebook your
      // histograms or do anything else dependent on the revised HLT config
     
    }
  } else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    //    coutr") << " HLT config extraction failure with process name " << processName_<<endl;
    //    LogError("MyAnalyzer") << " HLT config extraction failure with process name " << "HLT";
    std::cout << " HLT config extraction failure with process name " << "HLT"<<std::endl;
    // In this case, all access methods will return empty values!
  }

}

// ------------ method called when ending the processing of a run  ------------
void 
Ntupler::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Ntupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Ntupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Ntupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntupler);
