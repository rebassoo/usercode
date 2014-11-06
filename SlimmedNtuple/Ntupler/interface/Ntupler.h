// -*- C++ -*-
//
// Package:    Ntupler
// Class:      Ntupler
// 
/**\class Ntupler Ntupler.cc SlimmedNtuple/Ntupler/src/Ntupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Finn Rebassoo
//         Created:  Thu Sep 20 13:42:33 CDT 2012
// $Id: Ntupler.h,v 1.6 2013/04/23 22:15:42 rebassoo Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TTree.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/METReco/interface/PFMET.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

//For Reco to sim association
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"

#include "SimTracker/Records/interface/VertexAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimTracker/VertexAssociation/interface/VertexAssociatorBase.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include <FWCore/Framework/interface/ESHandle.h> 


//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
//#include "HepPDT/ParticleDataTable.hh"

typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;

//
// class declaration
//

class Ntupler : public edm::EDAnalyzer {
   public:
      explicit Ntupler(const edm::ParameterSet&);
      ~Ntupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  TTree * tree_;

  HLTConfigProvider hltConfig_;

  bool isMC;
  std::vector<edm::InputTag> isoValInputTags;
  
  TH1F * h_trueNumInteractions;
  TH1F * h_trueNumInteractions0;
  TH1F * h_pt_extra_track;
  edm::LumiReWeighting *LumiWeights;

  //For Reco to sim association
  edm::ESHandle<TrackAssociatorBase> theAssociator;

  uint * ev_;
  uint * run_;
  uint * lumiblock_;
  uint * experimentType_;
  uint * bunchCrossing_;
  uint * orbitNumber_;
  float * Tnpv_;
  float * pileupWeight_;

  std::vector<bool> * muon_id_;
  std::vector<float> * muon_pt_;
  std::vector<float> * muon_px_;
  std::vector<float> * muon_py_;
  std::vector<float> * muon_pz_;
  std::vector<float> * muon_et_;
  std::vector<float> * muon_energy_;
  std::vector<float> * muon_charge_;
  std::vector<float> * muon_eta_;
  std::vector<float> * muon_phi_;
  std::vector<float> * muon_track_pt_;
  std::vector<float> * muon_track_eta_;
  std::vector<float> * muon_track_phi_;

  std::vector<float> * gen_pt_;
  std::vector<float> * gen_px_;
  std::vector<float> * gen_py_;
  std::vector<float> * gen_pz_;
  std::vector<float> * gen_et_;
  std::vector<float> * gen_energy_;
  std::vector<float> * gen_charge_;
  std::vector<float> * gen_eta_;
  std::vector<float> * gen_phi_;
  std::vector<int>   * gen_pdgId_;

  std::vector<bool> * electron_id_;
  std::vector<float> * electron_pt_;
  std::vector<float> * electron_px_;
  std::vector<float> * electron_py_;
  std::vector<float> * electron_pz_;
  std::vector<float> * electron_et_;
  std::vector<float> * electron_energy_;
  std::vector<float> * electron_charge_;
  std::vector<float> * electron_eta_;
  std::vector<float> * electron_phi_;
  std::vector<float> * electron_track_pt_;
  std::vector<float> * electron_track_eta_;
  std::vector<float> * electron_track_phi_;
  
  std::vector<int> * trigger_prescaleValue_;
  std::vector<std::string> * trigger_name_;
  std::vector<int> * trigger_decision_;
  std::vector<int> * vertex_nextra_tracks_;
  std::vector<int> * vertexCandMuIndex_;
  std::vector<int> * vertexCandEIndex_;

  std::vector<int> * vertex_nextra_tracks_ee_;
  std::vector<int> * vertexCandE1Index_ee_;
  std::vector<int> * vertexCandE2Index_ee_;

  std::vector<int> * vertex_nextra_tracks_mumu_;
  std::vector<int> * vertexCandMu1Index_mumu_;
  std::vector<int> * vertexCandMu2Index_mumu_;

  std::vector<float> * met_et_;
  std::vector<float> * met_px_;
  std::vector<float> * met_py_;
  std::vector<float> * met_sumEt_;

  std::vector<int> * vertex_ntracks_;
  std::vector<float> * vertex_x_;
  std::vector<float> * vertex_y_;
  std::vector<float> * vertex_z_;
  std::vector<float> * vertexerr_x_;
  std::vector<float> * vertexerr_y_;
  std::vector<float> * vertexerr_z_;
  std::vector<float> * vertex_chi2_;
  std::vector<float> * vertex_ndof_;
  std::vector<int> * vertex_isFake_;
  std::vector<int> * vertex_isValid_;
  std::vector<std::vector<float>> * vertex_tracks_pt_;
  std::vector<std::vector<float>> * vertex_tracks_eta_;
  std::vector<std::vector<float>> * vertex_tracks_phi_;
  std::vector<std::vector<float>> * vertex_tracks_d0_;
  std::vector<std::vector<float>> * vertex_tracks_d0Err_;
  std::vector<std::vector<float>> * vertex_tracks_dz_;
  std::vector<std::vector<float>> * vertex_tracks_dzErr_;
  std::vector<std::vector<float>> * vertex_tracks_ptErr_;
  std::vector<std::vector<int>> * vertex_tracks_highPurity_;
  std::vector<std::vector<float>> * vertex_tracks_weight_;
  std::vector<std::vector<float>> * vertex_tracks_dxy_;
  std::vector<std::vector<float>> * vertex_tracks_dxyErr_;
  std::vector<std::vector<float>> * vertex_tracks_normChi2_;
  std::vector<std::vector<float>> * vertex_tracks_pxlhits_;
  std::vector<std::vector<float>> * vertex_tracks_silcnhits_;

  std::vector<float> * tracks_chi2_;
  std::vector<float> * tracks_ndof_;
  std::vector<float> * tracks_charge_;
  std::vector<float> * tracks_pt_;
  std::vector<float> * tracks_px_;
  std::vector<float> * tracks_py_;
  std::vector<float> * tracks_pz_;
  std::vector<float> * tracks_eta_;
  std::vector<float> * tracks_phi_;
  std::vector<float> * tracks_d0_;
  std::vector<float> * tracks_d0Err_;
  std::vector<float> * tracks_pxlhits_;
  std::vector<float> * tracks_silcnhits_;
  std::vector<float> * tracks_dz_;
  std::vector<float> * tracks_vx_;
  std::vector<float> * tracks_vy_;
  std::vector<float> * tracks_vz_;
  std::vector<float> * tracks_sim_pt_;
  std::vector<float> * tracks_sim_eta_;
  std::vector<float> * tracks_sim_phi_;
  std::vector<float> * tracks_sim_id_;
  std::vector<float> * tracks_sim_motherid_;
  std::vector<float> * tracks_gen_id_;
  std::vector<float> * tracks_gen_motherid_;



};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
