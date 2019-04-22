// -*- C++ -*-
//
// Package:    MyAnalysis/MuonKinFit
// Class:      MuonKinFit
//
/**\class MuonKinFit MuonKinFit.cc MyAnalysis/MuonKinFit/plugins/MuonKinFit.cc

 Description: Runs a kinematic fit on the first dimuon pair and outputs muons
 with a corrected pt

 Implementation:
     Requires access to the transient track builder and thus the magnetic field.
 Only uses the first two muons in the event. Based on
 https://github.com/UFLX2MuMu/Ntupliser/blob/49e4fd57ffbdc18a98bcea64db8c736090d42eaf/DiMuons/src/MuonHelper.cc#L29
*/
//
// Original Author:  Joosep Pata
//         Created:  Mon, 22 Apr 2019 17:59:44 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParametersError.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "TLorentzVector.h"

// Class that actually runs the kinematic fit
class KinematicVertexFitter {

public:
  KinematicVertexFitter(){};
  virtual ~KinematicVertexFitter(){};

  // The reconstructed muon trajectories
  reco::TransientTrack getTransientTrack(const reco::TrackRef &trackRef) {
    reco::TransientTrack transientTrack(trackRef, paramField);
    return transientTrack;
  }

  RefCountedKinematicTree Fit(const pat::MuonCollection &pat_muons) {

    // Creating a KinematicParticleFactory
    KinematicParticleFactoryFromTransientTrack pFactory;

    // Making particles
    std::vector<RefCountedKinematicParticle> muons_for_kinematic_fit;

    // Passing transient muon tracks to the fitter
    std::vector<reco::TransientTrack> muon_tracks;

    int nmuons = 0;
    for (pat::MuonCollection::const_iterator imu = pat_muons.begin();
         imu != pat_muons.end(); ++imu) {
      if (imu->track().isAvailable()) {
        muons_for_kinematic_fit.push_back(
            pFactory.particle(getTransientTrack(imu->track()), muon_mass, chi,
                              ndf, muon_mass_sigma));
        nmuons += 1;
      }
      if (nmuons > 2) {
        break;
      }
    }

    // Fitting.
    KinematicConstrainedVertexFitter kvFitterDiMuons;
    RefCountedKinematicTree X2MuMuKinFitTree =
        kvFitterDiMuons.fit(muons_for_kinematic_fit);

    // Return a kinematic fit tree. Careful on manipulating it: see
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideKinematicVertexFit
    return X2MuMuKinFitTree;
  }

private:
  OAEParametrizedMagneticField *paramField =
      new OAEParametrizedMagneticField("3_8T");
  // The mass of a muon and the insignificant mass sigma to avoid singularities
  // in the covariance matrix.
  float const muon_mass = 0.105658367;
  float muon_mass_sigma = muon_mass * 1.e-6;
  // initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction
  // is not considered
  float chi = 0.;
  float ndf = 0.;
};

// Result of the fit is two four-momenta for the muons
class DimuonKinematicFitResult {
public:
  const TLorentzVector mu1_tlv;
  const TLorentzVector mu2_tlv;
  DimuonKinematicFitResult(const TLorentzVector &_mu1_tlv,
                           const TLorentzVector &_mu2_tlv)
      : mu1_tlv(_mu1_tlv), mu2_tlv(_mu2_tlv) {}
};

// Produces muons with fitted pt
class MuonKinFit : public edm::stream::EDProducer<> {
public:
  explicit MuonKinFit(const edm::ParameterSet &);
  ~MuonKinFit();

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event &, const edm::EventSetup &) override;
  virtual void endStream() override;

  // source muons
  edm::EDGetTokenT<edm::View<pat::Muon>> srcCands_;
};

MuonKinFit::MuonKinFit(const edm::ParameterSet &iConfig)
    : srcCands_(consumes<edm::View<pat::Muon>>(
          iConfig.getParameter<edm::InputTag>("srcCands"))) {
  produces<pat::MuonCollection>();
}

MuonKinFit::~MuonKinFit() {}

// Given two muons, runs the fit
DimuonKinematicFitResult dimuon_kinematic_fit(const pat::Muon &mu1,
                                              const pat::Muon &mu2) {

  pat::MuonCollection muonsWithTrack;
  muonsWithTrack.push_back(mu1);
  muonsWithTrack.push_back(mu2);

  TLorentzVector mu1_tlv;
  TLorentzVector mu2_tlv;

  Double_t mu1_ptErr_kinfit = 0.0;
  Double_t mu2_ptErr_kinfit = 0.0;

  RefCountedKinematicVertex dimu_vertex;

  // Instatiate KinematicVertexFitter object
  KinematicVertexFitter kinfit;
  // Fit and retrieve the tree
  RefCountedKinematicTree kinfittree = kinfit.Fit(muonsWithTrack);

  if (kinfittree->isEmpty() == 1 || kinfittree->isConsistent() == 0) {
    std::cerr << "Kinematic Fit unsuccesful" << std::endl;
  } else {
    // accessing the tree components
    kinfittree->movePointerToTheTop();
    // We are now at the top of the decay tree getting the dimuon reconstructed
    // KinematicPartlcle
    RefCountedKinematicParticle dimu_kinfit = kinfittree->currentParticle();

    // getting the dimuon decay vertex
    // RefCountedKinematicVertex
    dimu_vertex = kinfittree->currentDecayVertex();

    // Now navigating down the tree
    bool child = kinfittree->movePointerToTheFirstChild();
    // TLorentzVector mu1_tlv;

    if (child) {
      RefCountedKinematicParticle mu1_kinfit = kinfittree->currentParticle();
      AlgebraicVector7 mu1_kinfit_par =
          mu1_kinfit->currentState().kinematicParameters().vector();
      AlgebraicSymMatrix77 mu1_kinfit_cov =
          mu1_kinfit->currentState().kinematicParametersError().matrix();
      mu1_ptErr_kinfit = sqrt(mu1_kinfit_cov(3, 3) + mu1_kinfit_cov(4, 4));
      mu1_tlv.SetXYZM(mu1_kinfit_par.At(3), mu1_kinfit_par.At(4),
                      mu1_kinfit_par.At(5), mu1_kinfit_par.At(6));
      std::cout << "Mu1 chi2 = " << mu1_kinfit->chiSquared() << std::endl;
      std::cout << "Mu1 ndf = " << mu1_kinfit->degreesOfFreedom() << std::endl;
      std::cout << "Covariant matrix" << std::endl;
      std::cout << mu1_kinfit_cov(3, 3) << std::endl;
      std::cout << " - " << mu1_kinfit_cov(4, 4) << std::endl;
      std::cout << " -      -    " << mu1_kinfit_cov(5, 5) << std::endl;
      std::cout << "Muon pt uncertainty = "
                << sqrt(mu1_kinfit_cov(3, 3) + mu1_kinfit_cov(4, 4))
                << std::endl;
    }

    // Now navigating down the tree
    bool nextchild = kinfittree->movePointerToTheNextChild();

    if (nextchild) {
      RefCountedKinematicParticle mu2_kinfit = kinfittree->currentParticle();
      AlgebraicVector7 mu2_kinfit_par =
          mu2_kinfit->currentState().kinematicParameters().vector();
      AlgebraicSymMatrix77 mu2_kinfit_cov =
          mu2_kinfit->currentState().kinematicParametersError().matrix();
      mu2_ptErr_kinfit = sqrt(mu2_kinfit_cov(3, 3) + mu2_kinfit_cov(4, 4));
      mu2_tlv.SetXYZM(mu2_kinfit_par.At(3), mu2_kinfit_par.At(4),
                      mu2_kinfit_par.At(5), mu2_kinfit_par.At(6));
    }

  } // end else - isEmpty()

  std::cout << "Kin Fitted muons 1 :" << mu1_tlv.Pt() << " err "
            << mu1_ptErr_kinfit
            << "  -- Pat muons : " << muonsWithTrack.at(0).pt() << std::endl;
  std::cout << "Kin Fitted muons 2 :" << mu2_tlv.Pt() << " err "
            << mu2_ptErr_kinfit
            << "  -- Pat muons : " << muonsWithTrack.at(1).pt() << std::endl;
  // std::cout << "Kin fit mass from kinfit: " << higgs_tlv.M()  << " - Kin fit
  // mass from tlv: " << (mu1_tlv+mu2_tlv).M()<< std::endl;

  const DimuonKinematicFitResult kinfit_res(mu1_tlv, mu2_tlv);
  return kinfit_res;
}

// ------------ method called to produce the data  ------------
void MuonKinFit::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  using namespace edm;

  Handle<edm::View<pat::Muon>> pIn;
  iEvent.getByToken(srcCands_, pIn);
  auto muons = *pIn;

  // Muons must have an associated track to run the kinematic fit
  // Otherwise, you will get a segfault
  pat::MuonCollection muonsSelected;
  for (const auto &mu : muons) {
    if (mu.track().isAvailable()) {
      muonsSelected.push_back(mu);
    }
  }

  // Runs the fit on the first two muons
  if (muonsSelected.size() >= 2) {
    const auto kinfit_res =
        dimuon_kinematic_fit(muonsSelected.at(0), muonsSelected.at(1));
    muonsSelected.at(0).addUserFloat("pt_dimuon_kinfit",
                                     kinfit_res.mu1_tlv.Pt());
    muonsSelected.at(1).addUserFloat("pt_dimuon_kinfit",
                                     kinfit_res.mu2_tlv.Pt());

    if (muonsSelected.size() >= 3) {
      std::cerr << "Warning: kinematic fit only calculated for the first two "
                   "muons, but we had "
                << muonsSelected.size() << std::endl;
    }
  }

  // Put output muons to event
  iEvent.put(std::make_unique<pat::MuonCollection>(muonsSelected));
}

// ------------ method called once each stream before processing any runs, lumis
// or events  ------------
void MuonKinFit::beginStream(edm::StreamID) {}

// ------------ method called once each stream after processing all runs, lumis
// and events  ------------
void MuonKinFit::endStream() {}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void MuonKinFit::fillDescriptions(
    edm::ConfigurationDescriptions &descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(MuonKinFit);
