#include "HHbbTauTau/TreeProduction/interface/PhysicsObjects.h"

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

vhtm::Event::Event() :
  run(0), 
  event(0), 
  lumis(0), 
  bunch(0), 
  orbit(0), 
  time(-1), 
  isdata(false), 
  isPhysDeclared(false),
  isBPTX0(false),
  isBSCMinBias(false),
  isBSCBeamHalo(false),
  isPrimaryVertex(false),
  isBeamScraping(false)
{
  nPU.clear();
  bunchCrossing.clear();
  trueNInt.clear();
}

vhtm::GenEvent::GenEvent() :
  processID(0),
  ptHat(-999)
{
  pdfWeights.clear();
}

vhtm::Electron::Electron() :
  eta(-999),
  phi(-999),
  pt(-999),
  ecalDriven(false),
  hasGsfTrack(false),
  trackPt(-999),
  trackPtError(-999),
  energy(-999),
  caloEnergy(-999),
  caloEnergyError(-999),
  charge(-999),
  pixHits(-1),
  trkHits(-1),
  nValidHits(-1),
  trkD0(-999),
  trkD0Error(-999),
  hoe(-999),
  hoeDepth1(-999),
  eop(-999),
  sigmaEtaEta(-999),
  sigmaIEtaIEta(-999),
  deltaPhiTrkSC(-999),
  deltaEtaTrkSC(-999),
  classif(-1),
  e1x5overe5x5(-999),
  e2x5overe5x5(-99), 
  isoEcal03(-999),
  isoHcal03(-999),
  isoTrk03(-999),
  isoEcal04(-999),
  isoHcal04(-999),
  isoTrk04(-999),
  isoRel03(-999),
  isoRel04(-999),
  vx(-999),
  vy(-999),
  vz(-999),
  scEn(-999),
  scEta(-999),
  scPhi(-999),
  scET(-999),
  scRawEnergy(-999),
  vtxDist3D(-999),
  vtxIndex(-1),
  vtxDistZ(-999),
  relIso(-999), 
  pfRelIso(-999), 
  chargedHadronIso(-999),
  neutralHadronIso(-999),
  photonIso(-999),
  missingHits(-1),
  dB(-999),
  edB(-999),
  dB3d(-999),
  edB3d(-999),
  scE1E9(-999),
  scS4S1(-999),
  sckOutOfTime(-999),
  scEcalIso(-999),
  scHEEPEcalIso(-999),
  scHEEPTrkIso(-999),
  nBrems(-1),
  fbrem(-999),
  dist_vec(-999),
  dCotTheta(-99),
  hasMatchedConv(false),
  mva(-999),
  mvaPOGTrig(-999),
  mvaPOGNonTrig(-999),
  mvaPreselection(false),
  isTriggerElectron(false),
  isoMVA(-999),
  pfRelIso03v1(-999),
  pfRelIso03v2(-999),
  pfRelIsoDB03v1(-999),
  pfRelIsoDB03v2(-999),
  pfRelIsoDB03v3(-999),
  pfRelIso04v1(-999),
  pfRelIso04v2(-999),
  pfRelIsoDB04v1(-999),
  pfRelIsoDB04v2(-999),
  pfRelIsoDB04v3(-999),
  pfRelIso03(-999),
  pfRelIso04(-999),
  pfRelIsoDB03(-999),
  pfRelIsoDB04(-999),
  selbit(0),
  fidFlag(0) 
{}

vhtm::GenParticle::GenParticle() :
  index(-1),
  eta(-999),
  phi(-999),
  p(-999),
  px(-999),
  py(-999),
  pz(-999),
  pt(-999),
  energy(-999),
  pdgId(-999),
  vx(-999),
  vy(-999),
  vz(-999),
  status(-999),
  charge(-999),
  numDaught(-1),
  numMother(-1),
  motherIndex_1(-1),
  motherIndex_2(-1)
{}

vhtm::GenJet::GenJet() :
  eta(-999),
  phi(-999),
  p(-999),
  pt(-999),
  energy(-999),
  emf(-999),
  hadf(-999) 
{}

vhtm::MET::MET() :
  met(-999),
  metphi(-999),
  sumet(-999),
  metuncorr(-999),
  metphiuncorr(-999),
  sumetuncorr(-999) 
{}

vhtm::GenMET::GenMET() :
  met(-999),
  metphi(-999),
  sumet(-999) 
{}

vhtm::Tau::Tau():
  eta(-999),
  phi(-999),
  pt(-999),
  energy(-999),
  charge(-999),
  mass(-999),
  leadTrkPt(-999),
  leadTrkPtError(-999),
  leadTrkEta(-999),
  leadTrkPhi(-999),
  leadTrkCharge(-999),
  leadTrkD0(-999),
  leadTrkD0Error(-999),
  leadTrkDz(-999),
  leadTrkDzError(-999),
  vtxIndex(-1),
  vtxDxy(-999),
  vtxDz(-999),
  leadChargedParticlePt(-999),
  leadNeutralParticlePt(-999),
  leadParticlePt(-999),
  numChargedHadronsSignalCone(-1),
  numNeutralHadronsSignalCone(-1),
  numPhotonsSignalCone(-1),
  numParticlesSignalCone(-1),
  numChargedHadronsIsoCone(-1),
  numNeutralHadronsIsoCone(-1),
  numPhotonsIsoCone(-1),
  numParticlesIsoCone(-1),
  ptSumPFChargedHadronsIsoCone(-999),
  ptSumPFNeutralHadronsIsoCone(-999),
  ptSumPhotonsIsoCone(-999),

  decayModeFinding(-1),
  decayModeFindingNewDMs(-1),
  decayModeFindingOldDMs(-1),
  chargedIsoPtSum(-1),
  neutralIsoPtSum(-1),
  puCorrPtSum(-1),

  ChHadCandPt3Prong_1track(-1),
  ChHadCandPt3Prong_2track(-1),
  ChHadCandPt3Prong_3track(-1),

  ChHadCandEta3Prong_1track(-10),
  ChHadCandEta3Prong_2track(-10),
  ChHadCandEta3Prong_3track(-10),

  ChHadCandPhi3Prong_1track(-10),
  ChHadCandPhi3Prong_2track(-10),
  ChHadCandPhi3Prong_3track(-10),
  
  ChHadCandPt1Prong(-1),
  ChHadCandEta1Prong(-10),
  ChHadCandPhi1Prong(-10),

  CombinedIsolationDeltaBetaCorrRaw3Hits(-1),
  CombinedIsolationDeltaBetaCorrRaw(-1),
  byIsolationMVA3newDMwLTraw(-1),
  byIsolationMVA3newDMwoLTraw(-1),
  byIsolationMVA3oldDMwLTraw(-1),
  byIsolationMVA3oldDMwoLTraw(-1),

  byLooseCombinedIsolationDeltaBetaCorr(-1),
  byLooseCombinedIsolationDeltaBetaCorr3Hits(-1),

  byMediumCombinedIsolationDeltaBetaCorr(-1),
  byMediumCombinedIsolationDeltaBetaCorr3Hits(-1),


  byTightCombinedIsolationDeltaBetaCorr(-1),
  byTightCombinedIsolationDeltaBetaCorr3Hits(-1),

  byLooseIsolationMVA3newDMwLT(-1),
  byLooseIsolationMVA3newDMwoLT(-1),
  byLooseIsolationMva3oldDMwLT(-1),
  byLooseIsolationMVA3oldDMwoLT(-1),

  byMediumIsolationMVA3newDMwLT(-1),

  byMediumIsolationMVA3newDMwoLT(-1),
  byMediumIsolationMva3oldDMwLT(-1),
  byMediumIsolationMVA3oldDMwoLT(-1),

  byTightIsolationMVA3newDMwLT(-1),
  byTightIsolationMVA3newDMwoLT(-1),
  byTightIsolationMva3oldDMwLT(-1),
  byTightIsolationMVA3oldDMwoLT(-1),


      // discriminators against electrons/muons
  againstMuonLoose(-1),      
  againstMuonLoose2(-1),     
  againstMuonLoose3(-1),     
  againstMuonLooseMVA(-1),   
  againstMuonTight(-1),
  
  againstMuonTight2(-1),    
  againstMuonTight3(-1),    
  againstMuonTightMVA(-1),  
  againstElectronLoose(-1), 
  againstElectronMedium(-1),
  againstElectronTight(-1),
  pfElectronMVA(-1),        
                            

      // ElectronIDMVA, electron faking tau
  againstElectronMVALooseMVA5(-1),
  againstElectronMVAMediumMVA5(-1),
  againstElectronMVATightMVA5(-1), 

  byVLooseCombinedIsolationDeltaBetaCorr(-1),
     
   
      // MVA based isolation
  byLooseIsolationMVA(-1),                   
  byMediumIsolationMVA(-1),                  
  byTightIsolationMVA(-1),                  


  jetPt(-999),
  jetEta(-999),
  jetPhi(-999),
  emFraction(-999),
  maximumHCALPFClusterEt(-999),
  ecalStripSumEOverPLead(-999),
  bremsRecoveryEOverPLead(-999),
  hcalTotOverPLead(-999),
  hcalMaxOverPLead(-999),
  hcal3x3OverPLead(-999),
  etaetaMoment(-999),
  phiphiMoment(-999),
  etaphiMoment(-999),
  vx(-999), vy(-999), vz(-999),
  zvertex(-999), ltsipt(-999),

  TauVertexPos_x(-999),
  TauVertexPos_y(-999),
  TauVertexPos_z(-999),

  PVPos_x(-999),
  PVPos_y(-999),
  PVPos_z(-999),

  FlightLengh(-999),
  FlightLenghtSig(-999),
  DXY(-999),
  DXY_err(-999),
  DXY_sig(-999),

  NumChHad(-1),
  selbit(0) 
{

  for (int i = 0; i < kMaxPFChargedCand; ++i) {
    sigChHadCandPt[i] = sigChHadCandEta[i] = sigChHadCandPhi[i] = 
    isoChHadCandPt[i] = isoChHadCandEta[i] = isoChHadCandPhi[i] = -10;

    // ChHadCandPt[i]=-10;
    //ChHadCandEta[i]=-20;
    //ChHadCandPhi[i]=-20;
    //    ChHadCandEnergy[i]=-20;

  
  }
  for (int i = 0; i < kMaxPFNeutralCand; ++i) {
    sigNeHadCandPt[i] = sigNeHadCandEta[i] = sigNeHadCandPhi[i] =
    isoNeHadCandPt[i] = isoNeHadCandEta[i] = isoNeHadCandPhi[i] =
    sigGammaCandPt[i] = sigGammaCandEta[i] = sigGammaCandPhi[i] =
    isoGammaCandPt[i] = isoGammaCandEta[i] = isoGammaCandPhi[i] = -10;
  }
}
vhtm::Muon::Muon() :
  isTrackerMuon(false),
  isPFMuon(false),
  eta(-999),
  phi(-999),
  pt(-999),
  ptError(-999),
  p(-999),
  energy(-999),
  charge(-999),
  trkD0(-999),
  trkD0Error(-999),
  trkDz(-999),
  trkDzError(-999),
  globalChi2(-999),
  trkIso(-999),
  ecalIso(-999),
  hcalIso(-999),
  hoIso(-999),
  relIso(-999),
  passID(false),
  vtxDist3D(-999),
  vtxIndex(-1),
  vtxDistZ(-999),
  pixHits(-1),
  trkHits(-1),
  muoHits(-1),
  matches(-1),
  trackerLayersWithMeasurement(-1),
  pfRelIso(-999),
  vx(-999),
  vy(-999),
  vz(-999),
  dB(-999), 
  edB(-999),
  dB3d(-999),
  edB3d(-999),
  isGlobalMuonPromptTight(false),
  isAllArbitrated(false),
  nChambers(-1),
  nMatches(-1),
  nMatchedStations(-1),
  stationMask(0),
  stationGapMaskDistance(0),
  stationGapMaskPull(0),
  muonID(-1), 
  idMVA(-999),
  isoRingsMVA(-999),
  isoRingsRadMVA(-999),
  idIsoCombMVA(-999),
  pfRelIso03v1(-999),
  pfRelIso03v2(-999),
  pfRelIsoDB03v1(-999),
  pfRelIsoDB03v2(-999),
  pfRelIso04v1(-999),
  pfRelIso04v2(-999),
  pfRelIsoDB04v1(-999),
  pfRelIsoDB04v2(-999),
  selbit(0)
{}

vhtm::Jet::Jet() :
  eta(-999),
  phi(-999),
  pt(-999),
  pt_raw(-999),
  energy(-999),
  energy_raw(-999),
  jecUnc(-999),
  resJEC(-999),
  partonFlavour(-1),
  puIdMVA(-999),
  puIdFlag(-1),
  puIdBits(-1),
  chargedEmEnergyFraction(-999),
  chargedHadronEnergyFraction(-999),
  chargedMuEnergyFraction(-999),
  electronEnergyFraction(-999),
  muonEnergyFraction(-999),
  neutralEmEnergyFraction(-999),
  neutralHadronEnergyFraction(-999),
  photonEnergyFraction(-999),
  chargedHadronMultiplicity(-1),
  chargedMultiplicity(-1),
  electronMultiplicity(-1),
  muonMultiplicity(-1),
  neutralHadronMultiplicity(-1),
  neutralMultiplicity(-1),
  photonMultiplicity(-1),
  nConstituents(-1),
  trackCountingHighEffBTag(-999),
  trackCountingHighPurBTag(-999),
  simpleSecondaryVertexHighEffBTag(-999),
  simpleSecondaryVertexHighPurBTag(-999),
  jetProbabilityBTag(-999),
  jetBProbabilityBTag(-999),
  combinedSecondaryVertexBTag(-999),
  combinedSecondaryVertexMVABTag(-999),
  combinedInclusiveSecondaryVertexBTag(-999),
  combinedMVABTag(-999),
  passLooseID(-1),
  passTightID(-1),
  selbit(0) {}

vhtm::Vertex::Vertex() :
  x(-999),
  y(-999),
  z(-999),
  xErr(-999),
  yErr(-999),
  zErr(-999),
  rho(-999),
  chi2(-999),
  ndf(-1),
  ntracks(-1),
  ntracksw05(-1),
  isfake(true),
  isvalid(false),
  sumPt(-999),
  selbit(0) {}

vhtm::Track::Track()
  : eta(-999),
    etaError(-999),
    theta(-999),
    thetaError(-999),
    phi(-999),
    phiError(-999),
    p(-999),
    pt(-999),
    ptError(-999),
    qoverp(-999),
    qoverpError(-999),
    charge(-999),
    nValidHits(-1),
    nLostHits(-1),
    validFraction(-999),
    nValidTrackerHits(-1),
    nValidPixelHits(-1),
    nValidStripHits(-1),
    trackerLayersWithMeasurement(-1),
    pixelLayersWithMeasurement(-1),
    stripLayersWithMeasurement(-1),
    dxy(-999),
    dxyError(-999),
    dz(-999),
    dzError(-999),
    chi2(-999),
    ndof(-1),
    vx(-999), vy(-999), vz(-999),
    selbit(0) {}

vhtm::TriggerObject::TriggerObject() :
  energy(-999),
  pt(-999),
  eta(-999),
  phi(-999)
{
  pathList.clear();
}
