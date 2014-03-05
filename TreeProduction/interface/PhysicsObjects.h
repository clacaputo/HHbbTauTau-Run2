#ifndef __VHTauTau_TreeMaker_PhysicsObjects_h
#define __VHTauTau_TreeMaker_PhysicsObjects_h

#include <vector>
#include <map>
#include <string>

#include "TObject.h"

namespace vhtm {
  class Event: public TObject {
  public:
    Event();
    ~Event() {}
  
    unsigned int run;
    unsigned int event;
    unsigned int lumis;
    unsigned int bunch;
    unsigned int orbit;
    double time;
    bool isdata;
  
    bool isPhysDeclared;
    bool isBPTX0;
    bool isBSCMinBias;
    bool isBSCBeamHalo;
    bool isPrimaryVertex;
    bool isBeamScraping;
    bool passHBHENoiseFilter;
  
    std::vector<int> nPU;
    std::vector<int> bunchCrossing;
    std::vector<int> trueNInt;
  
    ClassDef(Event,1) 
  };
  class GenEvent: public TObject {
  public:
    GenEvent();
    ~GenEvent() {}
  
    unsigned int processID;
    double ptHat;
    std::vector<double> pdfWeights;
  
    ClassDef(GenEvent,1) 
  };
  
  class Electron: public TObject {
  public:
    Electron();
    ~Electron() {}
    double eta;
    double phi;
    double pt;
      bool ecalDriven;
      bool hasGsfTrack;
    double trackPt;
    double trackPtError;
    double energy;
    double caloEnergy;
    double caloEnergyError;
       int charge;
       int pixHits;
       int trkHits;
       int nValidHits;
     float trkD0;
     float trkD0Error;
  
    // ID variables
    double hoe;
    double hoeDepth1;
    double eop;
    double sigmaEtaEta;
    double sigmaIEtaIEta;
    double deltaPhiTrkSC;
    double deltaEtaTrkSC;
       int classif;
    double e1x5overe5x5;
    double e2x5overe5x5; 
  
    // Iso variables
    double isoEcal03;
    double isoHcal03;
    double isoTrk03;
    double isoEcal04;
    double isoHcal04;
    double isoTrk04;
    double isoRel03;
    double isoRel04;
  
    // Vertex
    double vx;
    double vy;
    double vz;

    // SC associated with electron
    double scEn;
    double scEta;
    double scPhi;
    double scET;
    double scRawEnergy;
  
    // Vertex association variables
    double vtxDist3D;
       int vtxIndex;
    double vtxDistZ;
    double relIso;
    double pfRelIso;
  
     // PFlow isolation variable
     float chargedHadronIso;
     float neutralHadronIso;
     float photonIso;

       int missingHits;

    double dB;
    double edB;
  
    double dB3d;
    double edB3d;
  
    double scE1E9;
    double scS4S1;
    double sckOutOfTime;
    double scEcalIso;
    double scHEEPEcalIso;
    double scHEEPTrkIso;
  
       int nBrems;
     float fbrem;
  
     float dist_vec;
     float dCotTheta;
     float hasMatchedConv;
    
     float mva;
     float mvaPOGTrig;
     float mvaPOGNonTrig;
     bool mvaPreselection;
     bool isTriggerElectron;

     float isoMVA;
   
     float pfRelIso03v1;
     float pfRelIso03v2;
     float pfRelIsoDB03v1;
     float pfRelIsoDB03v2;
     float pfRelIsoDB03v3;

     float pfRelIso04v1;
     float pfRelIso04v2;
     float pfRelIsoDB04v1;
     float pfRelIsoDB04v2;
     float pfRelIsoDB04v3;

     float pfRelIso03;
     float pfRelIso04;
     float pfRelIsoDB03;
     float pfRelIsoDB04;

       int selbit;
       int fidFlag;

    ClassDef(Electron, 1) 
  };
  class GenParticle: public TObject {
  public:
    GenParticle();
    ~GenParticle() {}
    unsigned index;
    double eta;
    double phi;
    double p;
    double px;
    double py;
    double pz;
    double pt;
    double energy;
    int pdgId;
    double vx;
    double vy;
    double vz;
    unsigned status;
    double charge; 
    unsigned numDaught;
    unsigned numMother;
    unsigned motherIndex_1;
    unsigned motherIndex_2;
    ClassDef(GenParticle,2)
  };
  class GenJet: public TObject {
  public:
    GenJet();
    ~GenJet() {}
  
    double eta;
    double phi;
    double p;
    double pt;
    double energy;
    double emf;
    double hadf;
  
    ClassDef(GenJet,1) 
  };
  class MET: public TObject {
  public:
    MET();
    ~MET() {}
  
    double met;
    double metphi;
    double sumet;
    double metuncorr;
    double metphiuncorr;
    double sumetuncorr;
  
    ClassDef(MET,1) 
  };
  class Tau: public TObject {
  public:
    Tau();
    ~Tau() {}
  
    enum {
      kMaxPFChargedCand = 40,
      kMaxPFNeutralCand = 20
    };
    double eta;
    double phi;
    double pt;
    double energy;
    int charge;
    double mass;
  
    double leadTrkPt; 
    double leadTrkPtError; 
    double leadTrkEta;
    double leadTrkPhi;
    double leadTrkCharge;
    double leadTrkD0;
    double leadTrkD0Error;
    double leadTrkDz;
    double leadTrkDzError;

       int vtxIndex;
    double vtxDxy;
    double vtxDz;

    // Leading particle pT
    double leadChargedParticlePt;
    double leadNeutralParticlePt;
    double leadParticlePt;
  
    // Number of charged/neutral candidates and photons in different cones
       int numChargedHadronsSignalCone;
       int numNeutralHadronsSignalCone;
       int numPhotonsSignalCone;
       int numParticlesSignalCone;
       int numChargedHadronsIsoCone;
       int numNeutralHadronsIsoCone;
       int numPhotonsIsoCone;
       int numParticlesIsoCone;
    double ptSumPFChargedHadronsIsoCone;
    double ptSumPFNeutralHadronsIsoCone;
    double ptSumPhotonsIsoCone;
    
    double sigChHadCandPt[kMaxPFChargedCand];
    double sigChHadCandEta[kMaxPFChargedCand];
    double sigChHadCandPhi[kMaxPFChargedCand];
    double sigNeHadCandPt[kMaxPFNeutralCand];
    double sigNeHadCandEta[kMaxPFNeutralCand];
    double sigNeHadCandPhi[kMaxPFNeutralCand];
    double sigGammaCandPt[kMaxPFNeutralCand];
    double sigGammaCandEta[kMaxPFNeutralCand];
    double sigGammaCandPhi[kMaxPFNeutralCand];

    double isoChHadCandPt[kMaxPFChargedCand];
    double isoChHadCandEta[kMaxPFChargedCand];
    double isoChHadCandPhi[kMaxPFChargedCand];
    double isoNeHadCandPt[kMaxPFNeutralCand];
    double isoNeHadCandEta[kMaxPFNeutralCand];
    double isoNeHadCandPhi[kMaxPFNeutralCand];
    double isoGammaCandPt[kMaxPFNeutralCand];
    double isoGammaCandEta[kMaxPFNeutralCand];
    double isoGammaCandPhi[kMaxPFNeutralCand];


    // std::vector<double> ChHadCandPt;
    //std::vector<double> ChHadCandEta;
    //std::vector<double> ChHadCandPhi;
    
    //    double    ChHadCandEnergy[kMaxPFNeutralCand];

    // double ChHadCandPt[6];
    // double ChHadCandEta[6];
    //double ChHadCandPhi[6];


    /////////////////
    float decayModeFinding;
    float   decayModeFindingNewDMs; 
    float   decayModeFindingOldDMs;
    
    float  chargedIsoPtSum ;
    float neutralIsoPtSum;
    float puCorrPtSum;

  double ChHadCandPt3Prong_1track;
  double ChHadCandPt3Prong_2track;
  double ChHadCandPt3Prong_3track;

  double ChHadCandEta3Prong_1track;
  double ChHadCandEta3Prong_2track;
  double ChHadCandEta3Prong_3track;

  double ChHadCandPhi3Prong_1track;
  double ChHadCandPhi3Prong_2track;
  double ChHadCandPhi3Prong_3track;

  double ChHadCandPt1Prong;
  double ChHadCandEta1Prong;
  double ChHadCandPhi1Prong;

  float CombinedIsolationDeltaBetaCorrRaw3Hits;
  float CombinedIsolationDeltaBetaCorrRaw;
  float byIsolationMVA3newDMwLTraw;
  float byIsolationMVA3newDMwoLTraw;
  float byIsolationMVA3oldDMwLTraw;
  float byIsolationMVA3oldDMwoLTraw;

  float byLooseCombinedIsolationDeltaBetaCorr;
  float byLooseCombinedIsolationDeltaBetaCorr3Hits;

  float byMediumCombinedIsolationDeltaBetaCorr;
  float byMediumCombinedIsolationDeltaBetaCorr3Hits;


  float byTightCombinedIsolationDeltaBetaCorr;
  float byTightCombinedIsolationDeltaBetaCorr3Hits;

  float byLooseIsolationMVA3newDMwLT;
  float byLooseIsolationMVA3newDMwoLT;
  float byLooseIsolationMva3oldDMwLT;
  float byLooseIsolationMVA3oldDMwoLT;

  float byMediumIsolationMVA3newDMwLT;

  float byMediumIsolationMVA3newDMwoLT;
  float byMediumIsolationMva3oldDMwLT;
  float byMediumIsolationMVA3oldDMwoLT;

  float byTightIsolationMVA3newDMwLT;
  float byTightIsolationMVA3newDMwoLT;
  float byTightIsolationMva3oldDMwLT;
  float byTightIsolationMVA3oldDMwoLT;

    // discriminators against electrons/muons
  float againstMuonLoose;      
  float againstMuonLoose2;     
  float againstMuonLoose3;     
  float againstMuonLooseMVA;   
  float againstMuonTight;
  
  float againstMuonTight2;    
  float againstMuonTight3;    
  float againstMuonTightMVA;  
  float againstElectronLoose; 
  float againstElectronMedium;
  float againstElectronTight;
 float  pfElectronMVA;        
                            

      // ElectronIDMVA, electron faking tau
  float againstElectronMVALooseMVA5;
  float againstElectronMVAMediumMVA5;
  float againstElectronMVATightMVA5; 

  float byVLooseCombinedIsolationDeltaBetaCorr;
     
   
      // MVA based isolation
  float byLooseIsolationMVA;                   
  float byMediumIsolationMVA;                  
  float byTightIsolationMVA;                  



    //////////// 


     // kinematic variables for PFJet associated to PFTau
    double jetPt;
    double jetEta;
    double jetPhi;
  
     float emFraction;
     float maximumHCALPFClusterEt;
     float ecalStripSumEOverPLead;
     float bremsRecoveryEOverPLead;
     float hcalTotOverPLead;
     float hcalMaxOverPLead;
     float hcal3x3OverPLead;
  
     float etaetaMoment;
     float phiphiMoment;
     float etaphiMoment;
  
    double vx;
    double vy;
    double vz;
  
    double zvertex;
    double ltsipt;
    double    TauVertexPos_x;
    double      TauVertexPos_y;
    double   TauVertexPos_z;
    double     PVPos_x;
    double  PVPos_y;
    double  PVPos_z;
    
    double    FlightLengh;
    double  FlightLenghtSig;
    double  DXY;
    double  DXY_err;
    double      DXY_sig;
    

    int NumChHad;
    int selbit;
  
    ClassDef(Tau, 1) 
  };
  class Muon: public TObject {
  public:
    Muon();
    ~Muon() {}
      bool isTrackerMuon;
      bool isPFMuon;
    double eta;
    double phi;
    double pt;
    double ptError;
    double p;
    double energy;
       int charge;
    double trkD0;
    double trkD0Error;
    double trkDz;
    double trkDzError;
    double globalChi2;
    double trkIso;
    double ecalIso;
    double hcalIso;
    double hoIso;
    double relIso;
       int passID;
    double vtxDist3D;
       int vtxIndex;
    double vtxDistZ;
       int pixHits;
       int trkHits;
       int muoHits;
       int matches;
       int trackerLayersWithMeasurement; 
    double pfRelIso;

    double vx;
    double vy;
    double vz;

    double dB;  // PV2D
    double edB;   

    double dB3d; // PV3D
    double edB3d;
  
    // UW Recommendation
      bool isGlobalMuonPromptTight;
      bool isAllArbitrated;
       int nChambers;
       int nMatches;
       int nMatchedStations;
       unsigned int stationMask;
       unsigned int stationGapMaskDistance;
       unsigned int stationGapMaskPull;
       int muonID;
     
       float idMVA;
       float isoRingsMVA;
       float isoRingsRadMVA;
       float idIsoCombMVA;

       float pfRelIso03v1;
       float pfRelIso03v2;
       float pfRelIsoDB03v1;
       float pfRelIsoDB03v2;
       float pfRelIso04v1;
       float pfRelIso04v2;
       float pfRelIsoDB04v1;
       float pfRelIsoDB04v2;

       int selbit;

     ClassDef(Muon, 1)
  };
  class Jet: public TObject {
  public:
    Jet();
    ~Jet() {}
    double eta;
    double phi;
    double pt;
    double pt_raw;
    double energy;
    double energy_raw;
    double jecUnc;
    double resJEC;
       int partonFlavour;

     float puIdMVA;
       int puIdFlag;
       int puIdBits;

    double chargedEmEnergyFraction;
    double chargedHadronEnergyFraction;
    double chargedMuEnergyFraction;
    double electronEnergyFraction;
    double muonEnergyFraction;
    double neutralEmEnergyFraction;
    double neutralHadronEnergyFraction;
    double photonEnergyFraction;
       int chargedHadronMultiplicity;
       int chargedMultiplicity;
       int electronMultiplicity;
       int muonMultiplicity;
       int neutralHadronMultiplicity;
       int neutralMultiplicity;
       int photonMultiplicity;
       int nConstituents;
    double trackCountingHighEffBTag;
    double trackCountingHighPurBTag;
    double simpleSecondaryVertexHighEffBTag;
    double simpleSecondaryVertexHighPurBTag;
    double jetProbabilityBTag;
    double jetBProbabilityBTag;
    double combinedSecondaryVertexBTag;
    double combinedSecondaryVertexMVABTag;
    double combinedInclusiveSecondaryVertexBTag;
    double combinedMVABTag;
       int passLooseID;
       int passTightID;
  
       int selbit;

    ClassDef(Jet, 1)
  };
  class Vertex: public TObject {
  public:
    Vertex();
    ~Vertex() {}
  
    double x;
    double y;
    double z;
    double xErr;
    double yErr;
    double zErr;
    double rho;
    double chi2;
    double ndf;
       int ntracks;
       int ntracksw05;
      bool isfake;
      bool isvalid;
    double sumPt; // vector sum
  
       int selbit;

    ClassDef(Vertex, 1)
  };
  class GenMET: public TObject {
  public:
    GenMET();
    ~GenMET() {}
  
    double met;
    double metphi;
    double sumet;
  
    ClassDef(GenMET, 1)
  };
  class Track: public TObject {
  public:
    Track();
    ~Track() {}
  
    double eta;
    double etaError;
    double theta;
    double thetaError;
    double phi;
    double phiError;
    double p;
    double pt;
    double ptError;
    double qoverp;
    double qoverpError;
     float charge;

       int nValidHits;
       int nLostHits;
    double validFraction;
       int nValidTrackerHits;
       int nValidPixelHits;
       int nValidStripHits;
       int trackerLayersWithMeasurement;
       int pixelLayersWithMeasurement;
       int stripLayersWithMeasurement;
  
    double dxy;
    double dxyError;
    double dz;
    double dzError;
  
    double chi2;
       int ndof;
    double vx;
    double vy;
    double vz;
  
       int selbit; 
         
    ClassDef(Track, 1)
  };
  class TriggerObject: public TObject {
  public:
    TriggerObject();
    ~TriggerObject() {}
  
    double energy;
    double pt;
    double eta;
    double phi;
  
    std::map< std::string, unsigned int > pathList;
  
    ClassDef(TriggerObject, 1)
  };
}
#endif
