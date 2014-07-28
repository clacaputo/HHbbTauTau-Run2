/*!
 * \file GenMatch.h
 * \brief Definition of functions to match gen-level and reco-level particles.
 * \author Riccardo Manzoni (University of Milano-Bicocca, INFN Milano-Bicocca)
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-07-22 created
 *
 * Copyright 2014 Riccardo Manzoni <riccardo.manzoni@cern.ch>,
 *                Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "MCfinalState.h"

// PDG Id: e 11, mu 13, tau 15, Z 23, h 25, H 35, A 36, H+ 37

namespace mcmatching {

  // RM ties together the pdgId and the p4 of a gen particle
  struct matchedGenParticle {
    bool matched ;
    int pdgId ;
    TLorentzVector p4 ;
  };

  struct genjetflavour {
    bool isBjet ;
    std::vector<ntuple::GenParticle> leptonsInJet ;
  };

  struct MChiggses {
    TLorentzVector heavyH  ;
    TLorentzVector lightH1 ;
    TLorentzVector lightH2 ;
  };

  // RM returns a struct containing both the pdgId and the p4 of the gen particle
  // closest to the reco particles passes to it
  matchedGenParticle MatchToGen(analysis::Candidate particle, std::vector<ntuple::GenParticle> genparticles, int status = -1, float dR = 0.3) {

    matchedGenParticle myGenP ;

    myGenP.pdgId = 0 ;
    myGenP.matched = false ;
    myGenP.p4 = particle.momentum ;

    for ( const ntuple::GenParticle genp : genparticles) {
      if (genp.pt == 0) continue ;
      TLorentzVector genpP4 ;
      genpP4.SetPtEtaPhiM(genp.pt, genp.eta, genp.phi, genp.mass) ;
      if (particle.momentum.DeltaR(genpP4) < dR ){
        if (genp.Status != status && status > 0) continue ;
        dR = particle.momentum.DeltaR(genpP4) ;
        myGenP.pdgId   = genp.PdgId ;
        myGenP.p4      = genpP4     ;
        myGenP.matched = true       ;
      }
    }

    return myGenP ;
  };

  genjetflavour IsGenBJet(ntuple::Jet jet, std::vector<ntuple::GenParticle> genparticles, int status = -1, float dR = 0.5) {

    genjetflavour gjf ;

    for ( const ntuple::GenParticle genp : genparticles) {

      if (genp.pt == 0) continue ;
      if (abs(genp.PdgId) != 11 && abs(genp.PdgId) != 13 && abs(genp.PdgId) != 5) continue ;
      if (genp.Status != status && status > 0) continue ;

      TLorentzVector genpP4 ;
      genpP4.SetPtEtaPhiM(genp.pt, genp.eta, genp.phi, genp.mass) ;
      TLorentzVector recoP4 ;
      recoP4.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass) ;

      if ( recoP4.DeltaR(genpP4) < dR ){
        //std::cout << __LINE__ << "]\t" << recoP4.DeltaR(genpP4) << std::endl ;
        if (abs(genp.PdgId) != 5) gjf.isBjet = true ;
        //std::cout << __LINE__ << "]\t" << genp.PdgId << std::endl ;
        if (abs(genp.PdgId) == 11 || abs(genp.PdgId) == 13) gjf.leptonsInJet.push_back(genp) ;
      }
    }
    return gjf ;
  };


  MChiggses GetMChiggses( std::vector<ntuple::GenParticle> genparticles,
                          int heavyHpdgID = 35, int h1DaughtersPdgIDs = 15,
                          int h2DaughtersPdgIDs = 5 ){

    MChiggses mch ;
//     mch.heavyH .SetPtEtaPhiM(1.,1.,1.,1.) ;
//     mch.lightH1.SetPtEtaPhiM(1.,1.,1.,1.) ;
//     mch.lightH2.SetPtEtaPhiM(1.,1.,1.,1.) ;

    for ( const ntuple::GenParticle genp : genparticles) {
      if ( genp.PdgId == heavyHpdgID ) mch.heavyH .SetPtEtaPhiM(genp.pt, genp.eta, genp.phi, genp.mass) ;
      if ( genp.PdgId == h1DaughtersPdgIDs ){
        for ( unsigned int mom : genp.Mother_Indexes ) if ( genparticles.at(mom).PdgId == 25 ) mch.lightH1 .SetPtEtaPhiM(genparticles.at(mom).pt, genparticles.at(mom).eta, genparticles.at(mom).phi, genparticles.at(mom).mass) ;
      }
      if ( genp.PdgId == h2DaughtersPdgIDs ){
        for ( unsigned int mom : genp.Mother_Indexes ) if ( genparticles.at(mom).PdgId == 25 ) mch.lightH2 .SetPtEtaPhiM(genparticles.at(mom).pt, genparticles.at(mom).eta, genparticles.at(mom).phi, genparticles.at(mom).mass) ;
      }
    }

    return mch ;

  };

}
