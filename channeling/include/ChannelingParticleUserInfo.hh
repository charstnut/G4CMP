//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.hh,v 1.1 2010-10-18 15:56:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef ChannelingUserInfo_h
#define ChannelingUserInfo_h 1

#include "globals.hh"
#include "G4VUserTrackInformation.hh"
#include "G4ThreeVector.hh"

class ChannelingParticleUserInfo : public G4VUserTrackInformation
{
    
public:
    
    ChannelingParticleUserInfo();
    ~ChannelingParticleUserInfo();
    
    void SetChanneling(bool flag); 
    bool GetChanneling();
    
    void SetChannelingFactor(G4double);
    G4double GetChannelingFactor();
    G4double GetPreStepChannelingFactor();
    
    G4ThreeVector GetMomentumChanneled();
    void SetMomentumChanneled(G4ThreeVector);
    
    G4ThreeVector GetPositionChanneled();
    void SetPositionChanneled(G4ThreeVector);

    G4ThreeVector GetMomentumChanneledFirst();
    void SetMomentumChanneledFirst(G4ThreeVector);
    
    G4ThreeVector GetPositionChanneledFirst();
    void SetPositionChanneledFirst(G4ThreeVector);

private:
    
    G4double preStepChannelingFactor;
    
    bool channelingFlag; //Has been in channeling in the last step
    G4ThreeVector momentumChanneled; //Last position of the particle in the channel
    G4ThreeVector positionChanneled; //Last projection fof the particle momentum in the crystal reference system
    G4double channelingFactor; //Last value of density seen by channeled particle

    G4ThreeVector momentumChanneledFirst; //Last position of the particle in the channel
    G4ThreeVector positionChanneledFirst; //Last projection fof the particle momentum in the crystal reference system

};



#endif