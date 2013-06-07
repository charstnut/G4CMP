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

#ifndef XWrapperProcess_h
#define XWrapperProcess_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "ChannelingParticleUserInfo.hh"

class G4Material;

class XWrapperProcess : public G4VDiscreteProcess 
{
  public:

  XWrapperProcess(const G4String& processName ="XWrapperProcess" );
  XWrapperProcess(const G4String& processName, G4VDiscreteProcess*);

  virtual ~XWrapperProcess();

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& );
  virtual G4bool IsApplicable(const G4ParticleDefinition&);
  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual G4double PostStepGetPhysicalInteractionLength (const G4Track &track, 
							 G4double previousStepSize, 
							 G4ForceCondition *condition);

  void registerProcess(G4VDiscreteProcess* toRegister);  
  //  void ResetNumberOfInteractionLengthLeft();
  void StartTracking(G4Track* aTrack);
  //  void SubtractNumberOfInteractionLengthLeft(G4double previousStepSize);

  G4double channelingFactor;
                           
  protected:

 
  virtual G4double GetMeanFreePath(
                const G4Track&, G4double, G4ForceCondition* );
 

  private: 
  
  // hide assignment operator as private 
  XWrapperProcess(XWrapperProcess&);
  XWrapperProcess& operator=(const XWrapperProcess& right);

  //private data members
  G4VDiscreteProcess* registeredProcess;
  
};

#endif









