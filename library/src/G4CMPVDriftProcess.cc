/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPVDriftProcess.cc
/// \brief Implementation of the G4CMPVDriftProcess base class
//
// $Id$
//
// 20140325  Move time-step calculation here from TimeStepper and LukeScat
// 20140331  Add required process subtype code
// 20141216  Set true electron velocity (and force it) in SetNewKinematics()
// 20150109  Use G4CMP_SET_ELECTRON_MASS to enable dynamical mass, velocity
// 20150111  Add functionality to enforce minimum step length
// 20150112  Handle holes as well as electrons with FillParticleChange()
// 20160624  Use GetTrackInfo() accessor
// 20160829  Drop G4CMP_SET_ELECTRON_MASS code blocks; not physical
// 20161114  Use new G4CMPDriftTrackInfo
// 20170601  Inherit from new G4CMPVProcess, which provides G4CMPProcessUtils
// 20170620  Follow interface changes in G4CMPProcessUtils
// 20201231  FillParticleChange() should also reset valley index if requested

#include "G4CMPVDriftProcess.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhononPolarization.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessType.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "Randomize.hh"

#include "G4CMPFieldUtils.hh"
#include "G4CMPProcessUtils.hh"
#include "G4CMPVScatteringRate.hh"


// Constructor and destructor
// NOTE:  Initial values are arbitrary and non-physical
G4CMPVDriftProcess::G4CMPVDriftProcess(const G4String& processName,
                                       G4CMPProcessSubType stype)
  : G4CMPVProcess(processName, stype) {
  if (verboseLevel) G4cout << GetProcessName() << " is created " << G4endl;
}

G4CMPVDriftProcess::~G4CMPVDriftProcess() {;}


// Only applies to the known charge carriers

G4bool G4CMPVDriftProcess::IsApplicable(const G4ParticleDefinition& aPD) {
  return G4CMP::IsChargeCarrier(aPD);
}


// Overload base version to set a minimum step size, avoiding "stuck" tracks

G4double
G4CMPVDriftProcess::PostStepGetPhysicalInteractionLength(
                      const G4Track& track,
                      G4double previousStepSize,
                      G4ForceCondition* condition) {
  G4double trueLength =
    G4VDiscreteProcess::PostStepGetPhysicalInteractionLength(track,
                                                             previousStepSize,
                                                             condition);

  G4double minLength = G4CMPConfigManager::GetMinStepScale();
  minLength *= (IsElectron() ? theLattice->GetElectronScatter()
		: theLattice->GetHoleScatter());

  G4double ekin = GetKineticEnergy(track);
  G4CMPVScatteringRate* processRate = GetRateModel();
  G4double energyStepMFP = processRate ? EnergyStep(processRate->Threshold(ekin)) : DBL_MAX;
  if (energyStepMFP  <= 1e-9*m) energyStepMFP = 1e-9*m;

  if (verboseLevel > 1) {
    G4cout << GetProcessName() << "::PostStepGPIL: minLength " << minLength
	   << " trueLength " << trueLength << " energyStepMFP " << energyStepMFP << G4endl;
  }
  
  trueLength = energyStepMFP<trueLength ? energyStepMFP : trueLength;
  return minLength<trueLength ? trueLength : minLength;
}


// Fill ParticleChange energy and mass for electron charge carrier momentum

void 
G4CMPVDriftProcess::FillParticleChange(G4int ivalley, const G4ThreeVector& p) {
  // Compute kinetic energy from momentum for electrons or holes
  G4double energy =
    (IsElectron() ? theLattice->MapPtoEkin(ivalley, GetLocalDirection(p))
     : p.mag2()/(2.*GetCurrentTrack()->GetDynamicParticle()->GetMass()) );

  FillParticleChange(ivalley, energy, p);
}

// Fill ParticleChange mass for electron charge carrier with given energy

void G4CMPVDriftProcess::FillParticleChange(G4int ivalley, G4double Ekin,
					    const G4ThreeVector& v) {
  G4CMP::GetTrackInfo<G4CMPDriftTrackInfo>(GetCurrentTrack())->SetValleyIndex(ivalley);

  aParticleChange.ProposeMomentumDirection(v.unit());
  aParticleChange.ProposeEnergy(Ekin);

  if (IsElectron()) {		// Geant4 wants mc^2, not plain mass
    G4double meff = theLattice->GetElectronEffectiveMass(ivalley,
							 GetLocalDirection(v));
    aParticleChange.ProposeMass(meff*c_squared);
  }
}

G4double G4CMPVDriftProcess::EnergyStep(G4double Efinal) const {
  const G4Track* trk = GetCurrentTrack();

  G4double Emag = G4CMP::GetFieldAtPosition(*trk).mag();
  if (Emag <= 0.) return DBL_MAX;		// No field, no acceleration

  G4double Ekin = GetKineticEnergy(trk);
  if (Ekin > Efinal) return DBL_MAX;		// Already over threshold

  if (verboseLevel>1) {
    G4cout << "G4CMPTimeStepper::EnergyStep from " << Ekin/eV
	   << " to " << Efinal/eV << " eV" << G4endl;
  }

  // Add 20% rescaling to account for electron valley systematics
  return 1.2*(Efinal-Ekin)/Emag;
}