/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
//
// 20140313  Introduce multiple inheritance from G4CMPProcessUtils, also
//	     use proper TransformAxis() on vectors, *not* TransformPoint()
//	     Add wrapper function to compute individual time steps
// 20140314  Fetch propagation parameters from lattice, instead of hardwired
// 20140324  Migrate to use of volume-local field, do coordinate transforms
// 20140325  Move most of time-step calculation to G4CMPProcessUtils
// 20140331  Add required process subtype code
// 20140418  Remove explicit valley transforms, use lattice function
// 20140429  Adjust "effective mass" and energy based on end-of-step momentum
// 20150112  Follow renaming of "SetNewKinematics" to FillParticleChange, drop
//	     redundant IsApplicable()
// 20150122  Use verboseLevel instead of compiler flag for debugging
// 20160624  Use GetTrackInfo() accessor
// 20161114  Use new G4CMPDriftTrackInfo
// 20170602  Use G4CMPUtils for track identity functions
// 20170806  Swap GPIL and MFP functions to work with G4CMPVProcess base
// 20170905  Cache Luke and IV rate models in local LoadDataFromTrack()
// 20170908  Remove "/10." rescaling of field when computing steps
// 20170919  Use rate threshold interface to define alternate step lengths
// 20180831  Fix compiler warning with PostStepDoIt() arguments
// 20190906  Provide functions to externally set rate models, move process
//		lookup functionality to G4CMP(Track)Utils.
// 20200331  C. Stanford (G4CMP-195): Added charge trapping
// 20200331  G4CMP-196: Added impact ionization mean free path
// 20200426  G4CMP-196: Use static function in TrapIonization for MFPs
// 20200504  M. Kelsey (G4CMP-195):  Get trapping MFPs from process
// 20200520  "First report" flag must be thread-local.
// 20200804  Move field access to G4CMPFieldUtils
// 20210923  Ensure that rate calculations are initialized for track
// 20220504  E. Michaud -- Change DBL_MAX to minStep, remove negative MFPs
// 20220729  M. Kelsey -- EnergyStep() has incorrect units for output: should
//		be delta(E)/(q*V).
// 20220730  Drop trapping processes, as they have built-in MFPs, and don't
//		need TimeStepper for energy-dependent calculation.
// 20230527  Drop competitive MFP calculation; use this process to enforce
//	       a maximum allowed step length, to support mass recalculation.

#include "G4CMPTimeStepper.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftTrackInfo.hh"
#include "G4CMPFieldUtils.hh"
#include "G4CMPGeometryUtils.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4DynamicParticle.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4LatticePhysical.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4VPhysicalVolume.hh"
#include <math.h>

G4CMPTimeStepper::G4CMPTimeStepper()
  : G4CMPVDriftProcess("G4CMPTimeStepper", fTimeStepper) {;}

G4CMPTimeStepper::~G4CMPTimeStepper() {;}


// Return configured maximum step length, or DBL_MAX if not set
// NOTE: GPIL is overridden to return exactly this value, no random throw

G4double G4CMPTimeStepper::GetMeanFreePath(const G4Track&, G4double,
					   G4ForceCondition* cond) {
  *cond = NotForced;		// Should we make this forced instead?

  G4double maxStep = G4CMPConfigManager::GetMaximumStep();

  G4ThreadLocal static G4bool first=true;
  if (verboseLevel>1 && first) {
    G4cout << "G4CMPTimeStepper::GetMFP using maxStep " << maxStep << G4endl;
    first = false;
  }
  
  return (maxStep>0. ? maxStep : DBL_MAX);
}


// At end of step, recompute kinematics; important for electrons

G4VParticleChange* G4CMPTimeStepper::PostStepDoIt(const G4Track& aTrack,
						  const G4Step& aStep) {
  if (verboseLevel>1) {
    G4cout << "G4CMPTimeStepper stepLen " << aStep.GetStepLength()/mm << " mm"
	   << G4endl;
  }
  
  aParticleChange.Initialize(aTrack);

  // Adjust mass and kinetic energy using end-of-step momentum
  G4ThreeVector pfinal = GetGlobalMomentum(aTrack);
  FillParticleChange(GetValleyIndex(aTrack), pfinal);

  ClearNumberOfInteractionLengthLeft();		// All processes must do this!
  return &aParticleChange;
}
