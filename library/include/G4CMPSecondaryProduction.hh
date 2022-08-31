/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/include/G4CMPSecondaryProduction.hh
/// \brief Definition of the G4CMPSecondaryProduction process class.  This
///	class will be used to generate phonons and charge carriers as
///     secondaries, based on energy loss along step.
//
// $Id$
//
// 20150310  Michael Kelsey
// 20160825  Replace implementation with use of G4CMPEnergyPartition
// 20201207  Add flag to suspend parent track for secondary processing.
// 20210203  G4CMP-241 : Process must run after PostStepDoIt, not AlongStep.
// 20210303  G4CMP-243 : Consolidate nearby steps into one effective hit.
// 20210318  G4CMP-245 : Enforce clearance from crystal surfaces.
// 20210608  G4CMP-260 : Add function to identify steps with energy deposit.
// 20210610  G4CMP-262 : Handle step accumulation including track suspension,
//	       by keeping a map of accumulators by track ID
// 20220216  G4CMP-290 : Add start/end arguments to GeneratePositions().

#ifndef G4CMPSecondaryProduction_hh
#define G4CMPSecondaryProduction_hh 1

#include "G4CMPVProcess.hh"
#include "G4CMPStepAccumulator.hh"
#include "G4ThreeVector.hh"
#include <map>
#include <vector>

class G4CMPEnergyPartition;
class G4DynamicParticle;
class G4ParticleDefinition;
class G4Step;
class G4Track;
class G4VParticleChange;


class G4CMPSecondaryProduction : public G4CMPVProcess {
public:
  G4CMPSecondaryProduction();
  virtual ~G4CMPSecondaryProduction();

  // Applies to all charged, non-resonance particles
  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  // Generate secondaries based on already-computed energy loss
  virtual G4VParticleChange* PostStepDoIt(const G4Track& track,
					   const G4Step& stepData);

  // Overload G4CMPProcessUtils function to fill energy parameters
  virtual void LoadDataForTrack(const G4Track* track);

  // Configurable flag to suspend parent track and process secondaries
  void ProcessSecondariesFirst(G4bool val) { secondariesFirst = val; }
  G4bool ProcessSecondariesFirst() const { return secondariesFirst; }

  // Configurable parameter to specify how to consolidate steps
  void SetCombiningStepLength(G4double val) { combiningStepLength = val; }
  G4bool GetCombiningStepLength() const { return combiningStepLength; }

protected:
  // Step limit for PostStep (sets process Forced for all tracks)
  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

protected:
  G4bool DoAddStep(const G4Step& stepData) const;	// Accumulate step?
  G4bool HasEnergy(const G4Step& stepData) const;	// Does step deposit?
  G4bool DoSecondaries(const G4Step& stepData) const;	// Ready to process?

  void AddSecondaries();		// Convert accumulator with partitioner

  // Create secondaries along the specified trajectory
  void GeneratePositions(size_t npos, const G4ThreeVector& start,
			 const G4ThreeVector& end);

  // Adjust position to be enforced inside current volume
  G4ThreeVector SurfaceClearance(const G4ThreeVector& pos);

public:
  static size_t RandomIndex(size_t imax);	// Used to randomize secondaries

private:
  G4CMPEnergyPartition* partitioner;	// Creates secondary kinematics
  std::vector<G4Track*> theSecs;	// List of created secondaries
  std::vector<G4ThreeVector> posSecs;	// Positions along trajectory
  G4bool secondariesFirst;		// Process secondaries immediately
  G4double combiningStepLength;		// Steps within which to accumulate

  // Collection of accumulators for individual tracks in event
  std::map<G4int, G4CMPStepAccumulator> trackAccum;
  G4CMPStepAccumulator* accumulator;	// Sums multiple steps along track
  G4int currentEventID;			// Remember event for clearing accums

  // No copying allowed
  G4CMPSecondaryProduction(const G4CMPSecondaryProduction& right);
  G4CMPSecondaryProduction& operator=(const G4CMPSecondaryProduction& right);
};

#endif	/* G4CMPSecondaryProduction_hh */
