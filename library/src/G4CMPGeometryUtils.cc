/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id$
// File: G4CMPGeometryUtils.cc
//
// Description: Free standing helper functions for geometry based calculations.
//
// 20161107  Rob Agnese

#include "G4CMPGeometryUtils.hh"

#include "G4CMPGlobalLocalTransformStore.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4Step.hh"
#include "G4TransportationManager.hh"

G4ThreeVector G4CMP::GetLocalDirection(const G4VPhysicalVolume* pv,
                                       const G4ThreeVector& dir) {
  return G4CMPGlobalLocalTransformStore::ToLocal(pv).TransformAxis(dir);
}

G4ThreeVector G4CMP::GetLocalPosition(const G4VPhysicalVolume* pv,
                                      const G4ThreeVector& pos) {
  return G4CMPGlobalLocalTransformStore::ToLocal(pv).TransformPoint(pos);
}

G4ThreeVector G4CMP::GetGlobalDirection(const G4VPhysicalVolume* pv,
                                        const G4ThreeVector& dir) {
  return G4CMPGlobalLocalTransformStore::ToGlobal(pv).TransformAxis(dir);
}

G4ThreeVector G4CMP::GetGlobalPosition(const G4VPhysicalVolume* pv,
                                       const G4ThreeVector& pos) {
  return G4CMPGlobalLocalTransformStore::ToGlobal(pv).TransformPoint(pos);
}

void G4CMP::RotateToLocalDirection(const G4VPhysicalVolume* pv,
                                   G4ThreeVector& dir) {
  G4CMPGlobalLocalTransformStore::ToLocal(pv).ApplyAxisTransform(dir);
}

void G4CMP::RotateToLocalPosition(const G4VPhysicalVolume* pv,
                                  G4ThreeVector& pos) {
  G4CMPGlobalLocalTransformStore::ToLocal(pv).ApplyPointTransform(pos);
}

void G4CMP::RotateToGlobalDirection(const G4VPhysicalVolume* pv,
                                    G4ThreeVector& dir) {
  G4CMPGlobalLocalTransformStore::ToGlobal(pv).ApplyAxisTransform(dir);
}

void G4CMP::RotateToGlobalPosition(const G4VPhysicalVolume* pv,
                                   G4ThreeVector& pos) {
  G4CMPGlobalLocalTransformStore::ToGlobal(pv).ApplyPointTransform(pos);
}

G4ThreeVector G4CMP::GetSurfaceNormal(const G4Step& step) {
  // Get outward normal using G4Navigator method (more reliable than G4VSolid)
  G4int navID = G4ParallelWorldProcess::GetHypNavigatorID();
  std::vector<G4Navigator*>::iterator iNav =
    G4TransportationManager::GetTransportationManager()
      ->GetActiveNavigatorsIterator();

  G4bool goodNorm;
  G4ThreeVector surfNorm = iNav[navID]->GetGlobalExitNormal(
                                      step.GetPostStepPoint()->GetPosition(),
                                      &goodNorm);

  // FIXME:  Sometimes G4Navigator fails, but still returns "good"
  if (!goodNorm || surfNorm.mag()<0.99) {
    G4VPhysicalVolume* thePrePV = step.GetPreStepPoint()->GetPhysicalVolume();
    G4VPhysicalVolume* thePostPV = step.GetPostStepPoint()->GetPhysicalVolume();
    G4Exception("G4CMPProcessUtils::GetSurfaceNormal", "Boundary001",
                EventMustBeAborted, ("Can't get normal vector of surface between " +
                                    thePrePV->GetName() + " and " +
                                    thePostPV->GetName()+ ".").c_str());
  }

  return surfNorm;
}