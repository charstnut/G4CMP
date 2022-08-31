/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file library/src/G4CMPPartitionData.cc
/// \brief Definition of the G4CMPPartitionData "hit" (G4VHit) container
///   Simple container class used by G4CMPEnergyPartition to record
///   some of the internal quantities used to generate an event.  Stored
///   into the G4Event via G4CMPPartitionCollection (G4VHitsCollection) 
///   for access by user applications.
///   
// $Id$
//
// 20200218  Michael Kelsey (TAMU) <kelsey@slac.stanford.edu>
// 20200316  Add hit position; improve energy quantity calculations.
// 20210202  Add particle type (PDGcode) to be filled from original track
// 20210327  Add downsampling factors for diagnostics
// 20210820  Add true number of charge pairs and phonons before downsampling
// 20210820  Add estimate of NTL (Luke) emission energy
// 20211030  Add track and step information to support data analysis

#include "globals.hh"
#include "G4CMPPartitionData.hh"
#include "G4SystemOfUnits.hh"


// Initialize allocator

G4ThreadLocal G4Allocator<G4CMPPartitionData>* G4CMPPartitionData_Allocator=0;

// Constructor must initialize everything to zero

G4CMPPartitionData::G4CMPPartitionData()
  : G4VHit(), PDGcode(0), trackID(0), stepID(0), totalEnergy(0.), truedEdx(0.),
    trueNIEL(0.), lindhardYield(0.), FanoFactor(0.), chargeEnergy(0.),
    chargeFano(0.), chargeGenerated(0.), truePairs(0), numberOfPairs(0),
    lukeEnergyEst(0.),
    phononEnergy(0.), phononGenerated(0.), truePhonons(0), numberOfPhonons(0),
    samplingEnergy(0.), samplingCharges(0.), samplingPhonons(0.),
    samplingLuke(0.) {
  position[0]=position[1]=position[2]=position[3]=0.;
}

// Report contents for diagnostic purposes

void G4CMPPartitionData::Print() {	// FIXME: Base class is non-const
  G4cout << "G4CMPPartitionData (particle type " << PDGcode << ","
	 << " track/step " << trackID << "/" << stepID << "):"
         << "\n Total energy deposit " << totalEnergy/eV << " eV"
	 << "\n Position (" << position[0]/mm << "," << position[1]/mm
	 << "," << position[2]/mm << ") mm " << position[3]/ns << " ns"
	 << "\n Ionization energy (dE/dx) " << truedEdx/eV << " eV"
	 << "\n Non-ionizing energy " << trueNIEL/eV << " eV"
	 << "\n Lindhard yield " << lindhardYield
	 << "\n True charge energy " << chargeEnergy/eV << " eV"
	 << "\n Fano factor " << FanoFactor
	 << " fluctuated energy " << chargeFano/eV << " eV"
	 << "\n Generated charge energy " << chargeGenerated/eV << " eV"
	 << "\n Number of charge pairs " << truePairs
	 << " after sampling " << numberOfPairs
	 << "\n Estimated Luke energy " << lukeEnergyEst/eV << " eV"
	 << "\n True phonon energy " << phononEnergy/eV << " eV"
	 << "\n Generated phonon energy " << phononGenerated/eV << " eV"
	 << "\n Number of direct phonons " << truePhonons
	 << " after sampling " << numberOfPhonons
	 << "\n Downsampling energy " << samplingEnergy/eV << " eV"
	 << "\n Downsampling charges " << samplingCharges
	 << "\n Downsampling direct phonons " << samplingPhonons
	 << "\n Downsampling Luke phonons " << samplingLuke
	 << G4endl;
}
