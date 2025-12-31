#ifndef DATAIMPORTER_H
#define DATAIMPORTER_H

#include "globalvars.h"
#include <string>

// Random seed generation
unsigned long generateRandomSeed();

// Parameter loading and initialization
bool updateParametersFromFile(const std::string& filename);
void validateParameters();
void printParameterSummary();
void loadParametersAllocateArraysAndInitializeRandomSeed(const std::string filename);

// Data input functions
void locateHighAndLowResolutionCutoffHalfFill(const clipper::HKL_data<clipper::data64::F_sigF>& reflectionData);
bool loadDiffractionData(const std::string& filename);
void initializeGlobalDensityMapsAndStructureFactors();
bool loadReferenceModelData(const std::string& filename);
bool loadHistogramFatomData(const std::string& filename);
bool loadHistogramFmodelData(const std::string& filename);
bool loadHistogramTxtFile(const std::string& filename);
void loadSigmaWeightArray(const std::string& filename);

// Map initialization
void generateReferenceMapAndMask();
void generateMasksForOrigins(const std::vector<clipper::Coord_frac>& origins,
                             const clipper::Xmap<int>& sourceMask,
                             const clipper::Grid_sampling& grid,
    std::vector<clipper::Xmap<int>>& maskList);
void generatePhasesForOrigins(const std::vector<clipper::Coord_frac>& origins,
                              const clipper::Xmap<float>& sourceDensity,
                              const clipper::Grid_sampling& grid,
                              const clipper::HKL_info& hklTarget,
                              std::vector<clipper::HKL_data<clipper::data64::Phi_fom>>& phaseList);
bool determineInvertedSpaceGroup( int spgNum, clipper::Spacegroup& invertedSpaceGroup,
                                 const clipper::Spacegroup& originalSpaceGroup);
void createInvertedStructureFactors(
                const clipper::HKL_data<clipper::data64::F_phi>& referenceFactors,
                clipper::HKL_data<clipper::data64::F_phi>& fp_inverted);
void printMaskSummary(const clipper::Xmap<int>& mask,
                      const clipper::Xmap<int>& maskInverted,
                      const std::vector<clipper::Coord_frac>& origins);
void printPhaseSummary(const clipper::Xmap<float>& density,
                       const clipper::Xmap<float>& densityInverted,
                       const std::vector<clipper::Coord_frac>& origins);
void generateReferenceMaskForAlternativeOrigins();
void generateReferencePhaseForAlternativeOrigins();

// Initialization modes
void initializeDensityFromRandomMap();
void initializeDensityFromPDBStructure(const std::string& filename, float percentageAtoms = 0.02f);
void initializeDensityFromRandomPhase();
bool initializeDensityFromMtz(const std::string& filename);

#endif // DATAIMPORTER_H
