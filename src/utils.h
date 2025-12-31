#ifndef UTILS_H
#define UTILS_H

#include "globalvars.h"
#include <string>
#include <vector>

// ===== 存储反射信息，sort reflections =====
struct ReflectionInfo {
    clipper::HKL hkl;
    bool missing;
    float resolution;
    float f;
    float sigf;
    
    // 构造函数
    ReflectionInfo(const clipper::HKL& h, bool miss, float res, float f_val, float sigf_val)
        : hkl(h), missing(miss), resolution(res), f(f_val), sigf(sigf_val) {}
    
    // 比较函数（降序排序）
    bool operator<(const ReflectionInfo& other) const {
        return resolution > other.resolution;
    }
};

// ===== R-factor calculation functions =====
float calculateRfactor(const int observationStatus);
void calculateRfactorWithoutSolventFlattening(const int iteration);
void computeRworkAndRfreeWithSolventFlattening(const int iteration, const clipper::Xmap<float>& density);

// ===== Space group and phase analysis functions =====
bool isSpaceGroupAchiral(const int spaceGroupNumber);
bool isSpaceGroupChiral(const int spaceGroupNumber);
bool isSpaceGroupEnantiomorphic(const int spaceGroupNumber);
void catchAlternateOrigins(const int& spaceGroupNumber, 
                           const clipper::Grid_sampling& grid, 
                           std::vector<clipper::Coord_frac>& origins);

bool findSymEquivReflection(const clipper::HKL_data<clipper::data64::F_phi>& fp,
                           const clipper::HKL& hkl_target,
                           const clipper::Spacegroup& sg,
                           clipper::datatypes::F_phi<float>& result);

void calculateCrossSpectrum(const clipper::HKL_data<clipper::data64::F_phi>& fp1,
                            const clipper::HKL_data<clipper::data64::F_phi>& fp2,
                            const clipper::Spacegroup& spacegroup,
                            clipper::HKL_data<clipper::data64::F_phi>& crossSpectrum,
                            bool invertFp2 = false);

void findMaximumAtOrigins(const clipper::Xmap<float>& translationFunction,
                          const std::vector<clipper::Coord_frac>& origins,
                          const clipper::Grid_sampling& grid,
                          int& bestOriginIndex,
                          float& maxCorrelation,
                          bool isInverted);

void generateOrigins(const clipper::Spacegroup& sg,
                    std::vector<clipper::Coord_frac>& origins);

void findGlobalMaximum(const clipper::Xmap<float>& trans_func,
                      const clipper::Grid_sampling& grid);

clipper::Coord_frac refinePeak(const clipper::Xmap<float>& trans_func,
                              const clipper::Coord_frac& initial,
                              const clipper::Grid_sampling& grid);

void performPhaseTranslation(const clipper::HKL_data<clipper::data64::F_phi>& fp1, 
                             const clipper::HKL_data<clipper::data64::F_phi>& fp2,
                             int& originIndex, bool& isInverted);
float computePhaseAgreement(const clipper::HKL_data<clipper::data64::F_phi>& fp1, 
                            clipper::HKL_data<clipper::data64::F_phi>& fp2,
                            const int& originIndex, const bool& isInverted);
std::vector<double> calculatePhaseError(const clipper::HKL_data<clipper::data64::F_phi>& fp,
                                        const int& originIndex, const bool& isInverted);

// ===== Mask comparison and correlation functions =====
float computeMaskAgreement(const clipper::Xmap<int>& maskA, const clipper::Xmap<int>& maskB,
                           int& originIndex, bool& isInverted);
float calculateTrueMaskCorrelationToFindOrigineChoice(const clipper::Xmap<int>& mask, int& originIndex, bool& isInverted);
void checkAndDealWithInvertedDensityByRvalue(clipper::Xmap<float>& map);
void checkAndDealWithOppositeDensityByProteinDensity(clipper::Xmap<float>& map);
void checkAndDealWithOppositeDensityBySolventDensity(clipper::Xmap<float>& map,
                                                     const clipper::Xmap<int>& mask);

// ===== File output functions =====
void writeHistogramToTxtFile(const std::string& filename);
void writeDensityDistributionFrequencyToCsvFile(const std::string& filename,
                    const clipper::Xmap<float>& density, const clipper::Xmap<int>& mask);
void writeAsuGridToPdbFile(const std::string& filename, const clipper::Xmap<float>& map);
void writeUnitCellGridToPdbFile(const std::string& filename, const clipper::Xmap<float>& map);
void writeUnitCellGridToPdbFileWrong(const std::string& filename, const clipper::Xmap<float>& map);
void writeFloatArrayToTxtFile(const std::string& filename, const float* values);
void writeFloat2DArrayToFile(const std::string& filename, const float* values);

void writeDensityMapAsuToMapFile(const std::string& filename, const clipper::Xmap<float>& map);
void writeDensityMapUnitCellToMapFile(const std::string& filename, const clipper::Xmap<float>& map);
void writeDensityMapUnitCellToMapFile2(const std::string& filename, const clipper::Xmap<float>& map);
void writeDensityMapUnitCellToPdbFile(const std::string& filename, const clipper::Xmap<float>& density, float densityCutoff);

void writeDensityMaskAsuToMapFile(const std::string& filename, const clipper::Xmap<int>& mask);
void writeDensityMaskUnitCellToMapFile(const std::string& filename, const clipper::Xmap<int>& mask);
void writeDensityMaskAsuToPdbFile(const std::string& filename, const clipper::Xmap<int>& mask);
void writeDensityMaskUnitCellToPdbFile(const std::string& filename, const clipper::Xmap<int>& mask);

void writeStructureFactorsToMtzFile(const std::string& filename, const clipper::Xmap<float>& map);
void expandToFullStructureFactors(const clipper::HKL_data<clipper::data64::F_phi>& structureFactors,
                                  clipper::HKL_data<clipper::data64::F_phi>& fullSF);
void expandAndWriteToCIF(const std::string& filename,
                         const clipper::HKL_data<clipper::data64::F_phi>& structureFactors);
void writeFullStructureFactorsToCifFile(const std::string& filename,
                                        const clipper::HKL_data<clipper::data64::F_phi>& structureFactors);
void writeCIFFileForVisualization(const std::string& filename,
                                  const clipper::HKL_data<clipper::data64::F_phi>& structureFactors);

// ===== Density analysis functions =====
float calculateDensityCutoff(const clipper::Xmap<float>& density, float proteinFraction);
void calculateProteinDensityDifference(const int iteration,
                                       const clipper::Xmap<float>& density,
                                       const clipper::Xmap<float>& densityA,
                                       const clipper::Xmap<float>& densityB);
void calculateSolventDensityDifference(const int iteration,
                                       const clipper::Xmap<float>& density,
                                       const clipper::Xmap<float>& densityA,
                                       const clipper::Xmap<float>& densityB);

// ===== Histogram output functions =====
void writeBinBoundariesToFile(const std::string& filename);
void writeHistogramLogToFile(const std::string& filename);

// ===== Check ASU / Unit Cell =====
void checkAsuMap(const clipper::Xmap<float>& map);
void checkUnitCellMap(const clipper::Xmap<float>& map);

// ===== 生成所有原点选择和反演的元胞胞pdb文件, 以及用于Coot展示的pdb文件 =====
bool extractSpaceGroupAndCellFromPDB(const std::string& pdbFile,
                                     clipper::Spacegroup& spacegroup,
                                     clipper::Cell& cell);
std::string applyOriginShiftToAtom(const std::string& atomString,
                                   const clipper::Coord_frac& originShift,
                                   const clipper::Cell& cell,
                                   const bool insideUnitCell);
std::string applyInversionToAtom(const std::string& atomString,
                                 const clipper::Cell& cell,
                                 const bool insideUnitCell);
int countNumAtomModel(const std::string& pdbFile);
void inputAtomString(const std::string& pdbFile, std::string* atomStringArray);
void generateUnitCellFromASU(const std::string& outputFile,
                             const std::string* atomStringArray, int numAtoms,
                             const clipper::Spacegroup& spacegroup,
                             const clipper::Cell& cell);
std::string applySymopToAtom(const std::string& atomString, const clipper::Symop& symop,
                             const clipper::Cell& cell, int newAtomSerial);
void generateUnitCellForAllOrigTranAndInversion(const std::string& inputPdbFile);
void generateAllOrigTranAndInversionForCoot(const std::string& inputPdbFile);

void writeUniqReflectionsToTxtFile(const string mtzFilename, const string outputFilename);

// ===== Compute data-weighting sigma======
void computeSigmaObsForDataWeighting(const clipper::HKL_data<clipper::data64::F_sigF>& reflectionData);

void displayResolutionShell();

void computeDensityRelatedValues(const clipper::Xmap<float>& density, const clipper::Xmap<int>& mask);

void outputRvalueInResolutionShell(const std::string& FobsFileName, const std::string& FmodelFilename, const std::string& outputFilename) ;

#endif // UTILS_H
