#ifndef GLOBALVARS_H
#define GLOBALVARS_H

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <unistd.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include <cstring>
#include <memory>
#include <mpi.h>
#include <unordered_map>

using namespace std;
using namespace clipper;

namespace GlobalVars {

// =============================================
// PDB Information
// =============================================
extern std::string pdbCode;     
extern std::string pdbHist;     

// =============================================
// Crystallographic Data Structures
// =============================================
extern clipper::HKL_info hklTarget;                     
extern clipper::HKL_info hklHist;                       
extern clipper::HKL_data<clipper::data64::F_sigF> measuredFAll;      
extern clipper::HKL_data<clipper::data64::F_sigF> measuredFAllSigma; 
extern clipper::HKL_data<clipper::data64::F_phi> structureFactors;     
extern clipper::HKL_data<clipper::data64::F_phi> referenceFactors; 
extern clipper::HKL_data<clipper::data64::F_phi> histFatomFactors;
extern clipper::HKL_data<clipper::data64::F_phi> histFmodelFactors;
extern clipper::HKL_data<clipper::data32::Flag> reflectionStatus;  

// Electron density maps
extern clipper::Xmap<float> currentDensity;       
extern clipper::Xmap<float> referenceDensity;  
extern clipper::Xmap<float> tempDensity;  
extern clipper::Xmap<float> densityReserved;
extern clipper::Xmap<float> averagedDensityReserved;

// Maskss
extern clipper::Xmap<int> proteinMask;  // reconstructed protein mask
extern clipper::Xmap<int> histogramMask;  // protein mask of histogram-matching structure
extern clipper::Xmap<int> referenceMask;  // protein mask computed from PDB structure

// 存储所有等效原点的mask
extern std::vector<clipper::Xmap<int>> alternativeOriginMasks;           // 正向结构
extern std::vector<clipper::Xmap<int>> alternativeOriginMasksInverted;   // 反演结构
extern std::vector<clipper::Coord_frac> alternativeOrigins;              // 原点选择
extern std::vector<clipper::Coord_frac> alternativeOriginsInverted;      // 反演结构原点选择
// 存储所有等效原点的phases
extern std::vector<clipper::HKL_data<clipper::data64::Phi_fom>> alternativeOriginPhases;
extern std::vector<clipper::HKL_data<clipper::data64::Phi_fom>> alternativeOriginPhasesInverted;

// =============================================
// Random Number Generation
// =============================================
extern std::uniform_real_distribution<double> uniformDistribution; 
extern std::default_random_engine randomGenerator;                      

// =============================================
// Algorithm Parameters
// =============================================
// Basic parameters
extern int maxIterations;
extern int gridPointCountAsu;
extern int uniqueReflectionCount;
extern float gridSampling;
extern float highResolutionCutoff;
extern float lowResolutionCutoff;
extern float lowResolutionCutoffHalfFill;
extern float histogramResolution;
extern float solventContent;                     
extern float histogramSolventContent;

// HIO parameters
extern float sigmaWeightInitial;               
extern float sigmaWeightFinal;
extern string updateSpeed;
extern float currentSigmaWeight;                  
extern float observedSigmaWeight;                 
extern float hioBeta;
extern float chioAlpha;
extern float hprBeta;
extern float dmBeta;
extern float raarBeta;
extern float mchioBeta;
extern float hioDensityLimitInitial;              
extern float hioDensityLimitFinal;

// THIO, CHIO, HPR, with transition region
extern float shrinkProteinMaskPercentInitial;
extern float shrinkProteinMaskPercentFinal;

// Iteration control parameters
extern int solventFlatteningIterations;                
extern int densityLimitingIterations;               
extern int solventContentVariationIterations;            
extern int sigmaWeightVariationIterations;                           

// =============================================
// Histogram Analysis Parameters
// =============================================
constexpr int SMALL_BIN_COUNT = 10000000;      
extern int smallBinCounts[SMALL_BIN_COUNT]; 

constexpr int LARGE_BIN_COUNT = 300;
extern int largeBinCounts[LARGE_BIN_COUNT]; 
extern float currentBinBoundaries[LARGE_BIN_COUNT + 1];  
extern float standardBinBoundaries[LARGE_BIN_COUNT + 1]; 

constexpr int RESOLUTION_BIN_COUNT = 10;
// Density statistics
extern float densityMin;                      
extern float densityMax;                      
extern float densityAverage;
extern float proteinDensityAverage;

// Histogram matching coefficients
extern float matchingCoeffA[LARGE_BIN_COUNT];             
extern float matchingCoeffB[LARGE_BIN_COUNT];             

// =============================================
// Quality Metrics (using smart pointers)
// =============================================
extern std::unique_ptr<float[]> RfreeValues;                       
extern std::unique_ptr<float[]> RworkValues;
extern std::unique_ptr<float[]> observedSigmaArray;
extern std::unique_ptr<float[]> sigmaWeightArray;
extern std::unique_ptr<float[]> maskCorrelationValues;
extern std::unique_ptr<float[]> maskCorrelationValuesAverageElite;
extern std::unique_ptr<float[]> phaseErrorValues;
extern std::unique_ptr<float[]> phaseErrorValuesAverageElite;

extern std::unique_ptr<float[]> proteinDensityAbsolute;
extern std::unique_ptr<float[]> proteinDensityDeviation;
extern std::unique_ptr<float[]> solventDensityAbsolute;
extern std::unique_ptr<float[]> solventDensityDeviation;

// =============================================
// Parallel Processing (MPI) Parameters
// =============================================
extern int mpiRank;                       
extern int mpiSize;

// =============================================
// Genetic Algorithm Parameters
// =============================================
extern int eliteCount;         // Number of converged/elite members

// Default genetic algorithm parameters
constexpr int DEFAULT_CROSSOVER_POINTS = 11;
constexpr float DEFAULT_CROSSOVER_RATIO = 0.05f;
constexpr float DEFAULT_MUTATION_RATIO = 0.01f;

// =============================================
// Convergence Monitoring
// =============================================s
extern int convergenceFlag;                       

// =============================================
// Utility Variables
// =============================================
extern std::stringstream fileNameStream;   

// =============================================
// Runtime Configuration Flags
// =============================================
extern bool useDataWeightSigma;
extern bool useHybridInputOutput;
extern bool useHIO1;    // 正常的轮廓
extern bool useHIO2;    // 轮廓内小半径精修轮廓
extern bool useHIOold;  // 交换PA、PB，PBPA代替PAPB, 正常的轮廓
extern bool useTransitionHybridInputOutput;  // 不包含论文中严格的内外各5%，值包含轮廓向内5%
extern bool useTHIO1;   // 正常的轮廓内嵌套轮廓，过渡区在蛋白质低密度区域
extern bool useTHIO2;   // 轮廓内小半径精修轮廓得过渡区，过渡区在蛋白质低密度区域
extern bool useTHIOold; // 正常的轮廓内外嵌套两层，过渡区域shirnkPercent*2
extern bool useContinuousHybridInputOutput;
extern bool useCHIO1;   // 正常的轮廓内嵌套轮廓，过渡区在蛋白质低密度区域
extern bool useCHIO2;   // 轮廓内小半径精修轮廓得过渡区，过渡区在蛋白质低密度区域
extern bool useHybridProjectionReflection;
extern bool useHPR1;    // 正常的轮廓内嵌套轮廓，过渡区在蛋白质低密度区域
extern bool useHPR2;    // 轮廓内小半径精修轮廓得过渡区，过渡区在蛋白质低密度区域
extern bool useModifiedContinuousHybridInputOutput;  // 过渡区域反馈比CHIO、HPR更平缓, mchioBeta=0.5, 远小于CHIO的(1-0.4)/0.4=1.5
extern bool useMCHIO1;  // 正常的轮廓内嵌套轮廓，过渡区在蛋白质低密度区域
extern bool useMCHIO2;  // 轮廓内小半径精修轮廓得过渡区，过渡区在蛋白质低密度区域
extern bool useDifferenceMap;  // 蛋白区域，xn
extern bool useDM1;     // 正常的轮廓
extern bool useDM2;     // 轮廓内小半径精修轮廓
extern bool useModifiedDifferenceMap;  // 蛋白区域，PA*PB*xn
extern bool useMDM1;    // 正常的轮廓
extern bool useMDM2;    // 轮廓内小半径精修轮廓
extern bool useAveragedSuccessiveReflections;
extern bool useASR1;    // 正常的轮廓
extern bool useASR2;    // 轮廓内小半径精修轮廓
extern bool useRelaxedAveragedAlternatingReflections;
extern bool useRAAR1;   // 正常的轮廓
extern bool useRAAR2;   // 轮廓内小半径精修轮廓
extern bool useModifiedRelaxedAveragedAlternatingReflections;
extern bool useMRAAR1;  // 正常的轮廓
extern bool useMRAAR2;  // 轮廓内小半径精修轮廓
extern bool useHybridDifferenceMap;
extern bool useHDM1;    // 公式HDM_formula1
extern bool useHDM2;    // 公式HDM_formula2
extern bool useHDM3;    // 公式HDM_formula3
extern bool useHDM4;    // 公式HDM_formula4
extern bool useHDM5;    // 公式HDM_formula5
extern bool useHDM6;    // 公式HDM_formula6
extern bool useGeneticAlgorithm;
extern bool startFromPDB;               
extern bool startFromRandomMap;           
extern bool startFromMR;
extern bool startFromMtz;

// =============================================
// Chiral Space Groups
// =============================================
extern std::unordered_map<int, std::string> spaceGroupSymbols;
extern std::unordered_map<int, int> spaceGroupPairs;

// =============================================
// Memory Management Functions
// =============================================
void allocateArrays();
void deallocateArrays();

} // namespace GlobalVars

#endif // GLOBALVARS_H
