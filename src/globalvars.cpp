#include "globalvars.h"

namespace GlobalVars {

// =============================================
// PDB Information
// =============================================
std::string pdbCode;         
std::string pdbHist;    

// =============================================
// Crystallographic Data Structures
// =============================================
clipper::HKL_info hklTarget;                          
clipper::HKL_info hklHist;                            
clipper::HKL_data<clipper::data64::F_sigF> measuredFAll;      
clipper::HKL_data<clipper::data64::F_sigF> measuredFAllSigma; 
clipper::HKL_data<clipper::data64::F_phi> structureFactors;          
clipper::HKL_data<clipper::data64::F_phi> referenceFactors;     
clipper::HKL_data<clipper::data64::F_phi> histFatomFactors;
clipper::HKL_data<clipper::data64::F_phi> histFmodelFactors;
clipper::HKL_data<clipper::data32::Flag> reflectionStatus;       

// Density maps
clipper::Xmap<float> currentDensity;          
clipper::Xmap<float> referenceDensity;     
clipper::Xmap<float> tempDensity;
clipper::Xmap<float> densityReserved;
clipper::Xmap<float> averagedDensityReserved;     

// Masks
clipper::Xmap<int> proteinMask;               
clipper::Xmap<int> histogramMask; 
clipper::Xmap<int> referenceMask;

// 存储所有等效原点的mask
std::vector<clipper::Xmap<int>> alternativeOriginMasks;           // 正向结构
std::vector<clipper::Xmap<int>> alternativeOriginMasksInverted;   // 反演结构
std::vector<clipper::Coord_frac> alternativeOrigins;              // 原点选择
std::vector<clipper::Coord_frac> alternativeOriginsInverted;      // 反演结构原点选择

// 存储所有等效原点的phases
std::vector<clipper::HKL_data<clipper::data64::Phi_fom>> alternativeOriginPhases;
std::vector<clipper::HKL_data<clipper::data64::Phi_fom>> alternativeOriginPhasesInverted;

// =============================================
// Random Number Generation
// =============================================
std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);  
std::default_random_engine randomGenerator;  // depend on interpreter and system
//std::mt19937 randomGenerator;  // independent of interpreter and system

// =============================================
// Algorithm Parameters
// =============================================
// Basic iteration and resolution parameters
int maxIterations = 10000;
int gridPointCountAsu = 0;
int uniqueReflectionCount = 0;
float gridSampling = 1.0f;
float highResolutionCutoff = 2.0f;
//float highResolutionCutoffFill = 1.95f;  // highResolutionCutoff-0.05
float lowResolutionCutoff = 20.0f;
float lowResolutionCutoffHalfFill = 15.0f;  // lowResolutionCutoff-5.0
float histogramResolution = 2.0f;                
float solventContent = 0.65f;                  
float histogramSolventContent = 0.65f;   

// Density modification parameters
float sigmaWeightInitial = 4.0f;
float sigmaWeightFinal = 3.0f;
string updateSpeed = "normally";
float hioBeta = 0.75f;
float chioAlpha = 0.4f;  // CHIO
float hprBeta = 0.588f;  // HPR
float dmBeta = 0.75f;
float raarBeta = 1.0f;  // RAAR, 0.95 ~ 1.05
float mchioBeta = 0.5f;  // modified CHIO, less than (1-0.4)/0.4 = 1.5
float hioDensityLimitInitial = 2.5f;
float hioDensityLimitFinal = 0.1f;
float observedSigmaWeight = 0.0f;                 
float currentSigmaWeight = 4.0f;
float proteinDensityAverage = 0.05;

// THIO, CHIO, HPR, with transition region
float shrinkProteinMaskPercentInitial = 0.05f;
float shrinkProteinMaskPercentFinal = 0.0f;

// Iteration control parameters
int solventFlatteningIterations = 100;
int densityLimitingIterations = 100;
int solventContentVariationIterations = 800;
int sigmaWeightVariationIterations = 800;

// Binning parameters for density analysis
int smallBinCounts[SMALL_BIN_COUNT] = {};    
int largeBinCounts[LARGE_BIN_COUNT] = {};    
float currentBinBoundaries[LARGE_BIN_COUNT + 1] = {};     
float standardBinBoundaries[LARGE_BIN_COUNT + 1] = {};  

// Density statistics
float densityMin = 0.0f;                   
float densityMax = 1.0f;                   
float densityAverage = 0.5f;                   

// Histogram matching parameters
float matchingCoeffA[LARGE_BIN_COUNT] = {};          
float matchingCoeffB[LARGE_BIN_COUNT] = {};          

// =============================================
// Quality Metrics (Smart Pointers)
// =============================================
std::unique_ptr<float[]> RfreeValues;
std::unique_ptr<float[]> RworkValues;
std::unique_ptr<float[]> maskCorrelationValues;
std::unique_ptr<float[]> maskCorrelationValuesAverageElite;
std::unique_ptr<float[]> phaseErrorValues;
std::unique_ptr<float[]> observedSigmaArray;
std::unique_ptr<float[]> sigmaWeightArray;
std::unique_ptr<float[]> phaseErrorValuesAverageElite;

std::unique_ptr<float[]> proteinDensityAbsolute;
std::unique_ptr<float[]> proteinDensityDeviation;
std::unique_ptr<float[]> solventDensityAbsolute;
std::unique_ptr<float[]> solventDensityDeviation;

// =============================================
// MPI/Parallel Processing Variables
// =============================================
int mpiRank = 0;
int mpiSize = 1;

// =============================================
// Genetic Algorithm Parameters
// =============================================
int eliteCount = 0;

// =============================================
// Convergence Flags
// =============================================
int convergenceFlag = 0;                     

// =============================================
// Utility Variables
// =============================================
std::stringstream fileNameStream;     

// =============================================
// Runtime Configuration Flags
// =============================================
bool useDataWeightSigma = false;
bool useHybridInputOutput = false;
bool useHIO1 = false;    // 正常的轮廓
bool useHIO2 = false;    // 轮廓内小半径精修轮廓
bool useHIOold = false;  // 交换PA、PB，PBPA代替PAPB, 正常的轮廓
bool useTransitionHybridInputOutput = false;  // 不包含论文中严格的内外各5%，值包含轮廓向内5%
bool useTHIO1 = false;   // 正常的轮廓内嵌套轮廓，过渡区在蛋白质低密度区域
bool useTHIO2 = false;   // 轮廓内小半径精修轮廓得过渡区，过渡区在蛋白质低密度区域
bool useTHIOold = false; // 正常的轮廓内外嵌套两层，过渡区域shrinkPercent*2
bool useContinuousHybridInputOutput = false;
bool useCHIO1 = false;   // 正常的轮廓内嵌套轮廓，过渡区在蛋白质低密度区域
bool useCHIO2 = false;   // 轮廓内小半径精修轮廓得过渡区，过渡区在蛋白质低密度区域
bool useHybridProjectionReflection = false;
bool useHPR1 = false;    // 正常的轮廓内嵌套轮廓，过渡区在蛋白质低密度区域
bool useHPR2 = false;    // 轮廓内小半径精修轮廓得过渡区，过渡区在蛋白质低密度区域
bool useModifiedContinuousHybridInputOutput = false;  // 过渡区域反馈比CHIO、HPR更平缓, mchioBeta=0.5, 远小于CHIO的(1-0.4)/0.4=1.5
bool useMCHIO1 = false;  // 正常的轮廓内嵌套轮廓，过渡区在蛋白质低密度区域
bool useMCHIO2 = false;  // 轮廓内小半径精修轮廓得过渡区，过渡区在蛋白质低密度区域
bool useDifferenceMap = false;  // 蛋白区域，xn
bool useDM1 = false;     // 正常的轮廓
bool useDM2 = false;     // 轮廓内小半径精修轮廓
bool useModifiedDifferenceMap = false;  // 蛋白区域，PA*PB*xn
bool useMDM1 = false;    // 正常的轮廓
bool useMDM2 = false;    // 轮廓内小半径精修轮廓
bool useAveragedSuccessiveReflections = false;
bool useASR1 = false;    // 正常的轮廓
bool useASR2 = false;    // 轮廓内小半径精修轮廓
bool useRelaxedAveragedAlternatingReflections = false;
bool useRAAR1 = false;   // 正常的轮廓
bool useRAAR2 = false;   // 轮廓内小半径精修轮廓
bool useModifiedRelaxedAveragedAlternatingReflections = false;
bool useMRAAR1 = false;   // 正常的轮廓
bool useMRAAR2 = false;   // 轮廓内小半径精修轮廓
bool useHybridDifferenceMap = false;
bool useHDM1 = false;    // 公式HDM_formula1
bool useHDM2 = false;    // 公式HDM_formula2
bool useHDM3 = false;    // 公式HDM_formula3
bool useHDM4 = false;    // 公式HDM_formula4
bool useHDM5 = false;    // 公式HDM_formula5
bool useHDM6 = false;    // 公式HDM_formula6
bool useGeneticAlgorithm = false;
bool startFromPDB = false;          
bool startFromRandomMap = false;      
bool startFromMR = false;
bool startFromMtz = false;


// =============================================
// Chiral Space Groups
// =============================================
// 在65个Sohncke空间群中，有22个组成11对对映体空间群，剩余的43个Sohncke空间群类型虽然不是手性的，但允许手性晶体结构。反演异构，但不改变空间群。
std::unordered_map<int, std::string> spaceGroupSymbols = {
    {76, "P41"},      {78, "P43"},
    {91, "P4122"},    {95, "P4322"},
    {92, "P41212"},   {96, "P43212"},
    {144, "P31"},     {145, "P32"},
    {151, "P3112"},   {153, "P3212"},
    {152, "P3121"},   {154, "P3221"},
    {169, "P61"},     {170, "P65"},
    {171, "P62"},     {172, "P64"},
    {178, "P6122"},   {179, "P6522"},
    {180, "P6222"},   {181, "P6422"},
    {212, "P4332"},   {213, "P4132"}
};  // 使用 auto itsName = spaceGroupSymbols.find(spaceGroupNumber)->second;

// 创建空间群配对映射表
std::unordered_map<int, int> spaceGroupPairs = {
    {76, 78},    {78, 76},
    {91, 95},    {95, 91},
    {92, 96},    {96, 92},
    {144, 145},  {145, 144},
    {151, 153},  {153, 151},
    {152, 154},  {154, 152},
    {169, 170},  {170, 169},
    {171, 172},  {172, 171},
    {178, 179},  {179, 178},
    {180, 181},  {181, 180},
    {212, 213},  {213, 212}
};  // 使用 auto itsNumber = spaceGroupPairs.find(spaceGroupNumber)->second;

// =============================================
// Memory Management Functions
// =============================================
void allocateArrays() {
    try {
        // C++11 compatible allocation using new and unique_ptr reset
        RfreeValues.reset(new float[maxIterations]());  // 括号初始化直接清零
        RworkValues.reset(new float[maxIterations]());  // 括号初始化直接清零
        observedSigmaArray.reset(new float[maxIterations]());
        sigmaWeightArray.reset(new float[maxIterations]());
        maskCorrelationValues.reset(new float[maxIterations]());
        maskCorrelationValuesAverageElite.reset(new float[maxIterations]());
        phaseErrorValues.reset(new float[maxIterations*(RESOLUTION_BIN_COUNT+1)]());
        phaseErrorValuesAverageElite.reset(new float[maxIterations*(RESOLUTION_BIN_COUNT+1)]());
        proteinDensityAbsolute.reset(new float[maxIterations]());
        proteinDensityDeviation.reset(new float[maxIterations]());
        solventDensityAbsolute.reset(new float[maxIterations]());
        solventDensityDeviation.reset(new float[maxIterations]());
    
    } catch (const std::bad_alloc& e) {
        std::cerr << "Memory allocation failed: " << e.what() << std::endl;
        throw;
    }
}

void deallocateArrays() {
    // Smart pointers automatically handle cleanup
    RfreeValues.reset();
    RworkValues.reset();
    observedSigmaArray.reset();
    sigmaWeightArray.reset();
    maskCorrelationValues.reset();
    maskCorrelationValuesAverageElite.reset();
    phaseErrorValues.reset();
    phaseErrorValuesAverageElite.reset();
    proteinDensityAbsolute.reset();
    proteinDensityDeviation.reset();
    solventDensityAbsolute.reset();
    solventDensityDeviation.reset();
}

} // namespace GlobalVars
