#include "dataimporter.h"
#include "functions.h"
#include "utils.h"
#include <clipper/clipper-ccp4.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <ctime>
#include <unistd.h>
#include <set>
#include <tuple>

using namespace GlobalVars;

unsigned long generateRandomSeed() {
    unsigned long i = clock();
    unsigned long j = time(nullptr);
    unsigned long n = getpid();
    
    // Mix the values
    i = i - j; i = i - n; i = i ^ (n >> 13);
    j = j - n; j = j - i; j = j ^ (i << 8);
    n = n - i; n = n - j; n = n ^ (j >> 13);
    i = i - j; i = i - n; i = i ^ (n >> 12);
    j = j - n; j = j - i; j = j ^ (i << 16);
    n = n - i; n = n - j; n = n ^ (j >> 5);
    i = i - j; i = i - n; i = i ^ (n >> 3);
    j = j - n; j = j - i; j = j ^ (i << 10);
    n = n - i; n = n - j; n = n ^ (j >> 15);
    
    return n;
}

bool updateParametersFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open parameter file " << filename << std::endl;
        return false;
    }

    std::string line;
    int lineNumber = 0;
    int parametersLoaded = 0;

    while (std::getline(file, line)) {
        lineNumber++;
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key;
        if (!(iss >> key)) continue;

        try {
            if (key == "PDB_CODE") {
                std::string value;
                if (iss >> value) {
                    pdbCode = value;
                    parametersLoaded++;
                }
            } else if (key == "PDB_HIST") {
                std::string value;
                if (iss >> value) {
                    pdbHist = value;
                    parametersLoaded++;
                }
            } else if (key == "numIter") {
                if (iss >> maxIterations) parametersLoaded++;
            } else if (key == "resoCutoff") {
                if (iss >> highResolutionCutoff) parametersLoaded++;
            } else if (key == "resoCutoff_hist") {
                if (iss >> histogramResolution) parametersLoaded++;
            } else if (key == "resoCutoffLow") {
                if (iss >> lowResolutionCutoff) parametersLoaded++;
            } else if (key == "resoCutoffLowHalfFill") {
                if (iss >> lowResolutionCutoffHalfFill) parametersLoaded++;
            } else if (key == "solvCont") {
                if (iss >> solventContent) parametersLoaded++;
            } else if (key == "solvCont_hist") {
                if (iss >> histogramSolventContent) parametersLoaded++;
            } else if (key == "sigmaWeight_initial") {
                if (iss >> sigmaWeightInitial) parametersLoaded++;
            } else if (key == "sigmaWeight_final") {
                if (iss >> sigmaWeightFinal) parametersLoaded++;
            } else if (key == "hioBeta") {
                if (iss >> hioBeta) parametersLoaded++;
            } else if (key == "chioAlpha") {
                if (iss >> chioAlpha) parametersLoaded++;
            } else if (key == "hprBeta") {
                if (iss >> hprBeta) parametersLoaded++;
            } else if (key == "dmBeta") {
                if (iss >> dmBeta) parametersLoaded++;
            } else if (key == "raarBeta") {
                if (iss >> raarBeta) parametersLoaded++;
            } else if (key == "mchioBeta") {
                if (iss >> mchioBeta) parametersLoaded++;
            } else if (key == "hioDensLimitIni") {
                if (iss >> hioDensityLimitInitial) parametersLoaded++;
            } else if (key == "hioDensLimitFina") {
                if (iss >> hioDensityLimitFinal) parametersLoaded++;
            } else if (key == "proteinDensAvg") {
                if (iss >> proteinDensityAverage) parametersLoaded++;
            } else if (key == "Obs_sigma") {
                float value;
                if (iss >> value) {
                    useDataWeightSigma = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "hybrid_input_output") {
                float value;
                if (iss >> value) {
                    useHybridInputOutput = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "HIO_k1") {
                float value;
                if (iss >> value) {
                    useHIO1 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "HIO_k2") {
                float value;
                if (iss >> value) {
                    useHIO2 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "HIO_old") {
                float value;
                if (iss >> value) {
                    useHIOold = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "transition_hybrid_input_output") {
                float value;
                if (iss >> value) {
                    useTransitionHybridInputOutput = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "THIO_t1") {
                float value;
                if (iss >> value) {
                    useTHIO1 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "THIO_t2") {
                float value;
                if (iss >> value) {
                    useTHIO2 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "THIO_old") {
                float value;
                if (iss >> value) {
                    useTHIOold = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "continuous_hybrid_input_output") {
                float value;
                if (iss >> value) {
                    useContinuousHybridInputOutput = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "CHIO_t1") {
                float value;
                if (iss >> value) {
                    useCHIO1 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "CHIO_t2") {
                float value;
                if (iss >> value) {
                    useCHIO2 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "hybrid_projection_reflection") {
                float value;
                if (iss >> value) {
                    useHybridProjectionReflection = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "HPR_t1") {
                float value;
                if (iss >> value) {
                    useHPR1 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "HPR_t2") {
                float value;
                if (iss >> value) {
                    useHPR2 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "modified_continuous_hybrid_input_output") {
                float value;
                if (iss >> value) {
                    useModifiedContinuousHybridInputOutput = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "MCHIO_t1") {
                float value;
                if (iss >> value) {
                    useMCHIO1 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "MCHIO_t2") {
                float value;
                if (iss >> value) {
                    useMCHIO2 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "difference_map") {
                float value;
                if (iss >> value) {
                    useDifferenceMap = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "DM_k1") {
                float value;
                if (iss >> value) {
                    useDM1 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "DM_k2") {
                float value;
                if (iss >> value) {
                    useDM2 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "modified_difference_map") {
                float value;
                if (iss >> value) {
                    useModifiedDifferenceMap = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "MDM_k1") {
                float value;
                if (iss >> value) {
                    useMDM1 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "MDM_k2") {
                float value;
                if (iss >> value) {
                    useMDM2 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "averaged_successive_reflections") {
                float value;
                if (iss >> value) {
                    useAveragedSuccessiveReflections = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "ASR_k1") {
                float value;
                if (iss >> value) {
                    useASR1 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "ASR_k2") {
                float value;
                if (iss >> value) {
                    useASR2 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "relaxed_averaged_alternating_reflections") {
                float value;
                if (iss >> value) {
                    useRelaxedAveragedAlternatingReflections = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "RAAR_k1") {
                float value;
                if (iss >> value) {
                    useRAAR1 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "RAAR_k2") {
                float value;
                if (iss >> value) {
                    useRAAR2 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "modified_relaxed_averaged_alternating_reflections") {
                float value;
                if (iss >> value) {
                    useModifiedRelaxedAveragedAlternatingReflections = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "MRAAR_k1") {
                float value;
                if (iss >> value) {
                    useMRAAR1 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "MRAAR_k2") {
                float value;
                if (iss >> value) {
                    useMRAAR2 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "hybrid_difference_map") {
                float value;
                if (iss >> value) {
                    useHybridDifferenceMap = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "HDM_f1") {
                float value;
                if (iss >> value) {
                    useHDM1 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "HDM_f2") {
                float value;
                if (iss >> value) {
                    useHDM2 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "HDM_f3") {
                float value;
                if (iss >> value) {
                    useHDM3 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "HDM_f4") {
                float value;
                if (iss >> value) {
                    useHDM4 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "HDM_f5") {
                float value;
                if (iss >> value) {
                    useHDM5 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "HDM_f6") {
                float value;
                if (iss >> value) {
                    useHDM6 = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "genetic_algorithm") {
                float value;
                if (iss >> value) {
                    useGeneticAlgorithm = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "start_from_pdb") {
                float value;
                if (iss >> value) {
                    startFromPDB = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "start_from_randmap") {
                float value;
                if (iss >> value) {
                    startFromRandomMap = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "start_from_mr") {
                float value;
                if (iss >> value) {
                    startFromMR = (value != 0);
                    parametersLoaded++;
                }
            } else if (key == "start_from_mtz") {
                float value;
                if (iss >> value) {
                    startFromMtz = (value != 0);
                    parametersLoaded++;
                }
            } else {
                std::cerr << "Warning: Unknown parameter '" << key 
                         << "' at line " << lineNumber << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: Invalid value for " << key 
                     << " at line " << lineNumber << std::endl;
        }
    }

    file.close();
    
    std::cout << "Loaded " << parametersLoaded << " parameters from " 
              << filename << std::endl;

    return parametersLoaded > 0;
}

void loadParametersAllocateArraysAndInitializeRandomSeed(const std::string paramFile) {
    // Update parameters
    if (!updateParametersFromFile(paramFile)) {
        throw std::runtime_error("Failed to load parameters");
    }
    // Validate parameters
    validateParameters();
    
    // Allocate memory for global arrays related to maxIterations
    allocateArrays();
    
    // Calculate iteration phase boundaries
    solventFlatteningIterations = static_cast<int>(maxIterations * 0.05);
    densityLimitingIterations = static_cast<int>(maxIterations * 0.15);
    solventContentVariationIterations = maxIterations - densityLimitingIterations - solventFlatteningIterations;
    sigmaWeightVariationIterations = solventContentVariationIterations;
    
    // set grid size  # gridSize = highResolutionCutoff/gridSampling/2
    gridSampling = highResolutionCutoff / 2.0f;  // grid_size = 1.0A
    
//    // set density limit for density modification
//    if (highResolutionCutoff < 2.5f or solventContent > 0.68) {  // 高溶剂含量也不用你限制hioDensity
//        hioDensityLimitInitial = 2.5f;
//        hioDensityLimitFinal = 0.0f;
//    } else {  // resolution~limit: 2.0~2.5, 2.5A~1.5, 3A~1.0, 3.5A~0.5, 3.7A~0.3, 4A~0.3,
//        hioDensityLimitInitial = std::max(0.30f, 1.5f - 1.0f * (highResolutionCutoff - 2.5f));
//        hioDensityLimitFinal = 0.1f;
//    }
    
    // Initialize random number generator
    try {
        std::random_device rd;
        randomGenerator.seed(rd());
    } catch (const std::exception& e) {
        std::cerr << "Error initializing random generator: " << e.what() << std::endl;
        throw;
    }
}

void validateParameters() {
    if (pdbCode.empty()) {
        throw std::invalid_argument("PDB_CODE not specified");
    }
    
    if (maxIterations <= 0) {
        throw std::invalid_argument("maxIterations must be positive");
    }
    
    if (highResolutionCutoff <= 0 || lowResolutionCutoff <= highResolutionCutoff) {
        throw std::invalid_argument("Invalid resolution cutoffs");
    }
    
    if (solventContent < 0.0f || solventContent > 1.0f) {
        throw std::invalid_argument("Solvent content must be between 0 and 1");
    }
    
    if (hioBeta < 0.0f || hioBeta > 2.5f) {
        throw std::invalid_argument("HIO beta must be between 0 and 2.5");
    }
    
    std::cout << "Parameter validation passed.\n";
}

void printParameterSummary() {
    std::cout << "PARAMETERS SUMMARY\n" << std::endl;
    std::cout << "Structure: " << pdbCode;
    if (!pdbHist.empty()) {
        std::cout << " (Histogram: " << pdbHist << ")";
    }
    std::cout << "\n";
    
    std::cout << "Iterations: " << maxIterations << "\n";
    std::cout << "Resolution: " << std::fixed << std::setprecision(3)
    << highResolutionCutoff << " Å (High) / " << lowResolutionCutoff << " Å (Low)\n";
    std::cout << "Reflections will be fully replaced, " << std::fixed << std::setprecision(3)
    << lowResolutionCutoff << " ~ inf " << " Å" << std::endl;
    std::cout << "Reflections will be half replaced, "  << std::fixed << std::setprecision(3)
    << lowResolutionCutoffHalfFill << " ~ " << lowResolutionCutoff << " Å" << std::endl;
    std::cout << "Solvent Content: " << std::fixed << std::setprecision(1)
    << solventContent * 100 << "%\n";
    
    std::cout << "HIO Parameters: β=" << std::setprecision(2) << hioBeta
    << ", Density Limits=[" << hioDensityLimitInitial
    << " → " << hioDensityLimitFinal << "] e-/Å^3\n";
    std::cout << "Difference Map Parameters: β=" << std::setprecision(2) << dmBeta << std::endl;
    std::cout << "Weighted Average Density, Sigma: Initial=" << sigmaWeightInitial
    << ", Final=" << sigmaWeightFinal << " Å\n";
    
    std::cout << "Initialization Mode: ";
    if (startFromPDB) std::cout << "PDB Structure";
    else if (startFromRandomMap) std::cout << "Random Map";
    else if (startFromMR) std::cout << "Molecular Replacement";
    else std::cout << "Random Phases";
    std::cout << "\n";
    
    if (useHybridInputOutput) {
        std::cout << "Hybrid Input Output: Enabled\n";
    } else if (useHybridDifferenceMap) {
        std::cout << "Hybrid Difference Map: Enabled\n";
    } else if (useDifferenceMap) {
        std::cout << "Difference Map: Enabled\n";
    }
    if (useGeneticAlgorithm) {
        std::cout << "Genetic Algorithm: Enabled\n";
    }
    std::cout << "\n---------------------------------------------\n" << std::endl;
}

void locateHighAndLowResolutionCutoffHalfFill(const clipper::HKL_data<clipper::data64::F_sigF>& reflectionData) {
    
    std::vector<ReflectionInfo> sortedReflections;
    sortedReflections.reserve(reflectionData.num_obs()); // 预分配内存
    
    clipper::HKL_data_base::HKL_reference_index ih;
    for (ih = reflectionData.first(); !ih.last(); ih.next()) {
        bool missing = reflectionData[ih].missing();
        float resolution = 1.0f / std::sqrt(ih.invresolsq());
        float f = reflectionData[ih].f();
        float sigf = reflectionData[ih].sigf();
        sortedReflections.push_back(ReflectionInfo(ih.hkl(), missing, resolution, f, sigf));
    }
    // 按分辨率从大到小排序（降序排列函数已在构造体ReflectionInfo中定义），低分辨率在前，高分辨率在后
    std::sort(sortedReflections.begin(), sortedReflections.end());
    
    int obsReflCount = 0;
    for (int irefl = 0; irefl < sortedReflections.size(); ++irefl) {
        if (!sortedReflections[irefl].missing) {
            obsReflCount++;
            if (obsReflCount == 5) {
                lowResolutionCutoff = sortedReflections[irefl].resolution;
            }
            if (obsReflCount == 300) {
                lowResolutionCutoffHalfFill = sortedReflections[irefl].resolution;
                break;
            }
        }
    }
//    obsReflCount = 0;
//    // 那些未观测到的分辨率高于highResolutionCutoff的，也要用计算值填充一些（不超过2000），有助于降低相位偏差
//    const int fillCount = (int)std::min(5000.0f, sortedReflections.size()*0.05f);
//    for (int irefl = sortedReflections.size()-1; irefl >= 0; --irefl) {
//        if (!sortedReflections[irefl].missing and sortedReflections[irefl].resolution > highResolutionCutoff) {
//            obsReflCount++;
//            if (obsReflCount == fillCount) {
//                highResolutionCutoffFill = highResolutionCutoff - (sortedReflections[irefl].resolution - highResolutionCutoff);
//                break;
//            }
//        }
//    }
}


bool loadDiffractionData(const std::string& filename) {
    try {
        // Initialize MTZ file reader
        clipper::CCP4MTZfile mtzReader;
        mtzReader.open_read(filename.c_str());

        // Import basic crystallographic information
        clipper::HKL_info hklInput;
        mtzReader.import_hkl_info(hklInput);

        // Initialize and import structure factor data
        // 当用CCP4 import cif2mtz 之后，mtz中包含了完整的hkl_list, 只有Fp=?时，inputData[hkl].missing()=1，即使Fp=0，inputData[hkl].missing()=0
        clipper::HKL_data<clipper::data64::F_sigF> inputData(hklInput, hklInput.cell());
        //mtzReader.import_hkl_data(inputData, "/*/*/[FOBS,SIGFOBS]"); //"/*/*/[FP,SIGFP]");
        
        // Try importing with FOBS/SIGFOBS first, fallback to FP/SIGFP
        bool imported = false;
        try {
            mtzReader.import_hkl_data(inputData, "/*/*/[FOBS,SIGFOBS]");
            imported = true;
        } catch (...) {
            try {
                mtzReader.import_hkl_data(inputData, "/*/*/[FP,SIGFP]");
                imported = true;
            } catch (...) {
                // Both attempts failed
            }
        }
        if (!imported) {
            mtzReader.close_read();
            throw std::runtime_error(filename + "MTZ文件缺少必需的强度/方差列 (FOBS/SIGFOBS 或 FP/SIGFP)");
        }
        mtzReader.close_read();

        // Pick highResolutionCutoffFill and lowResolutionCutoffHalfFill
        //locateHighAndLowResolutionCutoffHalfFill(inputData);
        
        // Initialize target data structures at specified resolution
        // 如果MTZ文件包含的hkl_list不严谨，hklInput的hkl_list就不严谨，但以下hklTarget包含给定分辨率下完整的hkl_list
        clipper::Resolution effectiveResolution(highResolutionCutoff);
        hklTarget.init(hklInput.spacegroup(), hklInput.cell(), effectiveResolution);
        hklTarget.generate_hkl_list(); //有可能hklTarget.num_reflections() <= hklInput.num_reflections()
        
        // Initialize working data arrays, 空间群没有正确传递给measuredFAll
        measuredFAll.init(hklTarget, hklTarget.cell());
        measuredFAllSigma.init(hklTarget, hklTarget.cell());
        reflectionStatus.init(hklTarget, hklTarget.cell());
        
        //measuredFAll.init(hklTarget.spacegroup(), hklTarget.cell(), hklTarget.hkl_sampling());
        //measuredFAllSigma.init(hklTarget.spacegroup(), hklTarget.cell(), hklTarget.hkl_sampling());

        // clipper库的HKL_data初始化设计存在隐式依赖缺陷. 仅复制晶胞参数和反射列表，未显式绑定空间群对象（需通过hklTarget.spacegroup()主动获取）
        //cout << "hklTarget.spacegroup().symbol_hm() = " << hklTarget.spacegroup().symbol_hm() << endl;
        //cout << "measuredFAll.spacegroup().symbol_hm() = " << measuredFAll.spacegroup().symbol_hm() << endl; //Unknown
        
        int cnt = 0;
        uniqueReflectionCount = 0; // global variable
        // Copy data to working arrays with resolution cutoff
        clipper::HKL_info::HKL_reference_index ih;
        for (ih = hklTarget.first(); !ih.last(); ih.next()) {
            const auto hkl = ih.hkl();
            float resolution = 1.0f / std::sqrt(ih.invresolsq());
            // 没观测到的，以及测量误差大的，统统标记为缺失数据。
            // inputData使用hklInput初始化，当用hklTarget的下标ih索引inputData时，要转换成ih.hkl()。
            // 由于hklInput和hklTarget略有不同，ih.hkl()指标可能超出hklInput索引范围
            if (!inputData[hkl].missing() and resolution >= highResolutionCutoff
                and resolution <= lowResolutionCutoff
                and inputData[hkl].f() > 2.0 * inputData[hkl].sigf()) {
                
                measuredFAll[ih] = inputData[hkl];
                measuredFAllSigma[ih] = inputData[hkl];
                reflectionStatus[ih].flag() = 1; // Work set
                uniqueReflectionCount++;
                cnt++;
                // Randomly assign 1% to free set
                if (static_cast<float>(rand()) / RAND_MAX < 0.01) {
                    reflectionStatus[ih].flag() = 2; // Free set
                }
            } else {
                // Set missing data, 即使measuredFAll[ih].f()赋零，measuredFAll[ih].missing()就不是零了！！
                //measuredFAll[ih].f() = 0.0f;  // measuredFAll[ih].missing() != 0
                //measuredFAll[ih].sigf() = 0.0f;
                //measuredFAllSigma[ih].f() = 0.0f;  // measuredFAllSigma[ih].missing() != 0
                //measuredFAllSigma[ih].sigf() = 0.0f;
                reflectionStatus[ih].flag() = 0; // Work set
            }
        }
        
        // 此处扩展highResolutionCutoff，在fillAndReplacedObservedData时扩展填充数据的范围，有助于降低观测值的相位偏差
        //highResolutionCutoff = highResolutionCutoff - 3.0f*(highResolutionCutoffFill - highResolutionCutoff);
        
        // Log successful import
        const clipper::Spacegroup& spacegroup = hklTarget.spacegroup();
        const clipper::Cell& cell = hklTarget.cell();
        const clipper::Resolution& resolution = hklTarget.resolution();
        
        std::cout << "\n---------------------------------------------\n"
                  << "SUCCESSFUL DATA IMPORT\n"
                  << "\n"
                  << "File:          " << filename << "\n"
                  << "Reflections:   " << hklTarget.num_reflections() << "\n"
                  << "Resolution Observed:    " << std::fixed << std::setprecision(3)
                  << resolution.limit() << " Å\n"
                  << "Observed data completeness: " << std::setprecision(2) << (float)cnt/hklTarget.num_reflections()*100.0 << " %\n\n"
                  << "SPACEGROUP\n"
                  << "Number:        " << spacegroup.spacegroup_number() << "\n"
                  << "Symbol:        " << spacegroup.symbol_hm() << "\n"
                  << "SymOps:        " << spacegroup.num_symops() << " operations\n\n"
                  << "UNIT CELL\n"
                  << "a,b,c (Å):     " << std::setprecision(3)
                  << std::setw(8) << cell.a() 
                  << std::setw(8) << cell.b() 
                  << std::setw(8) << cell.c() << "\n"
                  << "α,β,γ (°):    " << std::setprecision(2)
                  << std::setw(8) << cell.alpha_deg() 
                  << std::setw(8) << cell.beta_deg() 
                  << std::setw(8) << cell.gamma_deg() << "\n"
                  << "Volume (Å³):   " << std::setprecision(1) << cell.volume() << "\n"
                  << "---------------------------------------------\n"
                  << std::endl;
        
        printParameterSummary();
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error loading diffraction data: " << e.what() << std::endl;
        return false;
    }
}

void initializeGlobalDensityMapsAndStructureFactors() {
    clipper::Resolution effectiveResolution(highResolutionCutoff);
    clipper::Grid_sampling gridTarget(hklTarget.spacegroup(), hklTarget.cell(),
                                      effectiveResolution, gridSampling);
    
    // Initialize all density maps
    currentDensity.init(hklTarget.spacegroup(), hklTarget.cell(), gridTarget);
    referenceDensity.init(hklTarget.spacegroup(), hklTarget.cell(), gridTarget);
    tempDensity.init(hklTarget.spacegroup(), hklTarget.cell(), gridTarget);
    averagedDensityReserved.init(hklTarget.spacegroup(), hklTarget.cell(), gridTarget);
    
    // Set to zero with clipper operator=
    currentDensity = 0;
    referenceDensity = 0;
    tempDensity = 0;
    averagedDensityReserved = 0;
    
    // Initialize masks
    proteinMask.init(hklTarget.spacegroup(), hklTarget.cell(), gridTarget);
    referenceMask.init(hklTarget.spacegroup(), hklTarget.cell(), gridTarget);
    
    // Initialize structure factors
    structureFactors.init(hklTarget, hklTarget.cell());
    
    // Count grid points in asymmetric unit
    gridPointCountAsu = 0;
    for (clipper::Xmap_base::Map_reference_index ix = currentDensity.first(); !ix.last(); ix.next()) {
        gridPointCountAsu++;
    }

    std::cout << "Grid dimensions: " << gridTarget.nu() << " x " << gridTarget.nv() << " x " << gridTarget.nw() << std::endl;
    std::cout << "ASU grid points: " << gridPointCountAsu << std::endl;
    std::cout << "ASU grid points: " << gridTarget.nu()*gridTarget.nv()*gridTarget.nw()/hklTarget.spacegroup().num_symops() << " (theoretically)"<< std::endl;
    std::cout << "Density maps initialized.\n" << std::endl;
}

bool loadReferenceModelData(const std::string& filename) {
    try {
        clipper::CCP4MTZfile mtzReader;
        mtzReader.open_read(filename.c_str());
        
        clipper::HKL_info hklInput;
        mtzReader.import_hkl_info(hklInput);
        //mtzReader.import_hkl_data(referenceFactors, "/*/*/[FMODEL,PHIFMODEL]"); //"/*/*/[FC,PHIC]");
        
        // Try importing with FOBS/SIGFOBS first, fallback to FC,PHIC
        bool imported = false;
        try {
            mtzReader.import_hkl_data(referenceFactors, "/*/*/[FMODEL,PHIFMODEL]");
            imported = true;
        } catch (...) {
            try {
                mtzReader.import_hkl_data(referenceFactors, "/*/*/[FC,PHIC]");
                imported = true;
            } catch (...) {
                // Both attempts failed
            }
        }
        if (!imported) {
            mtzReader.close_read();
            throw std::runtime_error(filename + "MTZ文件缺少必需的强度/方差列 (FMODEL,PHIFMODEL 或 FC,PHIC)");
        }
        mtzReader.close_read();
        
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Warning: Could not load reference model data: " << e.what() << std::endl;
        return false;
    }
}

bool loadHistogramFatomData(const std::string& filename) {
    try {
        clipper::CCP4MTZfile mtzReader;
        mtzReader.open_read(filename.c_str());
        
        clipper::HKL_info hklInput;
        mtzReader.import_hkl_info(hklInput);
        //mtzReader.import_hkl_data(histFatomFactors, "/*/*/[FMODEL,PHIFMODEL]"); //"/*/*/[FC,PHIC]");
        
        // Try importing with FOBS/SIGFOBS first, fallback to FC,PHIC
        bool imported = false;
        try {
            mtzReader.import_hkl_data(histFatomFactors, "/*/*/[FMODEL,PHIFMODEL]");
            imported = true;
        } catch (...) {
            try {
                mtzReader.import_hkl_data(histFatomFactors, "/*/*/[FC,PHIC]");
                imported = true;
            } catch (...) {
                // Both attempts failed
            }
        }
        if (!imported) {
            mtzReader.close_read();
            throw std::runtime_error(filename + "MTZ文件缺少必需的强度/方差列 (FMODEL,PHIFMODEL 或 FC,PHIC)");
        }
        mtzReader.close_read();
        
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Warning: Could not load histogram Fatom data: " << e.what() << std::endl;
        return false;
    }
}

bool loadHistogramFmodelData(const std::string& filename) {
    try {
        clipper::CCP4MTZfile mtzReader;
        mtzReader.open_read(filename.c_str());
        
        clipper::HKL_info hklInput;
        mtzReader.import_hkl_info(hklInput);
        
        clipper::HKL_data<clipper::data64::F_phi> inputData(hklInput, hklInput.cell());
        //mtzReader.import_hkl_data(inputData, "/*/*/[FMODEL,PHIFMODEL]"); //"/*/*/[FC,PHIC]");
        
        // Try importing with FOBS/SIGFOBS first, fallback to FC,PHIC
        bool imported = false;
        try {
            mtzReader.import_hkl_data(inputData, "/*/*/[FMODEL,PHIFMODEL]");
            imported = true;
        } catch (...) {
            try {
                mtzReader.import_hkl_data(inputData, "/*/*/[FC,PHIC]");
                imported = true;
            } catch (...) {
                // Both attempts failed
            }
        }
        if (!imported) {
            mtzReader.close_read();
            throw std::runtime_error(filename + "MTZ文件缺少必需的强度/方差列 (FMODEL,PHIFMODEL 或 FC,PHIC)");
        }
        mtzReader.close_read();
        
        hklHist.init(hklInput.spacegroup(), hklInput.cell(), clipper::Resolution(histogramResolution));
        hklHist.generate_hkl_list();
        
        histFmodelFactors.init(hklHist, hklHist.cell());
                
        // Copy data to working arrays with histogramResolution
        clipper::HKL_info::HKL_reference_index ih;
        for (ih = hklHist.first(); !ih.last(); ih.next()) {
            float resolution = 1.0f / std::sqrt(ih.invresolsq());
            // inputData使用hklInput初始化，当用hklTarget的下标ih索引inputData时，要转换成ih.hkl()
            const auto hkl = ih.hkl();
            if (!inputData[hkl].missing() and resolution >= histogramResolution) {
                histFmodelFactors[ih] = inputData[hkl];
            }
        }
        
        // clipper::HKL_data 通常不存储 F(0,0,0)，即使你设置了，在FFT时也可能被忽略。
        // Clipper 在 FFT 计算中默认会移除 F(0,0,0) 项，这是晶体学软件中的常见做法。
//        const float crystal_vol = hklInput.cell().volume();
//        const float F000 = crystal_vol * (1.0f - histogramSolventContent) * (0.43f - 0.33f);
//
//        clipper::HKL_data_base::HKL_reference_index ih;
//        for ( ih = histFmodelFactors.first(); !ih.last(); ih.next() ) {
//            if ((ih.hkl().h() == 0) && (ih.hkl().k() == 0) && (ih.hkl().l() == 0)){
//                histFmodelFactors[ih].f() = F000;
//                break;
//            }
//        }
        
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Warning: Could not load histogram Fmodel data: " << e.what() << std::endl;
        return false;
    }
}

bool loadHistogramTxtFile(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Cannot open reference histogram file: " << filename << std::endl;
        return false;
    }
    
    int binIndex;
    float boundaryValue;
    while (inputFile >> binIndex >> boundaryValue) {
        if (binIndex >= 0 && binIndex <= LARGE_BIN_COUNT) {
            standardBinBoundaries[binIndex] = boundaryValue;
        }
    }
    
    inputFile.close();
    std::cout << "Loaded reference histogram boundaries from: " << filename << std::endl;
    return true;
}

void loadSigmaWeightArray(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Cannot open sigma weight file: " << filename << std::endl;
        return;
    }

    std::string line;
    int iteration;
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        if (iss >> iteration) {
            if (iteration >= 0 && iteration < maxIterations) {
                iss >> observedSigmaArray[iteration];
            }
        }
    }
    inputFile.close();
    std::cout << "Loaded sigma weight arrays from: " << filename << std::endl;
}

//===== 可以删除 =====
void generateReferenceMapAndMask() {
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Generating true density map, mask and inverted mask ..." << std::endl;
    
    // Generate reference density map if reference factors are available
    if (!referenceFactors.first().last()) {
        referenceDensity.fft_from(referenceFactors);
        
        // Generate reference mask
        generateProteinSolventMask(referenceDensity, referenceMask, sigmaWeightFinal, solventContent);
        
        // Generate inverted reference mask
//        clipper::Xmap_base::Map_reference_index ix;
//        for (ix = referenceMask.first(); !ix.last(); ix.next()) {
//            referenceMaskInverted[ix] = (referenceMask[ix] == 1) ? 0 : 1;
//        }
    }
    
    std::cout << "True data generation completed" << std::endl;
}

// 对Xmap应用原点偏移
template<typename T>
clipper::Xmap<T> shiftMapByOrigin(
    const clipper::Xmap<T>& sourceMap,
    const clipper::Coord_frac& origin,
    const clipper::Grid_sampling& grid)
{
    clipper::Xmap<T> shiftedMap(
        sourceMap.spacegroup(),
        sourceMap.cell(),
        sourceMap.grid_sampling());
    
    clipper::Coord_grid origin_grid = origin.coord_grid(grid);
    
    for (auto ix = sourceMap.first(); !ix.last(); ix.next()) {
        clipper::Coord_grid cg = ix.coord();
        
        // 应用偏移
        clipper::Coord_grid cg_shifted(
            cg.u() - origin_grid.u(),
            cg.v() - origin_grid.v(),
            cg.w() - origin_grid.w()
        );
        
        // 处理周期性边界
        clipper::Coord_frac cf_shifted = cg_shifted.coord_frac(grid);
        cf_shifted = cf_shifted.lattice_copy_unit();  // 规范化到 [0, 1)
        clipper::Coord_grid cg_final = cf_shifted.coord_grid(grid);
        
        shiftedMap[ix] = sourceMap.get_data(cg_final);
    }
    
    return shiftedMap;
}

// 为所有原点生成偏移后的masks
void generateMasksForOrigins(
    const std::vector<clipper::Coord_frac>& origins,
    const clipper::Xmap<int>& sourceMask,
    const clipper::Grid_sampling& grid,
    std::vector<clipper::Xmap<int>>& maskList)
{
    maskList.clear();
    maskList.reserve(origins.size());
    
    for (size_t i = 0; i < origins.size(); ++i) {
        std::cout << "  Processing origin " << i << ": ("
                  << origins[i].u() << ", "
                  << origins[i].v() << ", "
                  << origins[i].w() << ")" << std::endl;
        
        clipper::Xmap<int> shiftedMask = shiftMapByOrigin(sourceMask, origins[i], grid);
        maskList.push_back(shiftedMask);
    }
}

// 为所有原点生成相位
void generatePhasesForOrigins(
    const std::vector<clipper::Coord_frac>& origins,
    const clipper::Xmap<float>& sourceDensity,
    const clipper::Grid_sampling& grid,
    const clipper::HKL_info& hklTarget,
    std::vector<clipper::HKL_data<clipper::data64::Phi_fom>>& phaseList)
{
    phaseList.clear();
    phaseList.reserve(origins.size());
    
    for (size_t i = 0; i < origins.size(); ++i) {
        std::cout << "  Processing origin " << i << ": ("
                  << origins[i].u() << ", "
                  << origins[i].v() << ", "
                  << origins[i].w() << ")" << std::endl;
        
        // 偏移密度图
        clipper::Xmap<float> shiftedDensity = shiftMapByOrigin(
            sourceDensity, origins[i], grid);
        
        // 计算结构因子
        clipper::HKL_data<clipper::data64::F_phi> shiftedSF(hklTarget);
        shiftedDensity.fft_to(shiftedSF);
        
        // 提取相位
        clipper::HKL_data<clipper::data64::Phi_fom> phase;
        phase.init(hklTarget, hklTarget.cell());
        
        for (auto ih = hklTarget.first(); !ih.last(); ih.next()) {
            phase[ih].phi() = shiftedSF[ih.hkl()].phi();
            phase[ih].fom() = 1.0f;
        }
        
        phaseList.push_back(phase);
    }
}

// 确定反演结构的空间群
bool determineInvertedSpaceGroup(
    int spgNum,
    clipper::Spacegroup& invertedSpaceGroup,
    const clipper::Spacegroup& originalSpaceGroup)
{
    if (isSpaceGroupAchiral(spgNum)) {
        // 非手性空间群：反演异构且空间群不变
        invertedSpaceGroup = originalSpaceGroup;
        //std::cout << "  Space group is achiral, inverted structure (same space group)" << std::endl;
        return true;
    }
    else if (isSpaceGroupChiral(spgNum)) {
        // 手性空间群：反演异构且需要改变空间群
        auto it = spaceGroupPairs.find(spgNum);
        if (it != spaceGroupPairs.end()) {
            int pairedNumber = it->second;
            clipper::Spgr_descr spgrDescr(pairedNumber);
            invertedSpaceGroup.init(spgrDescr);
            //std::cout << "  Space group is chiral, inverted structure (change space group to paired)" << std::endl;
            return true;
        } else {
            std::cerr << "  Error: Can't find space group pair for chiral space group " << spgNum << std::endl;
            return false;
        }
    }
    else {
        // 中心对称或有镜面对称：反演同构，不需要考虑反演
        std::cout << "  Space group has inversion symmetry, skipping inverted structure" << std::endl;
        return false;
    }
}

// 创建反演结构因子
void createInvertedStructureFactors(
    const clipper::HKL_data<clipper::data64::F_phi>& referenceFactors,
    clipper::HKL_data<clipper::data64::F_phi>& fp_inverted)
{
    for (auto ih = fp_inverted.first(); !ih.last(); ih.next()) {
        fp_inverted[ih].f() = referenceFactors[ih.hkl()].f();
        fp_inverted[ih].phi() = -referenceFactors[ih.hkl()].phi();  // 相位取反
    }
}

// 统计信息输出
void printMaskSummary(
    const clipper::Xmap<int>& mask,
    const clipper::Xmap<int>& maskInverted,
    const std::vector<clipper::Coord_frac>& origins)
{
    // 统计正向mask
    int protein_voxels = 0;
    int total_voxels = 0;
    
    for (auto ix = mask.first(); !ix.last(); ix.next()) {
        total_voxels++;
        if (mask[ix] == 1) protein_voxels++;
    }
    
    // 统计反演mask
    int inverted_protein_voxels = 0;
    int inverted_total_voxels = 0;
    
    for (auto ix = maskInverted.first(); !ix.last(); ix.next()) {
        inverted_total_voxels++;
        if (maskInverted[ix] == 1) inverted_protein_voxels++;
    }
    
    float actual_solvent = 1.0f - float(protein_voxels) / float(total_voxels);
    float actual_solvent_inverted = 1.0f - float(inverted_protein_voxels) / float(inverted_total_voxels);
    
    std::cout << "\n  Summary:" << std::endl;
    std::cout << "  --------" << std::endl;
    std::cout << "  Total ASU voxels (forward): " << total_voxels << std::endl;
    std::cout << "  Total ASU voxels (inverted): " << inverted_total_voxels << std::endl;
    std::cout << "  Protein voxels (forward): " << protein_voxels << std::endl;
    std::cout << "  Protein voxels (inverted): " << inverted_protein_voxels << std::endl;
    std::cout << "  Target solvent content: " << solventContent << std::endl;
    std::cout << "  Actual solvent content (forward): " << actual_solvent << std::endl;
    std::cout << "  Actual solvent content (inverted): " << actual_solvent_inverted << std::endl;
    std::cout << "  Number of alternate origins: " << origins.size() << std::endl;
}

void printPhaseSummary(
    const clipper::Spacegroup& spacegroup,
    const clipper::Spacegroup& spacegroupInverted,
    const std::vector<clipper::Coord_frac>& origins,
    size_t forwardPhaseCount,
    size_t invertedPhaseCount)
{
    std::cout << "\n  Summary:" << std::endl;
    std::cout << "  --------" << std::endl;
    std::cout << "  Space group: " << spacegroup.symbol_hm() << std::endl;
    std::cout << "  Space group inverted: " << spacegroupInverted.symbol_hm() << std::endl;
    std::cout << "  Number of alternate origins: " << origins.size() << std::endl;
    std::cout << "  Forward phases generated: " << forwardPhaseCount << std::endl;
    std::cout << "  Inverted phases generated: " << invertedPhaseCount << std::endl;
}

// 生成参考Mask
void generateReferenceMaskForAlternativeOrigins() {
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Generating true masks for alternative origins ..." << std::endl;
    
    // 检查参考因子是否可用
    if (referenceFactors.first().last()) {
        std::cout << "  Warning: Reference factors are empty, skipping reference data generation" << std::endl;
        std::cout << "\nTrue data generation completed" << std::endl;
        return;
    }
    
    // 1. 生成参考密度图和默认mask
    referenceDensity.fft_from(referenceFactors);
    std::cout << "  True density generated from structure factors" << std::endl;
    
    generateProteinSolventMask(referenceDensity, referenceMask,
                              sigmaWeightFinal, solventContent);
    std::cout << "  True mask generated for default origin" << std::endl;
    
    // 2. 找到所有等效原点
    clipper::Grid_sampling grid = referenceDensity.grid_sampling();
    int spgNum = referenceDensity.spacegroup().spacegroup_number();
    
    std::vector<clipper::Coord_frac> origins;
    catchAlternateOrigins(spgNum, grid, origins);
    alternativeOrigins = origins;
    
    std::cout << "  Found " << origins.size() << " alternate origins for space group "
              << referenceDensity.spacegroup().symbol_hm() << std::endl;
    
    // 3. 为每个原点生成mask
    generateMasksForOrigins(origins, referenceMask, grid, alternativeOriginMasks);
    std::cout << "  Generated " << alternativeOriginMasks.size()
              << " alternate origin masks" << std::endl;
    
    // 4. 处理反演结构
    clipper::Spacegroup invertedSpaceGroup;
    if (!determineInvertedSpaceGroup(spgNum, invertedSpaceGroup,
                                     referenceDensity.spacegroup())) {
        // 不需要反演结构，只输出正向统计
        std::cout << "\n  Summary:" << std::endl;
        std::cout << "  --------" << std::endl;
        
        int protein_voxels = 0;
        int total_voxels = 0;
        for (auto ix = referenceMask.first(); !ix.last(); ix.next()) {
            total_voxels++;
            if (referenceMask[ix] == 1) protein_voxels++;
        }
        
        float actual_solvent = 1.0f - float(protein_voxels) / float(total_voxels);
        std::cout << "  Total ASU voxels: " << total_voxels << std::endl;
        std::cout << "  Protein voxels: " << protein_voxels << std::endl;
        std::cout << "  Target solvent content: " << solventContent << std::endl;
        std::cout << "  Actual solvent content: " << actual_solvent << std::endl;
        std::cout << "  Number of alternate origins: " << origins.size() << std::endl;
        
        std::cout << "\nTrue data generation completed" << std::endl;
        return;
    }
    
    std::cout << "  Generating inverted true masks..." << std::endl;
    
    // 5. 创建反演结构因子
    // 创建反演HKL信息
    clipper::HKL_info invertedHKL;
    clipper::Resolution effectiveResolution(highResolutionCutoff);
    invertedHKL.init(invertedSpaceGroup, referenceDensity.cell(), effectiveResolution);
    invertedHKL.generate_hkl_list();
    clipper::HKL_data<clipper::data64::F_phi> fp_inverted(invertedHKL);
    
    // 创建反演结构因子（相位取反）
    createInvertedStructureFactors(referenceFactors, fp_inverted);
    
    // 6. 生成反演密度图
    clipper::Grid_sampling gridSamplingInverted(
        invertedSpaceGroup, fp_inverted.base_cell(),
        effectiveResolution, gridSampling);
    
    clipper::Xmap<float> referenceDensityInverted(
        invertedSpaceGroup, fp_inverted.base_cell(),
        gridSamplingInverted);
    referenceDensityInverted.fft_from(fp_inverted);
    
    // 7. 生成反演mask
    clipper::Xmap<int> referenceMaskInverted(
        invertedSpaceGroup, fp_inverted.base_cell(),
        gridSamplingInverted);
    generateProteinSolventMask(referenceDensityInverted, referenceMaskInverted,
                              sigmaWeightFinal, solventContent);
    
    std::cout << "  Inverted true mask generated" << std::endl;
    
    // 8. 为反演结构找到等效原点
    int spgNumInverted = referenceDensityInverted.spacegroup().spacegroup_number();
    std::vector<clipper::Coord_frac> originsInverted;
    catchAlternateOrigins(spgNumInverted, gridSamplingInverted, originsInverted);
    alternativeOriginsInverted = originsInverted;
    
    std::cout << "  Found " << originsInverted.size()
              << " alternate origins inverted for space group "
              << referenceDensityInverted.spacegroup().symbol_hm() << std::endl;
    
    // 9. 为反演结构的所有原点生成masks
    generateMasksForOrigins(originsInverted, referenceMaskInverted,
                           gridSamplingInverted, alternativeOriginMasksInverted);
    
    std::cout << "  Generated " << alternativeOriginMasksInverted.size()
              << " alternate origin masks for inverted structure" << std::endl;
    
    // 10. 输出统计信息
    printMaskSummary(referenceMask, referenceMaskInverted, origins);
    
    std::cout << "\nTrue data generation completed" << std::endl;
}

// 生成参考Phase
void generateReferencePhaseForAlternativeOrigins() {
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Generating reference phases from reference density for alternative origins ..." << std::endl;
    
    // 检查参考因子是否可用
    if (referenceFactors.first().last()) {
        std::cout << "  Warning: Reference factors are empty, skipping reference phase generation" << std::endl;
        std::cout << "\nReference phase generation completed" << std::endl;
        return;
    }
    
    // 1. 生成参考密度图
    referenceDensity.fft_from(referenceFactors);
    std::cout << "  Reference density generated from structure factors" << std::endl;
    
    // 2. 找到所有等效原点
    clipper::Grid_sampling grid = referenceDensity.grid_sampling();
    int spgNum = referenceDensity.spacegroup().spacegroup_number();
    
    std::vector<clipper::Coord_frac> origins;
    catchAlternateOrigins(spgNum, grid, origins);
    alternativeOrigins = origins;
    
    std::cout << "  Found " << origins.size() << " alternate origins for space group "
              << referenceDensity.spacegroup().symbol_hm() << std::endl;
    
    // 3. 为每个原点生成相位
    generatePhasesForOrigins(origins, referenceDensity, grid,
                            hklTarget, alternativeOriginPhases);
    
    std::cout << "  Generated " << alternativeOriginPhases.size()
              << " alternate origin phases" << std::endl;
    
    // 4. 处理反演结构
    clipper::Spacegroup invertedSpaceGroup;
    if (!determineInvertedSpaceGroup(spgNum, invertedSpaceGroup,
                                     referenceDensity.spacegroup())) {
        // 不需要反演结构
        std::cout << "\n  Summary:" << std::endl;
        std::cout << "  --------" << std::endl;
        std::cout << "  Space group: " << referenceDensity.spacegroup().symbol_hm() << std::endl;
        std::cout << "  Number of alternate origins: " << origins.size() << std::endl;
        std::cout << "  Forward phases generated: " << alternativeOriginPhases.size() << std::endl;
        
        std::cout << "\nReference phase generation completed" << std::endl;
        return;
    }
    
    std::cout << "  Generating inverted reference phases..." << std::endl;
    
    // 5. 创建反演结构因子
    // 创建反演HKL信息
    clipper::HKL_info invertedHKL;
    clipper::Resolution effectiveResolution(highResolutionCutoff);
    invertedHKL.init(invertedSpaceGroup, referenceDensity.cell(), effectiveResolution);
    invertedHKL.generate_hkl_list();
    clipper::HKL_data<clipper::data64::F_phi> fp_inverted(invertedHKL);
    
    // 创建反演结构因子（相位取反）
    createInvertedStructureFactors(referenceFactors, fp_inverted);
    
    // 6. 生成反演密度图
    clipper::Grid_sampling gridSamplingInverted(
        invertedSpaceGroup, fp_inverted.base_cell(),
        effectiveResolution, gridSampling);
    
    clipper::Xmap<float> referenceDensityInverted(
        invertedSpaceGroup, fp_inverted.base_cell(),
        gridSamplingInverted);
    referenceDensityInverted.fft_from(fp_inverted);
    
    std::cout << "  Inverted reference density generated" << std::endl;
    
    // 7. 为反演结构找到等效原点
    int spgNumInverted = referenceDensityInverted.spacegroup().spacegroup_number();
    std::vector<clipper::Coord_frac> originsInverted;
    catchAlternateOrigins(spgNumInverted, gridSamplingInverted, originsInverted);
    alternativeOriginsInverted = originsInverted;
    
    std::cout << "  Found " << originsInverted.size()
              << " alternate origins inverted for space group "
              << referenceDensityInverted.spacegroup().symbol_hm() << std::endl;
    
    // 8. 为反演结构的所有原点生成相位
    generatePhasesForOrigins(originsInverted, referenceDensityInverted,
                            gridSamplingInverted, hklTarget,
                            alternativeOriginPhasesInverted);
    
    std::cout << "  Generated " << alternativeOriginPhasesInverted.size()
              << " alternate origin phases for inverted structure" << std::endl;
    
    // 9. 输出统计信息
    printPhaseSummary(referenceDensity.spacegroup(),
                     referenceDensityInverted.spacegroup(),
                     origins,
                     alternativeOriginPhases.size(),
                     alternativeOriginPhasesInverted.size());
    
    std::cout << "\nReference phase generation completed" << std::endl;
}

void initializeDensityFromRandomMap() {
    
    clipper::Xmap_base::Map_reference_index ix;
    for (ix = currentDensity.first(); !ix.last(); ix.next()) {
        currentDensity[ix] = static_cast<float>(rand()) / RAND_MAX;
    }
    
    currentDensity.fft_to(structureFactors);
    
    fillAndReplaceObservedData(0, structureFactors);  // iteration=0
    
    currentDensity.fft_from(structureFactors);
    
    // Generate protein mask
    generateProteinSolventMaskUnderControl(0, currentDensity, proteinMask, currentSigmaWeight, solventContent, updateSpeed);  //iteration=0
}

void initializeDensityFromPDBStructure(const std::string& filename, float percentageAtoms) {
    // Reset density map
    clipper::Xmap_base::Map_reference_index ix;
    for (ix = currentDensity.first(); !ix.last(); ix.next()) {
        currentDensity[ix] = 0.05 * static_cast<float>(rand()) / RAND_MAX;
    }

    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Cannot open PDB file: " << filename << std::endl;
        return;
    }

    std::string line;
    int atomsPlaced = 0;
    int totalAtoms = 0;
    
    while (std::getline(inputFile, line)) {
        if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
            totalAtoms++;
            
            float random = static_cast<float>(rand()) / RAND_MAX;
            if (random < percentageAtoms) {
                try {
                    // Parse coordinates
                    float x = std::stof(line.substr(30, 8));
                    float y = std::stof(line.substr(38, 8));
                    float z = std::stof(line.substr(46, 8));
                    
                    // Add small random displacement to avoid grid edge issues
                    x += (static_cast<float>(rand()) / RAND_MAX - 0.5f) * 0.01f;
                    y += (static_cast<float>(rand()) / RAND_MAX - 0.5f) * 0.01f;
                    z += (static_cast<float>(rand()) / RAND_MAX - 0.5f) * 0.01f;

                    // Convert to map coordinates
                    clipper::Coord_orth orthCoord(x, y, z);
                    clipper::Coord_frac fracCoord = orthCoord.coord_frac(currentDensity.cell());
                    clipper::Coord_grid gridCoord = fracCoord.coord_grid(currentDensity.grid_sampling());
                    clipper::Xmap_base::Map_reference_coord mapCoord(currentDensity, gridCoord);

                    currentDensity[mapCoord] += 1.0f;
                    atomsPlaced++;
                } catch (const std::exception& e) {
                    // Skip invalid coordinates
                    continue;
                }
            }
        }
    }
    
    inputFile.close();
    currentDensity.fft_to(structureFactors);
    
    std::cout << "Initialized from PDB structure: " << filename << std::endl;
    std::cout << "Placed " << atomsPlaced << " atoms out of " << totalAtoms
              << " total (" << std::fixed << std::setprecision(1)
              << (100.0f * percentageAtoms) << "% selection)" << std::endl;
    
    fillAndReplaceObservedData(0, structureFactors);  // iteration=0
    
    currentDensity.fft_from(structureFactors);
    
    // Generate protein mask
    generateProteinSolventMaskUnderControl(0, currentDensity, proteinMask, currentSigmaWeight, solventContent, updateSpeed);  //iteration=0
    
}

void initializeDensityFromRandomPhase() {
    
//    // 以下计算出的structureFactors 和把衍射振幅与相位直接赋值给clipper::HKL_data<clipper::data32::F_phi> structureFactors 在本质上没有区别
//    // Initialize random phases
//    clipper::HKL_data<clipper::data64::Phi_fom> startingPhases;
//    startingPhases.init(hklTarget, hklTarget.cell());
//
//    for (clipper::HKL_data_base::HKL_reference_index ih = hklTarget.first(); !ih.last(); ih.next()) {
//        startingPhases[ih].fom() = 1.0f;
//
//        if (ih.hkl_class().centric()) {
//            startingPhases[ih].phi() = (uniformDistribution(randomGenerator) < 0.5) ? 0.0f : clipper::Util::pi();
//        } else {
//            startingPhases[ih].phi() = uniformDistribution(randomGenerator) * clipper::Util::twopi();
//        }
//    }
//
//    // Create initial structure factors with observed amplitudes and random phases
//    clipper::HKL_data<clipper::data64::F_sigF> dataForPhasing;
//    dataForPhasing.init(hklTarget, hklTarget.cell());
//
//    for (clipper::HKL_data_base::HKL_reference_index ih = hklTarget.first(); !ih.last(); ih.next()) {
//        if (!measuredFAll[ih].missing()) {
//            dataForPhasing[ih] = measuredFAll[ih];
//        } else {
//            dataForPhasing[ih].f() = 0.0f;
//            dataForPhasing[ih].sigf() = 0.0f;
//        }
//    }
//    structureFactors.compute(dataForPhasing, startingPhases, clipper::data64::Compute_fphi_from_fsigf_phifom());
    
    // 以下为直接赋值，和上面用structureFactors.compute等效
    for (clipper::HKL_data_base::HKL_reference_index ih = hklTarget.first(); !ih.last(); ih.next()) {
        if (!measuredFAll[ih].missing()) {
            structureFactors[ih].f() = useDataWeightSigma ? measuredFAllSigma[ih].f() : measuredFAll[ih].f();
        } else {
            structureFactors[ih].f() = 0.0f;  // F_phi类型变量，即使f()赋值零，structureFactors[ih].missing() 也会变成0
        }
        if (ih.hkl_class().centric()) {
            structureFactors[ih].phi() = (uniformDistribution(randomGenerator) < 0.5) ? 0.0f : clipper::Util::pi();
        } else {
            structureFactors[ih].phi() = uniformDistribution(randomGenerator) * clipper::Util::twopi();
        }
    }

    currentDensity.fft_from(structureFactors);
    
    fillAndReplaceObservedData(0, structureFactors);  // iteration=0
    
    currentDensity.fft_from(structureFactors);
    
    // Generate protein mask
    generateProteinSolventMaskUnderControl(0, currentDensity, proteinMask, currentSigmaWeight, solventContent, updateSpeed);  //iteration=0
    
}

bool initializeDensityFromMtz(const std::string& filename) {
    try {
        // Initialize MTZ file reader
        clipper::CCP4MTZfile mtzReader;
        mtzReader.open_read(filename.c_str());
        
        // Import structure factor data
        //mtzReader.import_hkl_data(structureFactors, "/*/*/[FP,PHIC]"); //"/*/*/[FMODEL,PHIFMODEL]"); //"/*/*/[FC,PHIC]");

        // Try importing with FOBS/SIGFOBS first, fallback to FC,PHIC
        bool imported = false;
        try {
            mtzReader.import_hkl_data(structureFactors, "/*/*/[FP,PHIC]");
            imported = true;
        } catch (...) {
            try {
                mtzReader.import_hkl_data(structureFactors, "/*/*/[FC,PHIC]");
                imported = true;
            } catch (...) {
                try {
                    mtzReader.import_hkl_data(structureFactors, "/*/*/[FMODEL,PHIFMODEL]");
                    imported = true;
                } catch (...) {
                    // All attempts failed
                }
            }
        }
        if (!imported) {
            mtzReader.close_read();
            throw std::runtime_error(filename + "MTZ文件缺少必需的强度/方差列 (FMODEL,PHIFMODEL 或 FC,PHIC 或 FP,PHIC)");
        }
        mtzReader.close_read();
        
        std::cout << "Initialized density from " << filename << std::endl;
        
        // 针对负密度，flip the density
        for (clipper::HKL_data_base::HKL_reference_index ih = structureFactors.first(); !ih.last(); ih.next()) {
            structureFactors[ih].phi() += clipper::Util::pi(); // 等价于rho = -rho;
        }
        
        fillAndReplaceObservedData(0, structureFactors);  // iteration=0
        
        currentDensity.fft_from(structureFactors);
        
        // Generate protein mask
        generateProteinSolventMaskUnderControl(0, currentDensity, proteinMask, currentSigmaWeight, solventContent, updateSpeed);  //iteration=0
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error loading diffraction data: " << e.what() << std::endl;
        return false;
    }
}
