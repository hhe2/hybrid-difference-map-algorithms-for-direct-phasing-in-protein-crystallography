#include "dataimporter.h"
#include "functions.h"
#include "utils.h"
#include <mpi.h>
#include <ctime>

using namespace GlobalVars;

void initializeDataProcessing() {
    // Load parameters from parameters.txt
    loadParametersAllocateArraysAndInitializeRandomSeed("parameters.txt");
    
    // -------------------------------------------------------------------------
    // compute R value of the measured data vs Fmodel when pdb is available
    //outputRvalueInResolutionShell("./PDB_file/" + pdbCode + "-sf.mtz",
    //                              "./PDB_file/" + pdbCode + "_fmodel.mtz",
    //                              "./PDB_file/" + pdbCode + "_R_Fobs_Fmodel.txt");
    //
    // -------------------------------------------------------------------------
    
    // Load diffraction data
    if (!loadDiffractionData("./PDB_file/" + pdbCode + "-sf.mtz")) {
        throw std::runtime_error("Failed to load diffraction data");
    }
    // Initialize global density maps, masks and structureFactors
    initializeGlobalDensityMapsAndStructureFactors();
    displayResolutionShell();
    
    // -------------------------------------------------------------------------
    
    // 生成所有原点平移和反演的单胞pdb文件，对应当前空间群允许的原点选择以及反演结构
    //generateUnitCellForAllOrigTranAndInversion("./PDB_file/" + pdbCode + ".pdb");
    // 生成所有原点平移和反演的用于Coot展示pdb文
    //generateAllOrigTranAndInversionForCoot("./PDB_file/" + pdbCode + ".pdb");
    // 将衍射数据按照分辨率排序
    //string mtzFile = "./PDB_file/" + pdbCode + "-sf.mtz";
    //writeUniqReflectionsToTxtFile(mtzFile, "./PDB_file/" + pdbCode + "_unique_reflections.txt");
    
    // -------------------------------------------------------------------------
    
    // Load structure factors from mtz file for computing density histogram
    if (!pdbHist.empty()) {
        loadHistogramFatomData("./PDB_file/" + pdbHist + "_fatom.mtz");
        loadHistogramFmodelData("./PDB_file/" + pdbHist + "_fmodel.mtz");
    }
    computeHistogramBoundaries();  // reference histogram
    
    // Alternatively, a histogram txt file can be loaded directly
    // Load reference histogram from txt file
    //if (!pdbHist.empty()) {
    //    loadHistogramTxtFile("./PDB_file/" + pdbHist + "_histogram.txt");
    //}  // exit(EXIT_SUCCESS);
    
    // -------------------------------------------------------------------------
    
    // Load true phase data for validation
    loadReferenceModelData("./PDB_file/" + pdbCode + "_fmodel.mtz");
    
    generateReferenceMapAndMask();  // reference mask of the default origin choice
    generateReferenceMaskForAlternativeOrigins();  // true masks of alternative origin choices
    generateReferencePhaseForAlternativeOrigins(); // true phases of alternative origin choices
    /*
    // 输出不同原点选择方式对应的真实轮廓
    for(size_t i = 0; i < std::min(alternativeOrigins.size(), static_cast<size_t>(8)); ++i){
        fileNameStream.str(std::string());
        fileNameStream.clear();
        fileNameStream << "./PDB_file/" << pdbCode << "_mask_orig_" << i << ".map";
        writeDensityMaskUnitCellToMapFile(fileNameStream.str(), alternativeOriginMasks[i]);
        //fileNameStream << "./PDB_file/" << pdbCode << "_mask_orig_" << i << ".pdb";
        //writeDensityMaskUnitCellToPdbFile(fileNameStream.str(), alternativeOriginMasks[i]);
    }
    int spgNum = referenceDensity.spacegroup().spacegroup_number();
    if (isSpaceGroupAchiral(spgNum) or isSpaceGroupChiral(spgNum)) { // 43个非手性反演异构同空间群，11对22个手性反演异构不同空间群
        for(size_t i = 0; i < std::min(alternativeOriginsInverted.size(), static_cast<size_t>(8)); ++i){
            fileNameStream.str(std::string());
            fileNameStream.clear();
            fileNameStream << "./PDB_file/" << pdbCode << "_inverted_mask_orig_" << i << ".map";
            writeDensityMaskUnitCellToMapFile(fileNameStream.str(), alternativeOriginMasksInverted[i]);
            //fileNameStream << "./PDB_file/" << pdbCode << "_mask_orig_" << i << ".pdb";
            //writeDensityMaskUnitCellToPdbFile(fileNameStream.str(), alternativeOriginMasks[i]);
        }
    }
    */
    // -------------------------------------------------------------------------
    // fixed protein mask
    //for (auto ix=referenceMask.first(); !ix.last(); ix.next()) {
    //    proteinMask[ix] = referenceMask[ix];
    //}
    // -------------------------------------------------------------------------
    
    // Load pre-computed sigma weights if using true sigma
    if (useDataWeightSigma) {
        //方式一，计算每次迭代对应的权重参数observedSigmaArray
        std::cout << "Computing observed sigma for data weighting ..." << std::endl;
        computeSigmaObsForDataWeighting(measuredFAll);
        //方式二，加载计算好的数据
        //std::cout << "Loading pre-computed sigma weights" << std::endl;
        //std::string fileName = "./PDB_file/" + pdbCode + "_sigma_weight.txt";
        //loadSigmaWeightArray(fileName);
    }

    // -------------------------------------------------------------------------
    
    // Initialize starting density based on selected method
    if (startFromRandomMap) {
        std::cout << "Initializing from random density map" << std::endl;
        unsigned long seed = generateRandomSeed();
        srand(seed);
        
        // Save seed for reproducibility
        std::string fileName = "./output/" + std::to_string(mpiRank) + "_seed_init_dens_asu.txt";
        std::ofstream seedFile(fileName);
        seedFile << seed << std::endl;
        seedFile.close();
        
        initializeDensityFromRandomMap();
        
    } else if (startFromPDB) {
        std::cout << "Initializing from PDB structure" << std::endl;
        unsigned long seed = generateRandomSeed();
        srand(seed);
//        if (mpiRank<3) {
//            std::string fileName = "./PDB_file/" + pdbCode + "_orig_" + std::to_string(mpiRank) + ".pdb";
//            initializeDensityFromPDBStructure(fileName, 0.005f); // Use 2% of atoms
//        } else {
//            std::string fileName = "./PDB_file/" + pdbCode + "_orig_" + std::to_string(mpiRank) + "_inverted.pdb";
//            initializeDensityFromPDBStructure(fileName, 0.005f); // Use 2% of atoms
//        }
        std::string fileName = "./PDB_file/" + pdbCode + ".pdb";
        initializeDensityFromPDBStructure(fileName, 0.02f); // Use 2% of atoms
    } else if (startFromMR) {
        std::cout << "Initializing from molecular replacement" << std::endl;
        currentDensity.fft_from(referenceFactors);
        copyPhaseBToPhaseA(structureFactors, referenceFactors);
        
    } else if (startFromMtz) {
        initializeDensityFromMtz("./PDB_file/5_structure_factors.mtz"); // Use calculated mtz
                
    } else { // Start from random phase
        std::cout << "Initializing from observed data with random phase" << std::endl;
        unsigned long seed = generateRandomSeed();
        srand(seed);
        
        // Save seed for reproducibility
        fileNameStream.str(std::string());
        fileNameStream.clear();
        fileNameStream << "./output/" << mpiRank << "_seed_init_phase.txt";
        std::ofstream seedFile(fileNameStream.str());
        seedFile << seed << std::endl;
        seedFile.close();
        
        initializeDensityFromRandomPhase();
    }
    
    std::cout << "Data processing initialization complete" << std::endl;
    
    // -------------------------------------------------------------------------
    
    //checkGridPointsAsu();
    //checkUnitCellDensity(currentDensity);
    //writeAsuGridToPdbFile("asu_grid.pdb", currentDensity);
    //writeUnitCellGridToPdbFile("unitcell_grid.pdb", currentDensity);
    //writeDensityMapUnitCellToMapFile("initial_density.map", currentDensity);
    //writeDensityMapUnitCellToPdbFile("initial_density.pdb", currentDensity, 0.1);
    //writeFullStructureFactorsToCifFile("./PDB_file/" + pdbCode + "_full_structure_factors", referenceFactors);
    //在C++中，使用exit()会跳过栈解退和析构函数的执行，导致对象无法正确清理。
    //exit(EXIT_SUCCESS); //clipper库中的对象依赖于析构函数来释放资源，因此直接退出会导致泄漏。
}

void initializeGeneticAlgorithmParameters(GeneticAlgorithmParams& params) {
    params.populationSize = mpiSize;
    params.genomeLength = gridPointCountAsu;
    params.crossoverPointCount = DEFAULT_CROSSOVER_POINTS;
    params.crossoverGeneRatio = DEFAULT_CROSSOVER_RATIO;
    params.mutationGeneRatio = DEFAULT_MUTATION_RATIO;
    std::cout << "Initialized genetic algorithm parameters." << std::endl;
}

void executeHybridInputOutputLoop() {
    std::cout << "Starting HIO iteration loop" << std::endl;
    
    // Initialize genetic algorithm parameters
    GeneticAlgorithmParams gaParams;
    if (useGeneticAlgorithm) {
        initializeGeneticAlgorithmParameters(gaParams);
    }
    
    // HIO本质是实空间的密度修正算法。HIO Loop包括实空间做HIO密度修正，倒空间做傅立叶振幅修正。
    for (int iter = 0; iter < maxIterations; iter++) {
        // Calculate variable sigma weights
        currentSigmaWeight = calculateLinearSigmaWeight(iter);
        
        // Apply sigma weighting if enabled
        if (useDataWeightSigma) {
            observedSigmaWeight = calculateVariableSigmaWeightNonlinear(iter);
            applyVariableStructureFactorWeighting(observedSigmaWeight);
        }

        // Apply HIO constraints
        applyHybridInputOutput(iter, currentDensity);
        
        // Calculate R-factors
        calculateRfactorWithoutSolventFlattening(iter);
        
        // Calculate quality metrics, calculateTrueMaskCorrelationToFindOrigineChoice()中判定并赋值
        int originIndex;
        bool isInverted;
        maskCorrelationValues[iter] = calculateTrueMaskCorrelationToFindOrigineChoice(proteinMask, originIndex, isInverted);
        //phaseErrorValues[iter] = calculatePhaseError(structureFactors, originIndex, isInverted);
        float* const rowStart = phaseErrorValues.get() + iter * (RESOLUTION_BIN_COUNT + 1);
        const auto& rowData = calculatePhaseError(structureFactors, originIndex, isInverted);
        for (int bin = 0; bin <= RESOLUTION_BIN_COUNT; ++bin) {
            rowStart[bin] = static_cast<float>(rowData[bin]);
        }
        
        // Genetic algorithm processing
        if (useGeneticAlgorithm) {
            
            // Update convergence status
            if (iter % 10 == 0 && iter > 0 && eliteCount < mpiSize) {  // 频繁检测耗时，太不频繁又回错过Rwork突然降低的时段
                updateConvergenceFlagAndEliteCount(iter);  // 用到MPI_Gather()和MPI_Scatter()耗时
            }

            // Early convergence check
            if (eliteCount == mpiSize &&
                iter < maxIterations - densityLimitingIterations - solventFlatteningIterations) {
                std::cout << "All population members converged at iteration " << iter << std::endl;
                std::cout << "Skipping to final phase at iteration "
                          << (maxIterations - densityLimitingIterations - solventFlatteningIterations)
                          << std::endl;
                iter = maxIterations - densityLimitingIterations - solventFlatteningIterations;
            }
            
            // Apply genetic algorithm every 100 iterations
            if (iter % 100 == 99) {
                //std::cout << "Executing genetic algorithm at iteration " << iter << std::endl;
                executeGeneticAlgorithmGeneration(iter, gaParams, currentDensity);
            }
        }
                
        // Progress reporting
        if (iter % 10 == 0 || iter == maxIterations - 1) {
            std::cout << "Iteration: " << std::setw(5) << iter 
                      << " | R-work: " << std::setprecision(4) << RworkValues[iter]
                      << " | R-free: " << std::setprecision(4) << RfreeValues[iter]
                      << " | Phase Error: " << std::setprecision(1) << rowStart[RESOLUTION_BIN_COUNT] << "°"
                      << " | Mask IoU: " << std::setprecision(3) << maskCorrelationValues[iter];
            if (useGeneticAlgorithm) {
                std::cout << " | Elite: " << eliteCount << "/" << mpiSize;
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "HIO iteration loop completed" << std::endl;
}

void executeTransitionHybridInputOutputLoop() {
    std::cout << "Starting THIO iteration loop" << std::endl;
    
    // Initialize genetic algorithm parameters
    GeneticAlgorithmParams gaParams;
    if (useGeneticAlgorithm) {
        initializeGeneticAlgorithmParameters(gaParams);
    }
    
    // THIO本质是实空间的密度修正算法。THIO Loop包括实空间做THIO密度修正，倒空间做傅立叶振幅修正。
    for (int iter = 0; iter < maxIterations; iter++) {
        // Calculate variable sigma weights
        currentSigmaWeight = calculateLinearSigmaWeight(iter);
        
        // Apply sigma weighting if enabled
        if (useDataWeightSigma) {
            observedSigmaWeight = calculateVariableSigmaWeightNonlinear(iter);
            applyVariableStructureFactorWeighting(observedSigmaWeight);
        }

        // Apply THIO constraints
        applyTransitionHybridInputOutput(iter, currentDensity);
            
        // Calculate R-factors
        calculateRfactorWithoutSolventFlattening(iter);
        
        // Calculate quality metrics, calculateTrueMaskCorrelationToFindOrigineChoice()中判定并赋值
        int originIndex;
        bool isInverted;
        maskCorrelationValues[iter] = calculateTrueMaskCorrelationToFindOrigineChoice(proteinMask, originIndex, isInverted);
        //phaseErrorValues[iter] = calculatePhaseError(structureFactors, originIndex, isInverted);
        float* const rowStart = phaseErrorValues.get() + iter * (RESOLUTION_BIN_COUNT + 1);
        const auto& rowData = calculatePhaseError(structureFactors, originIndex, isInverted);
        for (int bin = 0; bin <= RESOLUTION_BIN_COUNT; ++bin) {
            rowStart[bin] = static_cast<float>(rowData[bin]);
        }
        
        // Genetic algorithm processing
        if (useGeneticAlgorithm) {
            
            // Update convergence status
            if (iter % 10 == 0 && iter > 0 && eliteCount < mpiSize) {  // 频繁检测耗时，太不频繁又回错过Rwork突然降低的时段
                updateConvergenceFlagAndEliteCount(iter);  // 用到MPI_Gather()和MPI_Scatter()耗时
            }

            // Early convergence check
            if (eliteCount == mpiSize &&
                iter < maxIterations - densityLimitingIterations - solventFlatteningIterations) {
                std::cout << "All population members converged at iteration " << iter << std::endl;
                std::cout << "Skipping to final phase at iteration "
                          << (maxIterations - densityLimitingIterations - solventFlatteningIterations)
                          << std::endl;
                iter = maxIterations - densityLimitingIterations - solventFlatteningIterations;
            }
            
            // Apply genetic algorithm every 100 iterations
            if (iter % 100 == 99) {
                //std::cout << "Executing genetic algorithm at iteration " << iter << std::endl;
                executeGeneticAlgorithmGeneration(iter, gaParams, currentDensity);
            }
        }
                
        // Progress reporting
        if (iter % 10 == 0 || iter == maxIterations - 1) {
            std::cout << "Iteration: " << std::setw(5) << iter
                      << " | R-work: " << std::setprecision(4) << RworkValues[iter]
                      << " | R-free: " << std::setprecision(4) << RfreeValues[iter]
                      << " | Phase Error: " << std::setprecision(1) << rowStart[RESOLUTION_BIN_COUNT] << "°"
                      << " | Mask IoU: " << std::setprecision(3) << maskCorrelationValues[iter];
            if (useGeneticAlgorithm) {
                std::cout << " | Elite: " << eliteCount << "/" << mpiSize;
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "THIO iteration loop completed" << std::endl;
}

void executeContinuousHybridInputOutputLoop() {
    std::cout << "Starting CHIO iteration loop" << std::endl;
    
    // Initialize genetic algorithm parameters
    GeneticAlgorithmParams gaParams;
    if (useGeneticAlgorithm) {
        initializeGeneticAlgorithmParameters(gaParams);
    }
    
    // CHIO本质是实空间的密度修正算法。CHIO Loop包括实空间做CHIO密度修正，倒空间做傅立叶振幅修正。
    for (int iter = 0; iter < maxIterations; iter++) {
        // Calculate variable sigma weights
        currentSigmaWeight = calculateLinearSigmaWeight(iter);
        
        // Apply sigma weighting if enabled
        if (useDataWeightSigma) {
            observedSigmaWeight = calculateVariableSigmaWeightNonlinear(iter);
            applyVariableStructureFactorWeighting(observedSigmaWeight);
        }

        // Apply CHIO constraints
        applyContinuousHybridInputOutput(iter, currentDensity);
            
        // Calculate R-factors
        calculateRfactorWithoutSolventFlattening(iter);
        
        // Calculate quality metrics, calculateTrueMaskCorrelationToFindOrigineChoice()中判定并赋值
        int originIndex;
        bool isInverted;
        maskCorrelationValues[iter] = calculateTrueMaskCorrelationToFindOrigineChoice(proteinMask, originIndex, isInverted);
        //phaseErrorValues[iter] = calculatePhaseError(structureFactors, originIndex, isInverted);
        float* const rowStart = phaseErrorValues.get() + iter * (RESOLUTION_BIN_COUNT + 1);
        const auto& rowData = calculatePhaseError(structureFactors, originIndex, isInverted);
        for (int bin = 0; bin <= RESOLUTION_BIN_COUNT; ++bin) {
            rowStart[bin] = static_cast<float>(rowData[bin]);
        }
        
        // Genetic algorithm processing
        if (useGeneticAlgorithm) {
            
            // Update convergence status
            if (iter % 10 == 0 && iter > 0 && eliteCount < mpiSize) {  // 频繁检测耗时，太不频繁又回错过Rwork突然降低的时段
                updateConvergenceFlagAndEliteCount(iter);  // 用到MPI_Gather()和MPI_Scatter()耗时
            }

            // Early convergence check
            if (eliteCount == mpiSize &&
                iter < maxIterations - densityLimitingIterations - solventFlatteningIterations) {
                std::cout << "All population members converged at iteration " << iter << std::endl;
                std::cout << "Skipping to final phase at iteration "
                          << (maxIterations - densityLimitingIterations - solventFlatteningIterations)
                          << std::endl;
                iter = maxIterations - densityLimitingIterations - solventFlatteningIterations;
            }
            
            // Apply genetic algorithm every 100 iterations
            if (iter % 100 == 99) {
                //std::cout << "Executing genetic algorithm at iteration " << iter << std::endl;
                executeGeneticAlgorithmGeneration(iter, gaParams, currentDensity);
            }
        }
                
        // Progress reporting
        if (iter % 10 == 0 || iter == maxIterations - 1) {
            std::cout << "Iteration: " << std::setw(5) << iter
                      << " | R-work: " << std::setprecision(4) << RworkValues[iter]
                      << " | R-free: " << std::setprecision(4) << RfreeValues[iter]
                      << " | Phase Error: " << std::setprecision(1) << rowStart[RESOLUTION_BIN_COUNT] << "°"
                      << " | Mask IoU: " << std::setprecision(3) << maskCorrelationValues[iter];
            if (useGeneticAlgorithm) {
                std::cout << " | Elite: " << eliteCount << "/" << mpiSize;
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "CHIO iteration loop completed" << std::endl;
}

void executeHybridProjectionReflectionLoop() {
    std::cout << "Starting HPR iteration loop" << std::endl;
    
    // Initialize genetic algorithm parameters
    GeneticAlgorithmParams gaParams;
    if (useGeneticAlgorithm) {
        initializeGeneticAlgorithmParameters(gaParams);
    }
    
    // HPR本质是实空间的密度修正算法。HPR Loop包括实空间做HPR密度修正，倒空间做傅立叶振幅修正。
    for (int iter = 0; iter < maxIterations; iter++) {
        // Calculate variable sigma weights
        currentSigmaWeight = calculateLinearSigmaWeight(iter);
        
        // Apply sigma weighting if enabled
        if (useDataWeightSigma) {
            observedSigmaWeight = calculateVariableSigmaWeightNonlinear(iter);
            applyVariableStructureFactorWeighting(observedSigmaWeight);
        }

        // Apply HPR constraints
        applyHybridProjectionReflection(iter, currentDensity);
            
        // Calculate R-factors
        calculateRfactorWithoutSolventFlattening(iter);
        
        // Calculate quality metrics, calculateTrueMaskCorrelationToFindOrigineChoice()中判定并赋值
        int originIndex;
        bool isInverted;
        maskCorrelationValues[iter] = calculateTrueMaskCorrelationToFindOrigineChoice(proteinMask, originIndex, isInverted);
        //phaseErrorValues[iter] = calculatePhaseError(structureFactors, originIndex, isInverted);
        float* const rowStart = phaseErrorValues.get() + iter * (RESOLUTION_BIN_COUNT + 1);
        const auto& rowData = calculatePhaseError(structureFactors, originIndex, isInverted);
        for (int bin = 0; bin <= RESOLUTION_BIN_COUNT; ++bin) {
            rowStart[bin] = static_cast<float>(rowData[bin]);
        }
        
        // Genetic algorithm processing
        if (useGeneticAlgorithm) {
            
            // Update convergence status
            if (iter % 10 == 0 && iter > 0 && eliteCount < mpiSize) {  // 频繁检测耗时，太不频繁又回错过Rwork突然降低的时段
                updateConvergenceFlagAndEliteCount(iter);  // 用到MPI_Gather()和MPI_Scatter()耗时
            }

            // Early convergence check
            if (eliteCount == mpiSize &&
                iter < maxIterations - densityLimitingIterations - solventFlatteningIterations) {
                std::cout << "All population members converged at iteration " << iter << std::endl;
                std::cout << "Skipping to final phase at iteration "
                          << (maxIterations - densityLimitingIterations - solventFlatteningIterations)
                          << std::endl;
                iter = maxIterations - densityLimitingIterations - solventFlatteningIterations;
            }
            
            // Apply genetic algorithm every 100 iterations
            if (iter % 100 == 99) {
                //std::cout << "Executing genetic algorithm at iteration " << iter << std::endl;
                executeGeneticAlgorithmGeneration(iter, gaParams, currentDensity);
            }
        }
                
        // Progress reporting
        if (iter % 10 == 0 || iter == maxIterations - 1) {
            std::cout << "Iteration: " << std::setw(5) << iter
                      << " | R-work: " << std::setprecision(4) << RworkValues[iter]
                      << " | R-free: " << std::setprecision(4) << RfreeValues[iter]
                      << " | Phase Error: " << std::setprecision(1) << rowStart[RESOLUTION_BIN_COUNT] << "°"
                      << " | Mask IoU: " << std::setprecision(3) << maskCorrelationValues[iter];
            if (useGeneticAlgorithm) {
                std::cout << " | Elite: " << eliteCount << "/" << mpiSize;
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "HPR iteration loop completed" << std::endl;
}

void executeModifiedContinuousHybridInputOutputLoop() {
    std::cout << "Starting MCHIO iteration loop" << std::endl;
    
    // Initialize genetic algorithm parameters
    GeneticAlgorithmParams gaParams;
    if (useGeneticAlgorithm) {
        initializeGeneticAlgorithmParameters(gaParams);
    }
    
    // MCHIO本质是实空间的密度修正算法。MCHIO Loop包括实空间做MCHIO密度修正，倒空间做傅立叶振幅修正。
    for (int iter = 0; iter < maxIterations; iter++) {
        // Calculate variable sigma weights
        currentSigmaWeight = calculateLinearSigmaWeight(iter);
        
        // Apply sigma weighting if enabled
        if (useDataWeightSigma) {
            observedSigmaWeight = calculateVariableSigmaWeightNonlinear(iter);
            applyVariableStructureFactorWeighting(observedSigmaWeight);
        }

        // Apply MCHIO constraints
        applyModifiedContinuousHybridInputOutput(iter, currentDensity);
            
        // Calculate R-factors
        calculateRfactorWithoutSolventFlattening(iter);
        
        // Calculate quality metrics, calculateTrueMaskCorrelationToFindOrigineChoice()中判定并赋值
        int originIndex;
        bool isInverted;
        maskCorrelationValues[iter] = calculateTrueMaskCorrelationToFindOrigineChoice(proteinMask, originIndex, isInverted);
        //phaseErrorValues[iter] = calculatePhaseError(structureFactors, originIndex, isInverted);
        float* const rowStart = phaseErrorValues.get() + iter * (RESOLUTION_BIN_COUNT + 1);
        const auto& rowData = calculatePhaseError(structureFactors, originIndex, isInverted);
        for (int bin = 0; bin <= RESOLUTION_BIN_COUNT; ++bin) {
            rowStart[bin] = static_cast<float>(rowData[bin]);
        }
        
        // Genetic algorithm processing
        if (useGeneticAlgorithm) {
            
            // Update convergence status
            if (iter % 10 == 0 && iter > 0 && eliteCount < mpiSize) {  // 频繁检测耗时，太不频繁又回错过Rwork突然降低的时段
                updateConvergenceFlagAndEliteCount(iter);  // 用到MPI_Gather()和MPI_Scatter()耗时
            }

            // Early convergence check
            if (eliteCount == mpiSize &&
                iter < maxIterations - densityLimitingIterations - solventFlatteningIterations) {
                std::cout << "All population members converged at iteration " << iter << std::endl;
                std::cout << "Skipping to final phase at iteration "
                          << (maxIterations - densityLimitingIterations - solventFlatteningIterations)
                          << std::endl;
                iter = maxIterations - densityLimitingIterations - solventFlatteningIterations;
            }
            
            // Apply genetic algorithm every 100 iterations
            if (iter % 100 == 99) {
                //std::cout << "Executing genetic algorithm at iteration " << iter << std::endl;
                executeGeneticAlgorithmGeneration(iter, gaParams, currentDensity);
            }
        }
                
        // Progress reporting
        if (iter % 10 == 0 || iter == maxIterations - 1) {
            std::cout << "Iteration: " << std::setw(5) << iter
                      << " | R-work: " << std::setprecision(4) << RworkValues[iter]
                      << " | R-free: " << std::setprecision(4) << RfreeValues[iter]
                      << " | Phase Error: " << std::setprecision(1) << rowStart[RESOLUTION_BIN_COUNT] << "°"
                      << " | Mask IoU: " << std::setprecision(3) << maskCorrelationValues[iter];
            if (useGeneticAlgorithm) {
                std::cout << " | Elite: " << eliteCount << "/" << mpiSize;
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "MCHIO iteration loop completed" << std::endl;
}

void executeDifferenceMapLoop() {
    std::cout << "Starting Difference Map iteration loop" << std::endl;
    
    // Initialize genetic algorithm parameters
    GeneticAlgorithmParams gaParams;
    if (useGeneticAlgorithm) {
        initializeGeneticAlgorithmParameters(gaParams);
    }
    
    // DM本质是实空间的密度修正算法。DM Loop包括实空间做DM密度修正，倒空间做傅立叶振幅修正。
    for (int iter = 0; iter < maxIterations; iter++) {
        // Calculate variable sigma weights
        currentSigmaWeight = calculateLinearSigmaWeight(iter);
        
        // Apply sigma weighting if enabled
        if (useDataWeightSigma) {
            observedSigmaWeight = calculateVariableSigmaWeightNonlinear(iter);
            applyVariableStructureFactorWeighting(observedSigmaWeight);
        }
        
        // Apply difference map algorithm
        applyDifferenceMap(iter, currentDensity);
                    
        // Calculate R-factors
        calculateRfactorWithoutSolventFlattening(iter);
        
        // Calculate quality metrics, calculateTrueMaskCorrelationToFindOrigineChoice()中判定并赋值
        int originIndex;
        bool isInverted;
        maskCorrelationValues[iter] = calculateTrueMaskCorrelationToFindOrigineChoice(proteinMask, originIndex, isInverted);
        //phaseErrorValues[iter] = calculatePhaseError(structureFactors, originIndex, isInverted);
        float* const rowStart = phaseErrorValues.get() + iter * (RESOLUTION_BIN_COUNT + 1);
        const auto& rowData = calculatePhaseError(structureFactors, originIndex, isInverted);
        for (int bin = 0; bin <= RESOLUTION_BIN_COUNT; ++bin) {
            rowStart[bin] = static_cast<float>(rowData[bin]);
        }
        
        // Genetic algorithm processing
        if (useGeneticAlgorithm) {
            
            // Update convergence status
            if (iter % 10 == 0 && iter > 0 && eliteCount < mpiSize) {  // 频繁检测耗时，太不频繁又回错过Rwork突然降低的时段
                updateConvergenceFlagAndEliteCount(iter);  // 用到MPI_Gather()和MPI_Scatter()耗时
            }

            // Early convergence check
            if (eliteCount == mpiSize &&
                iter < maxIterations - densityLimitingIterations - solventFlatteningIterations) {
                std::cout << "All population members converged at iteration " << iter << std::endl;
                std::cout << "Skipping to final phase at iteration "
                          << (maxIterations - densityLimitingIterations - solventFlatteningIterations)
                          << std::endl;
                iter = maxIterations - densityLimitingIterations - solventFlatteningIterations;
            }
            
            // Apply genetic algorithm every 100 iterations
            if (iter % 100 == 99) {
                //std::cout << "Executing genetic algorithm at iteration " << iter << std::endl;
                executeGeneticAlgorithmGeneration(iter, gaParams, currentDensity);
            }
        }
        
        // Progress reporting
        if (iter % 10 == 0 || iter == maxIterations - 1) {
            std::cout << "Iteration: " << std::setw(5) << iter 
                      << " | R-work: " << std::setprecision(4) << RworkValues[iter]
                      << " | R-free: " << std::setprecision(4) << RfreeValues[iter]
                      << " | Phase Error: " << std::setprecision(1) << rowStart[RESOLUTION_BIN_COUNT] << "°"
                      << " | Mask IoU: " << std::setprecision(3) << maskCorrelationValues[iter];
            if (useGeneticAlgorithm) {
                std::cout << " | Elite: " << eliteCount << "/" << mpiSize;
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "Difference Map iteration loop completed" << std::endl;
}

void executeModifiedDifferenceMapLoop() {
    std::cout << "Starting Modified Difference Map iteration loop" << std::endl;
    
    // Initialize genetic algorithm parameters
    GeneticAlgorithmParams gaParams;
    if (useGeneticAlgorithm) {
        initializeGeneticAlgorithmParameters(gaParams);
    }
    
    // MDM本质是实空间的密度修正算法。MDM Loop包括实空间做DM密度修正，倒空间做傅立叶振幅修正。
    for (int iter = 0; iter < maxIterations; iter++) {
        // Calculate variable sigma weights
        currentSigmaWeight = calculateLinearSigmaWeight(iter);
        
        // Apply sigma weighting if enabled
        if (useDataWeightSigma) {
            observedSigmaWeight = calculateVariableSigmaWeightNonlinear(iter);
            applyVariableStructureFactorWeighting(observedSigmaWeight);
        }
        
        // Apply modified difference map algorithm
        applyModifiedDifferenceMap(iter, currentDensity);
                    
        // Calculate R-factors
        calculateRfactorWithoutSolventFlattening(iter);
        
        // Calculate quality metrics, calculateTrueMaskCorrelationToFindOrigineChoice()中判定并赋值
        int originIndex;
        bool isInverted;
        maskCorrelationValues[iter] = calculateTrueMaskCorrelationToFindOrigineChoice(proteinMask, originIndex, isInverted);
        //phaseErrorValues[iter] = calculatePhaseError(structureFactors, originIndex, isInverted);
        float* const rowStart = phaseErrorValues.get() + iter * (RESOLUTION_BIN_COUNT + 1);
        const auto& rowData = calculatePhaseError(structureFactors, originIndex, isInverted);
        for (int bin = 0; bin <= RESOLUTION_BIN_COUNT; ++bin) {
            rowStart[bin] = static_cast<float>(rowData[bin]);
        }
        
        // Genetic algorithm processing
        if (useGeneticAlgorithm) {
            
            // Update convergence status
            if (iter % 10 == 0 && iter > 0 && eliteCount < mpiSize) {  // 频繁检测耗时，太不频繁又回错过Rwork突然降低的时段
                updateConvergenceFlagAndEliteCount(iter);  // 用到MPI_Gather()和MPI_Scatter()耗时
            }

            // Early convergence check
            if (eliteCount == mpiSize &&
                iter < maxIterations - densityLimitingIterations - solventFlatteningIterations) {
                std::cout << "All population members converged at iteration " << iter << std::endl;
                std::cout << "Skipping to final phase at iteration "
                          << (maxIterations - densityLimitingIterations - solventFlatteningIterations)
                          << std::endl;
                iter = maxIterations - densityLimitingIterations - solventFlatteningIterations;
            }
            
            // Apply genetic algorithm every 100 iterations
            if (iter % 100 == 99) {
                //std::cout << "Executing genetic algorithm at iteration " << iter << std::endl;
                executeGeneticAlgorithmGeneration(iter, gaParams, currentDensity);
            }
        }
        
        // Progress reporting
        if (iter % 10 == 0 || iter == maxIterations - 1) {
            std::cout << "Iteration: " << std::setw(5) << iter
                      << " | R-work: " << std::setprecision(4) << RworkValues[iter]
                      << " | R-free: " << std::setprecision(4) << RfreeValues[iter]
                      << " | Phase Error: " << std::setprecision(1) << rowStart[RESOLUTION_BIN_COUNT] << "°"
                      << " | Mask IoU: " << std::setprecision(3) << maskCorrelationValues[iter];
            if (useGeneticAlgorithm) {
                std::cout << " | Elite: " << eliteCount << "/" << mpiSize;
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "Modified Difference Map iteration loop completed" << std::endl;
}

void executeAveragedSuccessiveReflectionsLoop() {
    std::cout << "Starting Averaged Successive Reflections iteration loop" << std::endl;
    
    // Initialize genetic algorithm parameters
    GeneticAlgorithmParams gaParams;
    if (useGeneticAlgorithm) {
        initializeGeneticAlgorithmParameters(gaParams);
    }
    
    // ASR本质是实空间的密度修正算法。ASR Loop包括实空间做DM密度修正，倒空间做傅立叶振幅修正。
    for (int iter = 0; iter < maxIterations; iter++) {
        // Calculate variable sigma weights
        currentSigmaWeight = calculateLinearSigmaWeight(iter);
        
        // Apply sigma weighting if enabled
        if (useDataWeightSigma) {
            observedSigmaWeight = calculateVariableSigmaWeightNonlinear(iter);
            applyVariableStructureFactorWeighting(observedSigmaWeight);
        }
        
        // Apply averaged successive reflections algorithm
        applyAveragedSuccessiveReflections(iter, currentDensity);
                    
        // Calculate R-factors
        calculateRfactorWithoutSolventFlattening(iter);
        
        // Calculate quality metrics, calculateTrueMaskCorrelationToFindOrigineChoice()中判定并赋值
        int originIndex;
        bool isInverted;
        maskCorrelationValues[iter] = calculateTrueMaskCorrelationToFindOrigineChoice(proteinMask, originIndex, isInverted);
        //phaseErrorValues[iter] = calculatePhaseError(structureFactors, originIndex, isInverted);
        float* const rowStart = phaseErrorValues.get() + iter * (RESOLUTION_BIN_COUNT + 1);
        const auto& rowData = calculatePhaseError(structureFactors, originIndex, isInverted);
        for (int bin = 0; bin <= RESOLUTION_BIN_COUNT; ++bin) {
            rowStart[bin] = static_cast<float>(rowData[bin]);
        }
        
        // Genetic algorithm processing
        if (useGeneticAlgorithm) {
            
            // Update convergence status
            if (iter % 10 == 0 && iter > 0 && eliteCount < mpiSize) {  // 频繁检测耗时，太不频繁又回错过Rwork突然降低的时段
                updateConvergenceFlagAndEliteCount(iter);  // 用到MPI_Gather()和MPI_Scatter()耗时
            }

            // Early convergence check
            if (eliteCount == mpiSize &&
                iter < maxIterations - densityLimitingIterations - solventFlatteningIterations) {
                std::cout << "All population members converged at iteration " << iter << std::endl;
                std::cout << "Skipping to final phase at iteration "
                          << (maxIterations - densityLimitingIterations - solventFlatteningIterations)
                          << std::endl;
                iter = maxIterations - densityLimitingIterations - solventFlatteningIterations;
            }
            
            // Apply genetic algorithm every 100 iterations
            if (iter % 100 == 99) {
                //std::cout << "Executing genetic algorithm at iteration " << iter << std::endl;
                executeGeneticAlgorithmGeneration(iter, gaParams, currentDensity);
            }
        }
        
        // Progress reporting
        if (iter % 10 == 0 || iter == maxIterations - 1) {
            std::cout << "Iteration: " << std::setw(5) << iter
                      << " | R-work: " << std::setprecision(4) << RworkValues[iter]
                      << " | R-free: " << std::setprecision(4) << RfreeValues[iter]
                      << " | Phase Error: " << std::setprecision(1) << rowStart[RESOLUTION_BIN_COUNT] << "°"
                      << " | Mask IoU: " << std::setprecision(3) << maskCorrelationValues[iter];
            if (useGeneticAlgorithm) {
                std::cout << " | Elite: " << eliteCount << "/" << mpiSize;
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "Averaged Sucessive Reflections iteration loop completed" << std::endl;
}

void executeRelaxedAveragedAlternatingReflectionsLoop() {
    std::cout << "Starting Relaxed Averaged Alternating Reflections iteration loop" << std::endl;
    
    // Initialize genetic algorithm parameters
    GeneticAlgorithmParams gaParams;
    if (useGeneticAlgorithm) {
        initializeGeneticAlgorithmParameters(gaParams);
    }
    
    // RAAR本质是实空间的密度修正算法。RAAR Loop包括实空间做DM密度修正，倒空间做傅立叶振幅修正。
    for (int iter = 0; iter < maxIterations; iter++) {
        // Calculate variable sigma weights
        currentSigmaWeight = calculateLinearSigmaWeight(iter);
        
        // Apply sigma weighting if enabled
        if (useDataWeightSigma) {
            observedSigmaWeight = calculateVariableSigmaWeightNonlinear(iter);
            applyVariableStructureFactorWeighting(observedSigmaWeight);
        }
        
        // Apply relaxed averaged alternative reflections algorithm
        applyRelaxedAveragedAlternatingReflections(iter, currentDensity);
                    
        // Calculate R-factors
        calculateRfactorWithoutSolventFlattening(iter);
        
        // Calculate quality metrics, calculateTrueMaskCorrelationToFindOrigineChoice()中判定并赋值
        int originIndex;
        bool isInverted;
        maskCorrelationValues[iter] = calculateTrueMaskCorrelationToFindOrigineChoice(proteinMask, originIndex, isInverted);
        //phaseErrorValues[iter] = calculatePhaseError(structureFactors, originIndex, isInverted);
        float* const rowStart = phaseErrorValues.get() + iter * (RESOLUTION_BIN_COUNT + 1);
        const auto& rowData = calculatePhaseError(structureFactors, originIndex, isInverted);
        for (int bin = 0; bin <= RESOLUTION_BIN_COUNT; ++bin) {
            rowStart[bin] = static_cast<float>(rowData[bin]);
        }
        
        // Genetic algorithm processing
        if (useGeneticAlgorithm) {
            
            // Update convergence status
            if (iter % 10 == 0 && iter > 0 && eliteCount < mpiSize) {  // 频繁检测耗时，太不频繁又回错过Rwork突然降低的时段
                updateConvergenceFlagAndEliteCount(iter);  // 用到MPI_Gather()和MPI_Scatter()耗时
            }

            // Early convergence check
            if (eliteCount == mpiSize &&
                iter < maxIterations - densityLimitingIterations - solventFlatteningIterations) {
                std::cout << "All population members converged at iteration " << iter << std::endl;
                std::cout << "Skipping to final phase at iteration "
                          << (maxIterations - densityLimitingIterations - solventFlatteningIterations)
                          << std::endl;
                iter = maxIterations - densityLimitingIterations - solventFlatteningIterations;
            }
            
            // Apply genetic algorithm every 100 iterations
            if (iter % 100 == 99) {
                //std::cout << "Executing genetic algorithm at iteration " << iter << std::endl;
                executeGeneticAlgorithmGeneration(iter, gaParams, currentDensity);
            }
        }
        
        // Progress reporting
        if (iter % 10 == 0 || iter == maxIterations - 1) {
            std::cout << "Iteration: " << std::setw(5) << iter
                      << " | R-work: " << std::setprecision(4) << RworkValues[iter]
                      << " | R-free: " << std::setprecision(4) << RfreeValues[iter]
                      << " | Phase Error: " << std::setprecision(1) << rowStart[RESOLUTION_BIN_COUNT] << "°"
                      << " | Mask IoU: " << std::setprecision(3) << maskCorrelationValues[iter];
            if (useGeneticAlgorithm) {
                std::cout << " | Elite: " << eliteCount << "/" << mpiSize;
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "Relaxed Averaged Alternating Reflections iteration loop completed" << std::endl;
}

void executeModifiedRelaxedAveragedAlternatingReflectionsLoop() {
    std::cout << "Starting Modified Relaxed Averaged Alternating Reflections iteration loop" << std::endl;
    
    // Initialize genetic algorithm parameters
    GeneticAlgorithmParams gaParams;
    if (useGeneticAlgorithm) {
        initializeGeneticAlgorithmParameters(gaParams);
    }
    
    // MRAAR本质是实空间的密度修正算法。MRAAR Loop包括实空间做DM密度修正，倒空间做傅立叶振幅修正。
    for (int iter = 0; iter < maxIterations; iter++) {
        // Calculate variable sigma weights
        currentSigmaWeight = calculateLinearSigmaWeight(iter);
        
        // Apply sigma weighting if enabled
        if (useDataWeightSigma) {
            observedSigmaWeight = calculateVariableSigmaWeightNonlinear(iter);
            applyVariableStructureFactorWeighting(observedSigmaWeight);
        }
        
        // Apply modified relaxed averaged alternative reflections algorithm
        applyModifiedRelaxedAveragedAlternatingReflections(iter, currentDensity);
                    
        // Calculate R-factors
        calculateRfactorWithoutSolventFlattening(iter);
        
        // Calculate quality metrics, calculateTrueMaskCorrelationToFindOrigineChoice()中判定并赋值
        int originIndex;
        bool isInverted;
        maskCorrelationValues[iter] = calculateTrueMaskCorrelationToFindOrigineChoice(proteinMask, originIndex, isInverted);
        //phaseErrorValues[iter] = calculatePhaseError(structureFactors, originIndex, isInverted);
        float* const rowStart = phaseErrorValues.get() + iter * (RESOLUTION_BIN_COUNT + 1);
        const auto& rowData = calculatePhaseError(structureFactors, originIndex, isInverted);
        for (int bin = 0; bin <= RESOLUTION_BIN_COUNT; ++bin) {
            rowStart[bin] = static_cast<float>(rowData[bin]);
        }
        
        // Genetic algorithm processing
        if (useGeneticAlgorithm) {
            
            // Update convergence status
            if (iter % 10 == 0 && iter > 0 && eliteCount < mpiSize) {  // 频繁检测耗时，太不频繁又回错过Rwork突然降低的时段
                updateConvergenceFlagAndEliteCount(iter);  // 用到MPI_Gather()和MPI_Scatter()耗时
            }

            // Early convergence check
            if (eliteCount == mpiSize &&
                iter < maxIterations - densityLimitingIterations - solventFlatteningIterations) {
                std::cout << "All population members converged at iteration " << iter << std::endl;
                std::cout << "Skipping to final phase at iteration "
                          << (maxIterations - densityLimitingIterations - solventFlatteningIterations)
                          << std::endl;
                iter = maxIterations - densityLimitingIterations - solventFlatteningIterations;
            }
            
            // Apply genetic algorithm every 100 iterations
            if (iter % 100 == 99) {
                //std::cout << "Executing genetic algorithm at iteration " << iter << std::endl;
                executeGeneticAlgorithmGeneration(iter, gaParams, currentDensity);
            }
        }
        
        // Progress reporting
        if (iter % 10 == 0 || iter == maxIterations - 1) {
            std::cout << "Iteration: " << std::setw(5) << iter
                      << " | R-work: " << std::setprecision(4) << RworkValues[iter]
                      << " | R-free: " << std::setprecision(4) << RfreeValues[iter]
                      << " | Phase Error: " << std::setprecision(1) << rowStart[RESOLUTION_BIN_COUNT] << "°"
                      << " | Mask IoU: " << std::setprecision(3) << maskCorrelationValues[iter];
            if (useGeneticAlgorithm) {
                std::cout << " | Elite: " << eliteCount << "/" << mpiSize;
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "Modified Relaxed Averaged Alternating Reflections iteration loop completed" << std::endl;
}

void executeHybridDifferenceMapLoop(){
    
    // Initialize genetic algorithm parameters
    GeneticAlgorithmParams gaParams;
    if (useGeneticAlgorithm) {
        initializeGeneticAlgorithmParameters(gaParams);
    }
    
    std::cout << "Starting Hybrid Difference Map iteration loop" << std::endl;
    // HDM本质是实空间的密度修正算法。HDM Loop包括实空间做HDM密度修正，倒空间做傅立叶振幅修正。
    for (int iter = 0; iter < maxIterations; iter++) {
        // Calculate variable sigma weights
        currentSigmaWeight = calculateLinearSigmaWeight(iter);

        // Apply sigma weighting if enabled
        if (useDataWeightSigma) {
            observedSigmaWeight = calculateVariableSigmaWeightNonlinear(iter);
            applyVariableStructureFactorWeighting(observedSigmaWeight);
        }
       
        // Apply new difference map algorithm
        applyHybridDifferenceMap(iter, currentDensity);
        
        // Calculate R-factors
        calculateRfactorWithoutSolventFlattening(iter);
        
        // Calculate quality metrics, calculateTrueMaskCorrelationToFindOrigineChoice()中判定并赋值
        int originIndex;
        bool isInverted;
        maskCorrelationValues[iter] = calculateTrueMaskCorrelationToFindOrigineChoice(proteinMask, originIndex, isInverted);
        //phaseErrorValues[iter] = calculatePhaseError(structureFactors, originIndex, isInverted);
        float* const rowStart = phaseErrorValues.get() + iter * (RESOLUTION_BIN_COUNT + 1);
        const auto& rowData = calculatePhaseError(structureFactors, originIndex, isInverted);
        for (int bin = 0; bin <= RESOLUTION_BIN_COUNT; ++bin) {
            rowStart[bin] = static_cast<float>(rowData[bin]);
        }
                        
        // Genetic algorithm processing
        if (useGeneticAlgorithm) {
            
            // Update convergence status
            if (iter % 10 == 0 && iter > 0 && eliteCount < mpiSize) {  // 频繁检测耗时，太不频繁又回错过Rwork突然降低的时段
                updateConvergenceFlagAndEliteCount(iter);  // 用到MPI_Gather()和MPI_Scatter()耗时
            }
            
            // Early convergence check
            if (eliteCount == mpiSize &&
                iter < maxIterations - densityLimitingIterations - solventFlatteningIterations) {
                std::cout << "All population members converged at iteration " << iter << std::endl;
                std::cout << "Skipping to final phase at iteration "
                          << (maxIterations - densityLimitingIterations - solventFlatteningIterations)
                          << std::endl;
                iter = maxIterations - densityLimitingIterations - solventFlatteningIterations;
            }
            
            // Apply genetic algorithm every 100 iterations
            if (iter % 100 == 99) {
                //std::cout << "Executing genetic algorithm at iteration " << iter << std::endl;
                executeGeneticAlgorithmGeneration(iter, gaParams, currentDensity);
            }
        }
        
        // Progress reporting
        if (iter % 10 == 0 || iter == maxIterations - 1) {
            std::cout << "Iteration: " << std::setw(5) << iter
                      << " | R-work: " << std::setprecision(4) << RworkValues[iter]
                      << " | R-free: " << std::setprecision(4) << RfreeValues[iter]
                      << " | Phase Error: " << std::setprecision(1) << rowStart[RESOLUTION_BIN_COUNT] << "°"
                      << " | Mask IoU: " << std::setprecision(3) << maskCorrelationValues[iter];
            if (useGeneticAlgorithm) {
                std::cout << " | Elite: " << eliteCount << "/" << mpiSize;
            }
            std::cout << std::endl;
        }
    }
    
    std::cout << "Hybrid Difference Map iteration loop completed" << std::endl;
}

void generateOutputFiles() {
    std::cout << "Generating output files" << std::endl;
    
//    // Write density values
//    fileNameStream.str(std::string());
//    fileNameStream.clear();
//    fileNameStream << "./output/" << mpiRank << "_protein_density_absolute.txt";
//    writeFloatArrayToTxtFile(fileNameStream.str(), proteinDensityAbsolute.get());
//
//    fileNameStream.str(std::string());
//    fileNameStream.clear();
//    fileNameStream << "./output/" << mpiRank << "_solvent_density_absolute.txt";
//    writeFloatArrayToTxtFile(fileNameStream.str(), solventDensityAbsolute.get());
    
    // Write density feedback values
    fileNameStream.str(std::string());
    fileNameStream.clear();
    fileNameStream << "./output/" << mpiRank << "_protein_density_deviation.txt";
    writeFloatArrayToTxtFile(fileNameStream.str(), proteinDensityDeviation.get());
    
    fileNameStream.str(std::string());
    fileNameStream.clear();
    fileNameStream << "./output/" << mpiRank << "_solvent_density_deviation.txt";
    writeFloatArrayToTxtFile(fileNameStream.str(), solventDensityDeviation.get());
    
    // Write R-factors
    fileNameStream.str(std::string());
    fileNameStream.clear();
    fileNameStream << "./output/" << mpiRank << "_R_work.txt";
    writeFloatArrayToTxtFile(fileNameStream.str(), RworkValues.get());

    fileNameStream.str(std::string());
    fileNameStream.clear();
    fileNameStream << "./output/" << mpiRank << "_R_free.txt";
    writeFloatArrayToTxtFile(fileNameStream.str(), RfreeValues.get());

    // Write phase errors and mask correlations
    fileNameStream.str(std::string());
    fileNameStream.clear();
    fileNameStream << "./output/" << mpiRank << "_mean_phase_error.txt";
    writeFloat2DArrayToFile(fileNameStream.str(), phaseErrorValues.get());

    // Write mask correlations
    fileNameStream.str(std::string());
    fileNameStream.clear();
    fileNameStream << "./output/" << mpiRank << "_mask_correlation.txt";
    writeFloatArrayToTxtFile(fileNameStream.str(), maskCorrelationValues.get());

    // Write structure factors
    fileNameStream.str(std::string());
    fileNameStream.clear();
    fileNameStream << "./output/" << mpiRank << "_structure_factors.mtz";
    writeStructureFactorsToMtzFile(fileNameStream.str(), currentDensity);

    // Write CIF file for visualization
    fileNameStream.str(std::string());
    fileNameStream.clear();
    fileNameStream << "./output/" << mpiRank << "_phase.cif";
    writeCIFFileForVisualization(fileNameStream.str(), structureFactors);
    
    // Write density file for visualization
    fileNameStream.str(std::string());
    fileNameStream.clear();
    fileNameStream << "./output/" << mpiRank << "_density_asu.map";
    writeDensityMapAsuToMapFile(fileNameStream.str(), currentDensity);
    
    // Write density file for visualization
    fileNameStream.str(std::string());
    fileNameStream.clear();
    fileNameStream << "./output/" << mpiRank << "_density_unit_cell.map";
    writeDensityMapUnitCellToMapFile(fileNameStream.str(), currentDensity);

    // Write density mask file for visualization
    fileNameStream.str(std::string());
    fileNameStream.clear();
    fileNameStream << "./output/" << mpiRank << "_mask_asu.map";
    writeDensityMaskAsuToMapFile(fileNameStream.str(), proteinMask);
    
    // Write density mask file for visualization
    fileNameStream.str(std::string());
    fileNameStream.clear();
    fileNameStream << "./output/" << mpiRank << "_mask_unit_cell.map";
    writeDensityMaskUnitCellToMapFile(fileNameStream.str(), proteinMask);
    
    if (useGeneticAlgorithm) {
        if (mpiRank == 0) {  // 在rank0上计算的平均密度
            // Write phase errors of average elite density
            fileNameStream.str(std::string());
            fileNameStream.clear();
            fileNameStream << "./output/" << mpiRank << "_mean_phase_error_average_elite.txt";
            writeFloat2DArrayToFile(fileNameStream.str(), phaseErrorValuesAverageElite.get());
        }
    }
    std::cout << "Output files generated." << std::endl;
}

int main(int argc, char** argv) {
    try {
        // Initialize MPI
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

        std::cout << "Starting crystallographic direct phasing (Rank " << mpiRank 
                  << "/" << mpiSize << ")\n" << std::endl;
        if (mpiSize == 1) {
            useGeneticAlgorithm = false;
        }
            
        // Data processing and initialization
        initializeDataProcessing();

        clock_t start = clock();
                
        // Execute main iteration loop (default: Hybrid-Input Output)
        if (useHybridInputOutput) {
            executeHybridInputOutputLoop();
        } else if (useTransitionHybridInputOutput) {
            executeTransitionHybridInputOutputLoop();
        } else if (useContinuousHybridInputOutput) {
            executeContinuousHybridInputOutputLoop();
        } else if (useHybridProjectionReflection) {
            executeHybridProjectionReflectionLoop();
        } else if (useModifiedContinuousHybridInputOutput) {
            executeModifiedContinuousHybridInputOutputLoop();
        } else if (useDifferenceMap) {
            executeDifferenceMapLoop();
        } else if (useModifiedDifferenceMap) {
            executeModifiedDifferenceMapLoop();
        } else if (useAveragedSuccessiveReflections) {
            executeAveragedSuccessiveReflectionsLoop();
        } else if (useRelaxedAveragedAlternatingReflections) {
            executeRelaxedAveragedAlternatingReflectionsLoop();
        } else if (useModifiedRelaxedAveragedAlternatingReflections) {
            executeModifiedRelaxedAveragedAlternatingReflectionsLoop();
        } else if (useHybridDifferenceMap) {
            executeHybridDifferenceMapLoop();
        }

        clock_t end = clock();
        double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        std::cout << "程序运行耗时: " << duration << "秒" << std::endl;
        // Generate output files
        generateOutputFiles();

        // Cleanup and finalization
        deallocateArrays();
        
        std::cout << "Processing completed." << std::endl;
        
        // Finalize MPI
        MPI_Finalize();
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error in process " << mpiRank << ": " << e.what() << std::endl;
        
        // Clean up on error
        deallocateArrays();
        MPI_Finalize();
        
        return 1;
    }
}
