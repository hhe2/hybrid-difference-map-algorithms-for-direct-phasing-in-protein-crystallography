#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "globalvars.h"

// Map averaging filter for density smoothing
class MapAverageFilter : public clipper::MapFilterFn_base {
public:
    explicit MapAverageFilter(const clipper::ftype& sigma) : sigma_(sigma) {}
    ~MapAverageFilter() {}
    
    clipper::ftype operator()(const clipper::ftype& radius) const {
        return exp(-0.5 * pow(radius / sigma_, 2));
    }
private:
    clipper::ftype sigma_;
};

// Density statistics
struct DensityStats {
    float min;
    float max;
    float average;
};
DensityStats calculateDensityMinMax(const clipper::Xmap<float>& localDensity);

// ===== Sigma weighting functions =====
float calculateVariableSigmaWeightNonlinear(const int iteration);
void applyVariableStructureFactorWeighting(float sigmaWeight);

// ===== Data replacement functions =====
void fillAndReplaceObservedData(const int iteration, clipper::HKL_data<clipper::data64::F_phi>& fphi);

// ===== Mask generation functions =====
float calculateLinearSigmaWeight(const int iteration);
float sortAndFindThresholdDensity(const clipper::Xmap<float>& density,
                           const float solventFraction);
float findThresholdDensity(const clipper::Xmap<float>& density,
                           const float solventFraction);
float findThresholdDensityInsideMask(const clipper::Xmap<float>& density,
                                     const clipper::Xmap<int>& mask,
                                     const float shrinkFraction);
void computeWeightedAverageDensity(const clipper::Xmap<float>& density,
                                   clipper::Xmap<float>& averagedDensity,
                                   const float sigmaWeight);
void generateProteinSolventMask(const clipper::Xmap<float>& density,
                                clipper::Xmap<int>& mask,
                                const float sigmaWeight,
                                const float solventFraction);
void generateProteinSolventMaskUnderControl(const int iteration,
                                const clipper::Xmap<float>& density,
                                clipper::Xmap<int>& mask,
                                const float sigmaWeight,
                                const float solventFraction,
                                const string updateSpeed);
void generateRefinedProteinSolventMask(const clipper::Xmap<float>& density,
                                       const clipper::Xmap<int>& mask,
                                       clipper::Xmap<int>& refinedMask,
                                       const float sigmaWeight,
                                       const float shrinkFraction);
// ===== HIO and solvent flattening =====
void applyHybridInputOutput(const int iteration, clipper::Xmap<float>& density);

void applyTransitionHybridInputOutput(const int iteration, clipper::Xmap<float>& density);

void applyContinuousHybridInputOutput(const int iteration, clipper::Xmap<float>& density);

void applyHybridProjectionReflection(const int iteration, clipper::Xmap<float>& density);

void applyModifiedContinuousHybridInputOutput(const int iteration, clipper::Xmap<float>& density);
                           
void applySolventFlattening(clipper::Xmap<float>& density, const clipper::Xmap<int>& mask);

void applySolventFlipping(clipper::Xmap<float>& density, const clipper::Xmap<int>& mask);

// ===== Histogram operations =====
float excludeOutlierDensities(int targetGridPoints,
                            clipper::Xmap<float>& density,
                            clipper::Xmap<int>& mask);
                            
void adjustOutlierDensities(float minDensity, 
                          float maxDensity,
                          clipper::Xmap<float>& density);
                          
void combineSmallBinsIntoLargeBins(const std::string& region,
                                   const clipper::Xmap<float>& density,
                                   const clipper::Xmap<int>& mask);

void catchLargeBinBoundary(const std::string& region,
                           const clipper::Xmap<float>& density,
                           const clipper::Xmap<int>& mask);

void catchLargeBinBoundaryExcludeHeavyAtom(const std::string& region,
                           const clipper::Xmap<float>& density,
                           const clipper::Xmap<int>& mask);
                                  
void computeHistogramBoundaries();

// ===== Histogram matching =====
void calculateDensityStatistics(const std::string& region, 
                               const clipper::Xmap<float>& density,
                               const clipper::Xmap<int>& mask,
                               const bool printMessage=false );
                               
void assignDensityValuesToBins(const std::string& region, 
                             const clipper::Xmap<float>& density,
                             const clipper::Xmap<int>& mask);
                                  
void applyHistogramMatching(clipper::Xmap<float>& density,
                          const clipper::Xmap<int>& mask);

// ===== Difference map algorithm =====
void applyDifferenceMap(const int iteration, clipper::Xmap<float>& density);

void applyModifiedDifferenceMap(const int iteration, clipper::Xmap<float>& density);

void applyAveragedSuccessiveReflections(const int iteration, clipper::Xmap<float>& density);

void applyRelaxedAveragedAlternatingReflections(const int iteration, clipper::Xmap<float>& density);

void applyModifiedRelaxedAveragedAlternatingReflections(const int iteration, clipper::Xmap<float>& density);

void applyHybridDifferenceMap(const int iteration, clipper::Xmap<float>& density);

// ===== Convergence monitoring =====
void checkConvergenceWithRworkDetector(const int iteration);
void updateConvergenceFlagAndEliteCount(const int iteration);

void copyPhaseBToPhaseA(clipper::HKL_data<clipper::data64::F_phi>& fpA,
                        const clipper::HKL_data<clipper::data64::F_phi>& fpB);

// ===== Genetic algorithm data structures and functions =====

// Structure to hold density and R-factor information for genetic algorithm
struct GeneticPopulationMember {
    int mpiRank;                              // MPI process rank
    int originIndex;                          // Origin translation index
    bool isInverted;                          // Whether structure is inverted
    bool hasConverged;                        // Convergence flag
    float rWork;                              // R-work value
    float selectionProbability;               // Natural selection probability
    std::unique_ptr<float[]> densityData;     // Current density data
    std::unique_ptr<float[]> densityBuffer;   // Buffer for crossover/mutation
    
    // Constructor
    GeneticPopulationMember()
        : mpiRank(0), originIndex(0), isInverted(false),
          hasConverged(false), rWork(1.0f), selectionProbability(0.0f) {}
    
    // Allocate memory for density arrays
    void allocate(int genomeLength) {
        densityData.reset(new float[genomeLength]());
        densityBuffer.reset(new float[genomeLength]());
    }
};

// Genetic algorithm parameters structure
struct GeneticAlgorithmParams {
    int populationSize;           // Total population size (MPI size)
    int genomeLength;             // Length of genome (grid point count)
    int eliteCount;               // Number of converged/elite members
    int crossoverPointCount;      // Number of crossover points
    float crossoverGeneRatio;     // Ratio of genes to crossover
    float mutationGeneRatio;      // Ratio of genes to mutate
};

// ===== Core genetic algorithm functions =====

// Initialize genetic algorithm population
void initializeGeneticPopulation(
    std::vector<GeneticPopulationMember>& population,
    int populationSize,
    int genomeLength);

// Evaluate fitness and sort population by R-work
void evaluateAndSortPopulation(
    std::vector<GeneticPopulationMember>& population);

// Update elite count based on convergence criteria
void updateEliteCount(
    const int currentIteration,
    const std::vector<GeneticPopulationMember>& population);

// Align density origins to match best solution
void alignDensityOrigins(
    std::vector<GeneticPopulationMember>& population);

// Compute selection probabilities based on diversity and fitness
void computeSelectionProbabilities(
    std::vector<GeneticPopulationMember>& population);

// Perform genetic crossover between population members
void performGeneticCrossover(
    std::vector<GeneticPopulationMember>& population,
    const GeneticAlgorithmParams& params);

// Perform genetic mutation on population members
void performGeneticMutation(
    std::vector<GeneticPopulationMember>& population,
    const GeneticAlgorithmParams& params);

// ===== Helper functions =====

// Copy density map to array format
void copyDensityToArray(
    const clipper::Xmap<float>& densityMap,
    float* densityArray,
    int arrayLength);

// Copy array to density map format
void copyArrayToDensity(
    const float* densityArray,
    clipper::Xmap<float>& densityMap,
    int arrayLength);

// Calculate Euclidean distance between two density distributions
float calculateDensityDistance(
    const float* densityA,
    const float* densityB,
    int length);

// Apply origin translation to density
void translateDensityOrigin(
    clipper::Xmap<float>& densityMap,
    const int& originIndex,
    const bool& isInverted);

// Comparison function for sorting population by R-work
bool compareByRwork(
    const GeneticPopulationMember& memberA,
    const GeneticPopulationMember& memberB);

// ===== Main genetic algorithm execution =====
// Execute one generation of genetic algorithm
void executeGeneticAlgorithmGeneration(
    const int currentIteration,
    const GeneticAlgorithmParams& params,
    clipper::Xmap<float>& density);

//===== Calculate Phase Error for Elite/Average Density =====
void calculatePhaseErrorForEliteDensity(
    std::vector<GeneticPopulationMember>& population,
    const int currentIteration);

#endif // FUNCTIONS_H
