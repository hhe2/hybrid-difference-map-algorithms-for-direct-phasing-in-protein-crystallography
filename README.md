# hybrid-difference-map-algorithms-for-direct-phasing-in-protein-crystallography
Hybrid Difference Map Algorithms for Direct Phasing in Protein Crystallography

# Crystallographic Phase Reconstruction GUI

A modern web-based interface for iterative crystallographic phase retrieval from X-ray diffraction amplitude data using hybrid input-output (HIO), difference map (DM) and hybrid difference map (HDM) algorithms.


## Features

✨ **Modern Web Interface**
- Intuitive file upload with drag-and-drop support
- Real-time parameter configuration
- Live job monitoring with progress tracking
- Interactive results visualization

🚀 **Powerful Backend**
- MPI-based parallel processing
- RESTful API for job management
- Real-time output parsing
- Multiple job tracking and management

🔬 **Advanced Algorithms**
- Hybrid Input-Output (HIO) method
- Difference Map algorithm
- Histogram matching
- Solvent flattening
- Genetic algorithm optimization
- Molecular replacement initialization

📊 **Quality Metrics**
- R-work and R-free monitoring
- Phase error tracking
- Mask correlation coefficients
- Real-time convergence visualization

## Quick Start

### Automated Setup

```bash
# Clone or download the project
cd crystallography-phase-reconstruction

# Run setup script
chmod +x setup.sh
./setup.sh

# Start the application
./start_all.sh
```

The application will open in your browser at `http://localhost:3000`

### Manual Setup

#### 1. Compile C++ Program

```bash
make clean
make
```

#### 2. Set Up Python Backend

```bash
python3 -m venv venv
source venv/bin/activate
pip install flask flask-cors werkzeug
```

#### 3. Set Up React Frontend

```bash
npx create-react-app frontend
cd frontend
npm install lucide-react
```

#### 4. Start Services

**Terminal 1 - Backend:**
```bash
cd backend
python server.py
```

**Terminal 2 - Frontend:**
```bash
cd frontend
npm start
```

## Usage

### 1. Prepare Input Files

Place your crystallographic data in the `PDB_file/` directory:

```
PDB_file/
├── 2uxj-sf.mtz           # Diffraction amplitudes (required)
├── 2uxj_fmodel.mtz       # Reference phases (optional)
├── 2uxj.pdb              # Structure file (optional)
└── reference_fmodel.mtz  # Histogram reference (optional)
```

### 2. Upload Files

1. Open the web interface at `http://localhost:3000`
2. Navigate to "File Upload" tab
3. Upload your MTZ and PDB files
4. Wait for upload confirmation

### 3. Configure Parameters

Switch to "Parameters" tab and adjust:

**Structure Information:**
- PDB Code: `2uxj`
- Histogram Reference: `2uxj`

**Iteration Settings:**
- Number of iterations: `200` (typical: 100-500)
- Grid sampling: `1.0` (higher = more detail, slower)
- Resolution high: `2.25` Å
- Resolution low: `15.0` Å

**Density Modification:**
- Solvent content: `0.68` (typical: 0.4-0.7)
- HIO beta: `0.7` (typical: 0.6-0.9)
- Sigma weights: `4.0` → `3.0`

**Algorithm Options:**
- ☑ Initialize from PDB structure
- ☐ Enable genetic algorithm
- ☐ Initialize from random map
- ☐ Initialize from molecular replacement

**MPI Configuration:**
- Number of processes: `4` (set to available CPU cores)

### 4. Run Job

1. Click "Run Phase Reconstruction"
2. Switch to "Job Monitor" tab
3. Watch real-time progress:
   - Iteration counter
   - R-work and R-free values
   - Phase errors
   - Console output

### 5. Download Results

When complete, download output files:

| File | Description |
|------|-------------|
| `R_work.txt` | R-work values per iteration |
| `R_free.txt` | R-free values for validation |
| `mean_phase_error.txt` | Phase errors over iterations |
| `mask_correlation.txt` | Protein mask correlation |
| `structure_factors.mtz` | Final structure factors with phases |
| `phase.cif` | CIF file for visualization in Coot/PyMOL |

## Architecture

```
┌─────────────────────────────────────────────────────────┐
│                   React Frontend                         │
│  - File upload interface                                 │
│  - Parameter configuration                               │
│  - Real-time job monitoring                              │
└──────────────────┬──────────────────────────────────────┘
                   │ HTTP/REST API
┌──────────────────▼──────────────────────────────────────┐
│                   Flask Backend                          │
│  - File management                                       │
│  - Job queue and tracking                                │
│  - MPI process execution                                 │
│  - Output parsing                                        │
└──────────────────┬──────────────────────────────────────┘
                   │ subprocess/MPI
┌──────────────────▼──────────────────────────────────────┐
│              C++ Phase Reconstruction                    │
│  - MPI parallel processing                               │
│  - HIO/Difference Map algorithms                         │
│  - FFT transformations (FFTW)                            │
│  - Density modification (Clipper)                        │
└─────────────────────────────────────────────────────────┘
```

## API Reference

### Endpoints

#### File Management

**Upload File**
```http
POST /api/upload
Content-Type: multipart/form-data

Body:
  file: <binary>
  type: "diffractionMtz" | "histogramMtz" | "referenceMtz" | "pdbFile"

Response:
{
  "success": true,
  "filename": "2uxj-sf.mtz",
  "filepath": "/path/to/file"
}
```

**Download File**
```http
GET /api/download/<filename>

Response: <binary file>
```

#### Job Management

**Submit Job**
```http
POST /api/submit
Content-Type: application/json

Body:
{
  "params": {
    "pdbCode": "2uxj",
    "numIter": 200,
    "solvCont": 0.68,
    ...
  },
  "numProcesses": 4
}

Response:
{
  "success": true,
  "job_id": "job_1_1634567890",
  "message": "Job submitted successfully"
}
```

**Get Job Status**
```http
GET /api/status/<job_id>

Response:
{
  "job_id": "job_1_1634567890",
  "status": "running",
  "progress": 45,
  "current_iteration": 90,
  "total_iterations": 200,
  "r_work": [0.45, 0.42, 0.39, ...],
  "r_free": [0.47, 0.44, 0.41, ...],
  "phase_errors": [65.3, 58.2, 52.1, ...],
  "recent_output": ["Iteration: 90 | R-work: 0.3912...", ...]
}
```

**Cancel Job**
```http
POST /api/cancel/<job_id>

Response:
{
  "success": true,
  "message": "Job cancelled"
}
```

**List All Jobs**
```http
GET /api/jobs

Response:
{
  "jobs": [
    {
      "job_id": "job_1_1634567890",
      "status": "completed",
      ...
    }
  ]
}
```

#### System

**Health Check**
```http
GET /api/health

Response:
{
  "status": "healthy",
  "timestamp": "2025-10-05T12:34:56",
  "active_jobs": 2
}
```

**Validate Setup**
```http
POST /api/validate

Response:
{
  "valid": true,
  "issues": [],
  "mpi_available": true,
  "base_dir": "/path/to/project"
}
```

## Configuration

### params.txt Format

```ini
# Structure Information
PDB_CODE 2uxj
PDB_HIST 2uxj

# Iteration Parameters
numIter 200
Grid_sampling 1.0
resoCutoff 2.25
resoCutoff_low 15.0

# Density Modification
solvCont 0.68
hioBeta 0.7
sigmWeigAvgIni 4.0
sigmWeigAvgFina 3.0

# Algorithm Flags (0=false, 1=true)
genetic_algorithm 0
start_from_pdb 1
start_from_randmap 0
start_from_mr 0
```

### Environment Variables

```bash
# Optional: Override default paths
export CRYSTALLOGRAPHY_UPLOAD_DIR=/custom/path/PDB_file
export CRYSTALLOGRAPHY_OUTPUT_DIR=/custom/path/output
export CRYSTALLOGRAPHY_PROGRAM=/custom/path/program
```

## Algorithm Details

### Hybrid Input-Output (HIO)

The HIO algorithm iteratively refines phases by:
1. Applying known constraints in real space (solvent regions)
2. Replacing calculated amplitudes with observed amplitudes
3. Using feedback parameter β to modify density outside protein mask

```
ρ_new(x) = {
  ρ_calc(x)           if x in protein region
  ρ_old(x) - β·ρ_calc(x)   if x in solvent region
}
```

### Difference Map

Alternates between two constraint sets:
- Fourier space: observed amplitudes
- Real space: density positivity and solvent flatness

### Histogram Matching

Adjusts density histogram to match reference distribution from known structures at similar resolution.

### Solvent Flattening

Forces density in solvent regions to a constant low value, typically applied in final iterations.

## Performance Optimization

### CPU Utilization

```bash
# Use all available cores
numProcesses = $(nproc)

# Or specific number
numProcesses = 8
```

### Memory Management

- **Grid sampling 1.0**: ~2-4 GB per process
- **Higher resolution**: More memory required
- **Genetic algorithm**: Additional 20-30% memory

### Iteration Speed

Typical performance on modern hardware:
- 2.5 Å resolution: ~1-2 seconds per iteration
- 2.0 Å resolution: ~3-5 seconds per iteration
- 1.5 Å resolution: ~10-15 seconds per iteration

## Troubleshooting

### Common Issues

**1. Program won't start**
```bash
# Check executable permissions
chmod +x program

# Test MPI
mpirun -np 1 ./program

# Check libraries
ldd program
```

**2. Upload fails**
```bash
# Check directory permissions
chmod 755 PDB_file/
chmod 755 output/

# Check disk space
df -h
```

**3. Job stuck at 0%**
```bash
# Check MPI is working
mpirun -np 2 hostname

# Check logs
tail -f logs/backend.log
```

**4. High memory usage**
```bash
# Reduce grid sampling
Grid_sampling 0.8

# Use fewer MPI processes
numProcesses 2

# Reduce resolution
resoCutoff 2.5
```

**5. CORS errors in browser**
```bash
# Ensure Flask-CORS is installed
pip install flask-cors

# Check proxy in frontend/package.json
"proxy": "http://localhost:5000"
```

### Validation Checklist

- [ ] MPI installed and working (`mpirun --version`)
- [ ] Program compiled (`ls -l program`)
- [ ] Python backend running (`curl http://localhost:5000/api/health`)
- [ ] Frontend accessible (`curl http://localhost:3000`)
- [ ] Input files in PDB_file/ directory
- [ ] Output directory writable

## Development

### Project Structure

```
.
├── backend/
│   └── server.py              # Flask API server
├── frontend/
│   ├── public/
│   └── src/
│       ├── App.jsx            # Main React component
│       └── index.js
├── src/                       # C++ source
│   ├── main.cpp
│   ├── functions.cpp
│   ├── utils.cpp
│   ├── globalvars.cpp
│   └── dataimporter.cpp
├── include/                   # C++ headers
├── PDB_file/                  # Input data
├── output/                    # Results
├── logs/                      # Log files
├── Makefile                   # C++ compilation
├── params.txt                 # Default parameters
├── setup.sh                   # Automated setup
├── start_all.sh               # Start both services
└── README.md                  # This file
```

### Building from Source

```bash
# Clean build
make clean
make

# Debug build
make CXXFLAGS="-g -O0"

# Verbose output
make V=1
```

### Running Tests

```bash
# Test MPI execution
mpirun -np 2 ./program

# Test backend API
cd backend
python -m pytest tests/

# Test frontend
cd frontend
npm test
```

### Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Dependencies

### C++ Libraries
- **CCP4/Clipper**: Crystallographic computations
- **FFTW**: Fast Fourier transforms
- **MPI**: Parallel processing

### Python Packages
- **Flask**: Web framework
- **Flask-CORS**: Cross-origin resource sharing
- **Werkzeug**: WSGI utilities

### JavaScript Packages
- **React**: Frontend framework
- **lucide-react**: Icon library

## Citations

If you use this software in your research, please cite:

```bibtex
@software{crystallography_gui,
  title={Crystallographic Phase Reconstruction GUI},
  author={Your Name},
  year={2025},
  url={https://github.com/yourusername/crystallography-gui}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- CCP4 for crystallographic libraries
- FFTW developers for FFT implementation
- Open source community for frameworks and tools

## Support

- **Documentation**: See `docs/` directory
- **Issues**: GitHub Issues
- **Email**: support@example.com
- **Forum**: https://forum.example.com

## Roadmap

- [ ] Add support for SAD/MAD phasing
- [ ] Implement automated model building
- [ ] Add 3D visualization in browser
- [ ] Support for multiple crystal forms
- [ ] Batch processing interface
- [ ] Cloud deployment option
- [ ] Docker containerization

---

**Version**: 1.0.0  
**Last Updated**: October 2025  
**Maintainer**: Your Name <your.email@example.com>
