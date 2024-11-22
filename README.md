# NASA TurboFan Predictive Maintenance Analysis

## Project Overview
This repository contains MATLAB code and analysis for NASA's Turbofan Engine Degradation Simulation Dataset. The project focuses on predictive maintenance by optimizing sensor data acquisition and employing advanced modeling techniques to predict the Remaining Useful Life (RUL) of engines.

## Dataset
The project utilizes NASA's Turbofan Engine Degradation Simulation Dataset, which includes:
- Four different scenarios (`trainFD001.mat` to `trainFD004.mat`), each with unique operating conditions and fault modes.
- Sensor data for 21 variables across multiple engines with varying lengths of time series data:
  - Scenario 1: 100 engines
  - Scenario 2: 260 engines
  - Scenario 3: 100 engines
  - Scenario 4: 249 engines

## Key Challenges
- Non-uniform time series lengths between engine runs.
- High dimensionality with multicollinearity between sensor variables.
- Variables with unknown physical meanings, necessitating technical rather than domain-specific analysis.

## Project Structure
The repository is organized as follows:

- **Data Processing Scripts**:
  - `Data_Pretreatment.m`: Preprocesses and cleans raw data, calculates RUL, and applies centering and scaling techniques.
  - `Data_Visualization.m`: Visualizes sensor data, trends, and preprocessed results.
- **Modeling Scripts**:
  - `PLS_modelling_linear.mlx`: Implements Partial Least Squares (PLS) regression for RUL prediction.
- **Dataset Files**:
  - `trainFD001.mat` to `trainFD004.mat`: Training datasets for each scenario.

## Features
- Preprocessing to handle sensor variability and constant values.
- Data visualization for exploratory analysis.
- Partial Least Squares regression for predictive modeling.
- Scenario-specific sensor analysis to optimize model simplicity and accuracy.

## Methodology
1. **Data Preprocessing**:
   - Calculated RUL as `RUL = max(cycles) - current_cycle`.
   - Excluded constant sensors to reduce redundancy and improve efficiency.
   - Standardized sensor data using z-scores for fair contribution in modeling.
2. **Modeling**:
   - Applied PLS regression to identify latent variables and key sensors.
   - Cross-validation to optimize model parameters like latent variables and training cycles.
3. **Validation**:
   - Reserved 20% of data as a holdout set for testing.
   - Evaluated model performance using metrics like RMSE, R², and explained variance.

## Key Results
- **Optimal Training Cycles**: Varied by scenario (e.g., 55 cycles for Scenario 1, 115 for Scenario 4).
- **Latent Variables**: Ranged from 2 to 9 across scenarios, capturing >99% of RUL variance.
- **Critical Sensors**:
  - Sensors 4 and 11 were consistently influential across scenarios.
  - Scenario-specific sensors (e.g., Sensors 8 and 18 in Scenario 4) provided additional insights.

## Requirements
- MATLAB (recommended: R2019b or later).
- MATLAB Statistics and Machine Learning Toolbox.
- MATLAB Signal Processing Toolbox.

## Usage
1. **Data Preprocessing**:
   - Run `Data_Pretreatment.m` to prepare and standardize the dataset.
2. **Data Exploration**:
   - Use `Data_Visualization.m` to examine sensor behaviors and trends.
3. **Model Training**:
   - Execute `PLS_modelling_linear.mlx` to build and validate PLS models.

## Contributors
- Jere Arokivi [@Fl4yd](https://github.com/Fl4yd)
- Erik Brückner [@ikceurb](https://github.com/ikceurb)
- Nanami Kuramochi

## References
- NASA Prognostics Center of Excellence (PCoE) Dataset
- Relevant academic papers and project resources

## Future Work
- Explore nonlinear modeling techniques (e.g., Kernelized PLS) to address identified nonlinear patterns.
- Incorporate domain knowledge to refine variable selection and enhance model interpretability.
