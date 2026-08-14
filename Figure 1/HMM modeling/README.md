# HMM Modeling

Pipeline for detecting play states in rat social interactions using Hidden Markov Models. The workflow prepares behavioral features from tracking and vocalization data, fits HMMs in Python, and analyzes the resulting state sequences against manual play-bout annotations.

## Pipeline Order

### 1. `CreateHMMFIiles.m` (MATLAB script)
Batch wrapper that drives the full pipeline across all animal sessions. Loops over all animal folders in the raw data directory and:
- Calls `Create_HMM_inputs` to generate the `.npy` feature matrices for each session.
- Calls `Load_HMM_outputs_and_analyze` to process HMM results after the Python fitting step has been run.

### 2. `Create_HMM_inputs.m` (MATLAB function)
Prepares the input feature matrix for HMM fitting, for a given animal session.
- Loads synchronization, tracking, behavior annotation, and ultrasonic call detection data.
- Synchronizes video timestamps to audio timestamps using a linear regression model.
- Restricts analysis to "Partners session" time windows and loops over each partner.
- Computes 14 behavioral variables at 10 ms resolution from tracking data:
  - **Focal animal:** speed, angular speed, angular acceleration, acceleration.
  - **Partner:** speed, angular speed, angular acceleration, acceleration.
  - **Relative (between-animal):** distance, speed, angular speed, angular acceleration, acceleration.
  - **Vocalization:** call presence (binary).
- Smooths all variables with a Gaussian convolution kernel (0.5 s half-width) and z-scores them.
- Labels each time bin as play or non-play using manually annotated play bouts.
- Saves the z-scored feature matrix (without the play label) as a `.npy` file for Python, and a `.mat` file with timestamps and variable names.

### 3. `python hmm estimate.txt` (Python script)
Fits Hidden Markov Models to the feature matrices produced by `Create_HMM_inputs`.
- Uses the `ssm` library to fit sticky HMMs with diagonal Gaussian emissions.
- Fits models with K = 2, 3, 4, and 5 states (14-dimensional observations, 500 EM iterations).
- Saves the most likely state sequence and transition matrix for each K.
- Computes log-likelihood, AIC, and BIC for model comparison.

### 4. `Load_HMM_outputs_and_analyze.m` (MATLAB function)
Loads the HMM state sequences from Python and performs all downstream analyses. Loops over each partner session for a given animal:
- **Confusion matrix (section 6):** Compares HMM 2-state output against manual play-bout labels (TP, TN, FP, FN), auto-swapping state labels if needed.
- **HMM onset/offset extraction (section 7):** Detects state transition times for both 2-state and 3-state models.
- **Behavior classification (section 8):** For each annotated behavior type, determines whether it falls within, at the transition boundary, or outside HMM play states; saves proportions.
- **Behavior peri-event analysis (section 9):** Builds PSTH-style matrices of behavior occurrence aligned to HMM state onsets and offsets (2- and 3-state models).
- **Variable peri-event analysis (section 10):** Builds peri-event matrices of each continuous behavioral variable aligned to HMM state onsets/offsets; computes mean responses and response indices.
- **State re-assignment (sections 11-11.B):** Visualizes the 3-state assignment ordered by HMM state duration; auto-reassigns state labels so the state with the most play-behavior overlap gets the highest label, then allows manual correction.
- **GLM play prediction (sections 12-14):** Fits a binomial GLM (logistic regression) predicting play-bout presence from the 14 behavioral variables; estimates play probability for each HMM state.
- **Call analysis (section 16):** Builds peri-event matrices of ultrasonic call occurrence and acoustic properties (frequency, slope, sinuosity, etc.) aligned to HMM state onsets/offsets.

## Saved Outputs (per animal and partner)
| File | Contents |
|------|----------|
| `confusion_matrix` | HMM vs. manual play-bout classification performance |
| `play_behavior_struct` | Proportion of each behavior type within/outside/at transitions of HMM states |
| `behavior_onset_offset_struct` | Peri-event behavior matrices aligned to HMM state boundaries |
| `variable_onset_struct` | Peri-event behavioral-variable matrices and mean responses |
| `prediction_struct` | GLM model, predicted play probabilities, state assignments |
| `call_struct` | Peri-event call matrices and acoustic properties |

## Dependencies
- MATLAB (with Statistics and Machine Learning Toolbox for `fitglm`, `zscore`)
- `readNPY` / `writeNPY` (npy-matlab)
- `play_bout` function
- Python with `numpy` and `ssm`
