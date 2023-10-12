# Adaptive-GRT
Repository of code, data, and analysis scripts for running Adaptive General Recognition Theory experiments (Glavan et al.)

## Import the AGRT library to use the following classes and functions:
- class **AGRTHandler**: PsychoPy handler object that initializes and manages the internal adaptive algorithm. Similar to PsychoPy's StairHandler.
- class **GRTHandler**: PsychoPy handler object that implements the complete identification paradigm trial structure. Similar to PsychoPy's TrialHandler
- function **RunAdaptiveBlock**: Creates an instance of AGRTHandler with the given parameters and runs the requested number of trials. Returns stimuli values corresponding to the requested level of accuracy.
- function **RunGRTBlock**: Creates an instance of GRTHandler with the given parameters and runs the requested number of trials. Returns None.
- function **RunAdaptiveGRTExperiment**: Uses the above functions to first run an adaptive block and then a main GRT block that uses the adapted stimuli.
- function **GRTSubjectModel**: Example trial function used in the simulation study from Glavan et al.

## Additional Files in this repository:
- **AGRT_Exp_x_x_x.py**: Runs the human subjects experiments from Glavan et al.
- **AGRT Human Analysis.R**: Analysis of the human subjects experiments from Glavan et al.
- **AGRT Quick Accuracy Check.R**: Script that prints the accuracy of human subjects in the pilot condition. Used to inform the selection of stimuli.
- **test_simulation.py**: Runs the simulation study from Glavan et al.
- **AGRT Python Simulation Study.R**: Analysis of the simulation study from Glavan et al.
