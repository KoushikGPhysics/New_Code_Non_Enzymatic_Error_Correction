Code for “Non-enzymatic error correction in self-replicators without extraneous energy supply”

This repository contains all MATLAB scripts used to generate the results and figures for the manuscript:

/src
    main_simulation.m          – Runs the full template-directed polymerization model
    kinetics_functions.m       – Contains rate equations and proofreading steps
    energy_model.m             – Computes ΔG terms and kinetic discrimination
    error_rate_calculation.m   – Computes copying error η under various parameters
    plotting_scripts.m         – Generates all manuscript figures

/data
    input_parameters.mat       – Model parameters used in the manuscript
    example_output.mat         – Sample output structure for testing

