// ============================================================
// Initialization file for theta-eigenform-decomposition--magma 
// (c) 2011 Jonathan Hanke
// Released under GNU Public License (GPL) version 2.0
// ============================================================

// Temporarily set the Theta_Decomposition Directory to the current path
THETA_DIR_ABS_PATH := GetEnvironmentValue("THETA_DECOMP_MAGMA_PATH");      // Set this environment variable!
OLD_PATH := GetPath();
SetPath(THETA_DIR_ABS_PATH);



// Load the example quadratic forms and related information
load "example_forms.m";


// Load simpler quadratic form utilities 
//   (Level, Theta series, and Eisenstein Series by averaging)
load "qf_utils.m";


// Load the main normalized cusp constant/decomposition routine
//   (compute_cusp_const)
load "decomposition_cusp_constant.m";


// Load the cusp dimension routine 
//   (computes the dimension of the full space of cusp forms to measure its complexity)
load "cusp_dim.m";


// Load the specialized 290-Theorem interface routines 
load "specialized_290-Theorem_interface.m";


// Load the automated examples (all commented out by default)
load "automated_examples.m";




// Restore the old path
SetPath(OLD_PATH);
