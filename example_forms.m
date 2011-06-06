
// ======================================================
// Examples for theta-eigenform-decomposition--magma 
// (c) 2011 Jonathan Hanke
// Released under GNU Public License (GPL) version 2.0
// ======================================================



// Some convenient MAGMA space defintitions
Q4 := RMatrixSpace(Rationals(),4,4);
Z4 := RMatrixSpace(IntegerRing(),4,4);




// --------------------------------------------------------------------------------
// Define some example Quadratic forms, and load all 290-Theorem Quadratic forms
// --------------------------------------------------------------------------------

// A few smaller interesting quadratic forms used for testing (Note: The Hessian matrix is given -- it's twice the Gram matrix!)
KneserForm := Z4 ! [2,0,0,0,0,6,0,0,0,0,10,0,0,0,0,14];
//KneserHALF := Z4 ! [1,0,0,0,0,3,0,0,0,0,5,0,0,0,0,7];
Form5331 := Z4 ! [ 2, 0, -1, -6, 0, 4, 1, 0, -1, 1, 10, 13, -6, 0, 13, 58 ];  //  920.1
Form6414 := Z4 ! [ 2, 0, 0, 0, 0, 4, 1, -1, 0, 1, 8, 3, 0, -1, 3, 62 ];       // 2331.9

// Load the files with the 290 Project quaternaries (NOTE: These matrices are half-integral!)
load "290_QF_files/e4r.m";
load "290_QF_files/104_auxiliary_quaternaries.m";
load "290_QF_files/290-cusp-all.m";
stein_data := Sort(data);

//Auxilliary_form_class_numbers := [#GenusRepresentatives(LatticeWithGram(2 * Q4 ! Auxiliary_Quaternaries[i])) : i in [1..#Auxiliary_Quaternaries]];
Auxilliary_form_class_numbers := [ 4, 6, 4, 6, 4, 18, 9, 8, 8, 17, 44, 10, 29, 29, 73, 12, 86, 8, 16, 36, 32, 21, 90, 61, 69, 33, 17, 5, 89, 48, 36, 46, 36, 14, 34, 36, 16, 38, 34, 43, 72, 70, 32, 98, 20, 48, 12, 12, 45, 33, 22, 55, 26, 99, 16, 88, 24, 41, 41, 10, 20, 37, 33, 22, 13, 51, 42, 15, 12, 50, 39, 45, 99, 62, 24, 99, 40, 35, 40, 10, 12, 45, 69, 12, 54, 79, 118, 35, 43, 29, 66, 23, 71, 41, 56, 35, 93, 54, 73, 28, 15, 38, 45, 14 ];



// time compute_cusp_constant(2* (Q4!Auxiliary_Quaternaries[6]));

// ========================================================================

/*
 const_range := [1..99];
 exact_constants := [RealField(10)!compute_cusp_constant(2* (Q4!Basic_Quaternaries[i])) : i in const_range];
 approx_constants := [RealField(10)!stein_data[i][4] : i in const_range];
 error_consts := [exact_constants[i] - approx_constants[i] : i in const_range];
*/

// ========================================================================


