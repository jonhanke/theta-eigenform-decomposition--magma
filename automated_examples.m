// ================================================================
// Some Automated Examples for theta-eigenform-decomposition--magma 
// (c) 2011 Jonathan Hanke
// Released under GNU Public License (GPL) version 2.0
// ================================================================





//========================================================================================================
//=======================================  Automated Examples  ===========================================
//========================================================================================================


//time compute_cusp_constant(KneserForm);
//time compute_cusp_constant(Form6414);




 /*
// Compute some Auxiliary cuspidal constants
time make_cusp_bounds_for_auxiliary_forms([27..40]);
quit;
 */


/*
// Compute some Basic  cuspidal constants
time make_cusp_bounds_for_basic_forms([2786..2900]);
quit;
*/



 /*
// Automated computation for the big form! =)
SetLogFile("10-27-2005____Magma__jon_cusp7.m____Form6414_cusp_constant_computation__take3.txt");
time compute_cusp_constant(Form6414);
UnsetLogFile();
quit;
 */

/*
// Automated computation for the Kneser form! =)
SetLogFile("9-4-2005____Magma__jon_cusp4.m____KneserForm_cusp_constant_computation.txt");
time compute_cusp_constant(KneserForm);
UnsetLogFile();
quit;
*/


/*
// Automated computation for all Auxiliaries
Auxiliary_cusp_bounds := make_cusp_bounds_for_list_of_forms(Auxiliary_Quaternaries, [1..#Auxiliary_Quaternaries]);
SetLogFile("Auxiliary_Logs/ALL_AUXILIARIES.txt");
print Auxiliary_cusp_bounds;
UnsetLogFile();
quit;
*/
