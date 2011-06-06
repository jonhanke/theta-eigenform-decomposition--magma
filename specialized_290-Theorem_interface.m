// ======================================================
// Specialized 290-Theorem individual and batch commands 
// for running theta-eigenform-decomposition--magma 
// (c) 2011 Jonathan Hanke
// Released under GNU Public License (GPL) version 2.0
// ======================================================








// ============================================================================================================
// ======================== These routines are called by runaux.py and runbasic.py ============================
// ============================================================================================================



// This makes the cuspidal bounds for the Auxiliary Escalator Forms
function make_cusp_bound_for_auxiliary_form(i, logfile, mndeg)

  // Some Basic Defintions
  Form_List := Auxiliary_Quaternaries;    // THIS IS A GLOBAL VARIABLE
  Q4 := RMatrixSpace(Rationals(),4,4);
  Z4 := RMatrixSpace(IntegerRing(),4,4);

  // Start logging output
  SetLogFile(logfile);

  // Print MAGMA Info
  print "MAGMA Info:";
  print "-----------";
  print "MAGMA Version:", GetVersion();
  print "MAGMA Current Memory Usage:", GetMemoryUsage();
  print "Current MAGMA process id:", Getpid();
  print "==============================================";

  // Find the cuspidal decomposition and bound
  print "";
  print " Starting to compute the cuspidal information for Auxiliary Form #", i;
  print "";
  print " Read in the vector:", Form_List[i];
  print "";
  QQ := Z4 ! (2 *(Q4 ! Form_List[i]));
  time cusp_bound := compute_cusp_constant(QQ, mndeg);
  print "";
  print "==============================================";

  // Print MAGMA Info again
  print "MAGMA Info:";
  print "-----------";
  print "MAGMA Version:", GetVersion();
  print "MAGMA Current Memory Usage:", GetMemoryUsage();
  print "==============================================";

  // Stop logging output
  UnsetLogFile();

  // Return the bound
  return cusp_bound;

end function;



// This makes the cuspidal bounds for the Basic Escalator Forms
function make_cusp_bound_for_basic_form(i, logfile, mndeg)

  // Some Basic Defintions
  Form_List := Basic_Quaternaries;    // THIS IS A GLOBAL VARIABLE
  Q4 := RMatrixSpace(Rationals(),4,4);
  Z4 := RMatrixSpace(IntegerRing(),4,4);

  // Start logging output
  SetLogFile(logfile);

  // Print MAGMA Info
  print "MAGMA Info:";
  print "-----------";
  print "MAGMA Version:", GetVersion();
  print "MAGMA Current Memory Usage:", GetMemoryUsage();
  print "Current MAGMA process id:", Getpid();
  print "==============================================";

  // Find the cuspidal decomposition and bound
  print "";
  print " Starting to compute the cuspidal information for Basic Form #", i;
  print "";
  print " Read in the vector:", Form_List[i];
  print "";
  QQ := Z4 ! (2 *(Q4 ! Form_List[i]));
  time cusp_bound := compute_cusp_constant(QQ, mndeg);
  print "";
  print "==============================================";

  // Print MAGMA Info again
  print "MAGMA Info:";
  print "-----------";
  print "MAGMA Version:", GetVersion();
  print "MAGMA Current Memory Usage:", GetMemoryUsage();
  print "";

  // Stop logging output
  UnsetLogFile();

  // Return the bound
  return cusp_bound;

end function;



// ============================================================================================================
// ============================================================================================================



// This makes the cuspidal bounds for the Auxiliary Forms
function make_cusp_bounds_for_auxiliary_forms(index_list:  OUT_STRING := "/tmp/Auxiliary_Logs/Auxilliary_")
     
  // Where to put the output?
  //OUT_STRING := "/tmp/Auxiliary_Logs/Auxilliary_";   // Moved to be an optional parameter
  //OUT_STRING := "Basic_Logs/Basic_";   


  // Some Basic Defintions
  Form_List := Auxiliary_Quaternaries;    // THIS IS A GLOBAL VARIABLE
  Q4 := RMatrixSpace(Rationals(),4,4);
  Z4 := RMatrixSpace(IntegerRing(),4,4);

  Cusp_list := [-1.0 : i in [1..#Form_List]];

  // Find all cusp decompositions and bounds
  for i in index_list do
     SetLogFile(OUT_STRING cat IntegerToString(i) cat "__cusp_info.txt");
     print " Starting to compute the cuspidal information for Auxiliary Form #", i;
     print "";
     print " Read in the vector:", Form_List[i];
     print "";
     QQ := Z4 ! (2 *(Q4 ! Form_List[i]));

     time cusp_bound := compute_cusp_constant(QQ);

     print "Adding the cusp bound", cusp_bound, "to the big list";
     Cusp_list[i] := cusp_bound;
     UnsetLogFile();
  end for;

  // Return them all
  return Cusp_list;

end function;


// This makes the cuspidal bounds for the Basic Forms
function make_cusp_bounds_for_basic_forms(index_list:  OUT_STRING := "/tmp/Basic_Logs/Basic_")
     
  // Where to put the output?
  //OUT_STRING := "/tmp/Basic_Logs/Basic_";   // Moved to be an optional parameter
  //OUT_STRING := "Basic_Logs/Basic_";   


  // Some Basic Defintions
  Form_List := Basic_Quaternaries;    // THIS IS A GLOBAL VARIABLE
  Q4 := RMatrixSpace(Rationals(),4,4);
  Z4 := RMatrixSpace(IntegerRing(),4,4);

  Cusp_list := [-1.0 : i in [1..#Form_List]];

  // Find all cusp decompositions and bounds
  for i in index_list do
     SetLogFile(OUT_STRING cat IntegerToString(i) cat "__cusp_info.txt");
     print " Starting to compute the cuspidal information for Basic Form #", i;
     print "";
     print " Read in the vector:", Form_List[i];
     print "";
     QQ := Z4 ! (2 *(Q4 ! Form_List[i]));

     time cusp_bound := compute_cusp_constant(QQ);

     print "Adding the cusp bound", cusp_bound, "to the big list";
     Cusp_list[i] := cusp_bound;
     UnsetLogFile();
  end for;

  // Return them all
  return Cusp_list;

end function;
