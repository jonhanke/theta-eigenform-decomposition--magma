// ======================================================
// theta-eigenform-decomposition--magma  Version 13
// (c) 2011 Jonathan Hanke
// Released under GNU Public License (GPL) version 2.0
// ======================================================











// ---------------------------------------------------------------------------------------------
// Main Routine -- Compute the Cuspidal Decomposition of the cuspidal part of the theta series
// ---------------------------------------------------------------------------------------------



///////////////////////////////////////////////////////////////////////
// Compute the cuspidal constant of any positive definite quaternary
// integer-valued quadratic form, given here by the matrix of 2*Q.
// (Note: This is really the same as the integer-matrix form, 
// but at a lower level, since theta_{2Q}(z) = theta_{Q}(2z).) =)
///////////////////////////////////////////////////////////////////////
function compute_cusp_constant(QQ: mndeg:=1)

// TIMING: Set the current time
t := Cputime();


print "";
print "     ======================";
print "     ===== Basic Info =====";
print "     ======================";
print "";


// DIAGNOSTIC
print "Using the matrix 2*Q =", QQ;
print "";

// Some useful definitions -- used in finding the determinant and level
Q4 := RMatrixSpace(Rationals(),4,4);
Z4 := RMatrixSpace(IntegerRing(),4,4);

// Make the appropriate (quadratic) character and level
N := LevelOfQuadForm(Q4!QQ);  // This is the level
D := Determinant(Z4!QQ);      // This is the determinant (since QQ is even-dim'l)
if IsSquare(D) then      // This deals with Magma's unhappiness about square discriminants... =(
  eps_top := 1;
else 
  eps_top := FundamentalDiscriminant(D);
end if;
eps := KroneckerCharacter(eps_top, RationalField());
print "The form Q has level", N, "and determinant", D/16;


M2 := ModularForms(DirichletGroup(N) ! eps);
S2 := CuspidalSubspace(M2);
Dim_S2 := Dimension(S2);
print "The dimension of the cuspidal subspace is ", Dim_S2;
precision := PrecisionBound(S2);
print "Need the precision: ", precision, "to uniquely identify a cusp form.";
print "";

// TIMING: Check the time taken up to this point, and reset the time.
print "TIMING: Time to setup the spaces of modular forms = ", Cputime(t), " seconds.";
t := Cputime();
print "";
print "";

// Return with a constant of zero if there is no cuspidal subspace! =)
if Dim_S2 eq 0 then
  print "The cuspidal space is 0-dim'l, so are no cusp forms!";
  print "\n This gives the overall cuspidal constant of:", 0;
  print "";
  print "";
  return 0;
end if;



// Make a bound on the degree of of newforms considered...based on the class number:
// ---------------------------------------------------------------------------------
//degree_bound := 2 * #GenusRepresentatives(LatticeWithGram(QQ));   // This was a good idea, but it failed. =(
//degree_bound := 100;   // This now gives a computationally meaningful bound! =)
degree_bound := 1000;   // This now gives a meaningless bound! =)



print "";
print "";
print "     =========================================";
print "     ===== Making newforms of each level =====";
print "     =========================================";
print "";


// Make a list of newforms, and their levels:
// ------------------------------------------
big_oldform_d_list := [d : d in Divisors(N div Conductor(eps))];
level_list := [* *];
form_list := [* *];
for r in [q : q in Divisors(N) | q mod Conductor(eps) eq 0] do
  // Make the new character (which carries the level)
  eps_r := DirichletGroup(r) ! eps;
  //  M_d := ModularForms(r);    // This uses only the trivial character
  M_r := ModularForms(eps_r);    // This uses the correct quadratic character, but it's too big to finish! =(
  //  SetPrecision(S_r, precision+1);   // Set the printing precision
  S_r := CuspidalSubspace(M_r);
  time NewNum_r := NumberOfNewformClasses(S_r);             // THIS TAKES A LONG TIME! =(
  print "There are ", NewNum_r, " classes of (cuspidal) newforms of level", r;
  for i in [1..NewNum_r] do
    time tmp_newform := Newform(S_r, i);
    if Degree(tmp_newform) le degree_bound then   // OLD COMMENT: Here we filter newforms with degree larger than the class number!
      level_list := Append(level_list, r);
      form_list := Append(form_list, tmp_newform);
      print "Class ", i , ": ", CoefficientField(tmp_newform);
    else 
      print "Skipping class ", i , "... since its degree = ", Degree(CoefficientField(tmp_newform)), ">", degree_bound;
      error "ERROR: We are discarding some large newform! =(";      // ERROR MESSAGE IF WE FAIL!!!      
    end if;
  end for;
//print "  form_list = ", form_list;
end for;
print " Using big_oldform_d_list =", big_oldform_d_list;
print " Using level_list =", level_list;

// TIMING: Check the time taken up to this point, and reset the time.
print "TIMING: Time to compute the newforms and their levels = ", Cputime(t), " seconds.";
t := Cputime();



print "";
print "";
print "     =================================================================";
print "     ===== Computing the newforms which contribute to each level =====";
print "     =================================================================";
print "";


// Precompute the allowed_j_lists for each d, and make d's which have some assoicated newform:
// -------------------------------------------------------------------------------------------
allowed_j_list_vec := [ ];
oldform_d_list := [ ];

for d in big_oldform_d_list do

  // Compute the allowed indices j for this new/oldform lift
  print " ";
  print "Starting the computation for the old/newform lift d =", d;
  allowed_j_list := [j : j in [1..#form_list] | N mod (d * level_list[j]) eq 0];
  print "The allowed_j_list is", allowed_j_list;

  // Add it to the vector if it's not empty
  if not IsEmpty(allowed_j_list) then
    oldform_d_list := Append(oldform_d_list, d);
    allowed_j_list_vec := Append(allowed_j_list_vec, allowed_j_list);
  end if;
  // *********  TO DO:  **********
  // If there are no allowed forms, then we should remove d from the list!!!
  // Also, we'll need to modify the level list too...
  // This will save us the trouble of protecting against this all the time!!! =)
  // *****************************

end for;
print " Using oldform_d_list =", oldform_d_list;

// TIMING: Check the time taken up to this point, and reset the time.
print "TIMING: Time to compute the allowed_j_list for all d = ", Cputime(t), " seconds.";
t := Cputime();



print "";
print "     =================================================================";
print "     ===== Computing the first 500 primes not dividing the level =====";
print "     =================================================================";
print "";


// Precompute a list of eligible coefficients to use which allows us to use the fewest number of primes:
// -----------------------------------------------------------------------------------------------------
possible_primes_500 := [];
p := 1;
while #possible_primes_500 lt 500 do                // Precompute 50 possible primes to try. =)
  p := NextPrime(p);
  if N mod p ne 0 then
    possible_primes_500 := Append(possible_primes_500, p);
  end if;
end while;
print "The list of possible primes p prime to N = ", N, " is ", possible_primes_500;



print "";
print "";
print "     =================================================================";
print "     ===== Computing coefficients of each newform for each level =====";
print "     =================================================================";
print "";


// Add a minimal set of coefficients (prime-by-prime) to span the (d=1) eigenspace:
// --------------------------------------------------------------------------------

// Set the number of Fourier coefficients to use
if mndeg lt 400 then
  coeff_bound := 10000;   // Use the bound of 10,000 since the theta function is easy to compute up to there! =)
else
  coeff_bound := 100000;   // Use the bound of 100,000 when we expect it will save us time (It's about a day to compute the theta function!)
end if;

newform__p_coeff_list := [ [* *] : j in [1.. #form_list] ];   // one list for each newform, then the primes for each are stored according to prime_seq. =)
Big_EigenRows_Matrix_list := [* *];
N_min := Conductor(eps);  // Conductor of the level character

good_mm_list_vec := [];
Big_EigenRows_Matrix_list := [* *];

// Run through all possible d's, and construct the minimal set of coefficients.
for ind in [1..#oldform_d_list] do  

  // Some initial conditions 
  r := Cputime();
  d := oldform_d_list[ind];                   // This should be d=1 
  allowed_j_list := allowed_j_list_vec[ind];  // This gives the full list!
  further_lift_num := N div (d * N_min);
  last_rank := 0;
  good_mm_list := [ ];
  old_eligible_coeff_seq := [ ];
  k := 1;  // prime index in possible_primes_50
  rank_is_full_flag := false;


  // Find the minimal set of coefficients
  while rank_is_full_flag eq false do

    p := possible_primes_500[k];


//    print "-------------------------------------";
    print "Computing d =", d, "coefficients using p =", p;
    print "==========================================";

    // Check if that prime hasn't been computed already
    if #newform__p_coeff_list[1] lt k then

      // Compute all p-power coefficients:
      // ---------------------------------
      for j in [1.. #form_list] do

        // Compute the p-th eigenvalue
        print "Computing the coefficient at ", p, "for the newform at j = ", j;
        time newform__p_coeff_list[j] := Append(newform__p_coeff_list[j], [ Coefficient(form_list[j], p) ]); // Keeps a list of p-powers, starting at p.  

        // Compute the p-power coefficients
        for r in [2 .. Floor(Log(coeff_bound)/Log(p))] do
          //print " Computing the coefficient for newform ", j, " at p^r = ", p, "^", r, " = ", p^r;
          if r eq 2 then
            p_r__coeff := newform__p_coeff_list[j][k][1]^2 - eps(p) * p;
          else
            p_r__coeff := newform__p_coeff_list[j][k][1] * newform__p_coeff_list[j][k][r-1] - eps(p) * p * newform__p_coeff_list[j][k][r-2];
          end if;

          // SANITY CHECK -- Check that our computation agrees with MAGMA. =)
          /*
          print " The relevant prime-power vector is: ", newform__p_coeff_list[j][k];
          print "a(p)^2 = ", newform__p_coeff_list[j][k][1]^2;
          print "eps(p) = ", eps(p);
          print " I computed ", p_r__coeff;
          print " Magma got  ", Coefficient(form_list[j], p^r);
          if p^r lt 100 then
            assert( p_r__coeff eq Coefficient(form_list[j], p^r) );
          end if;
         */

          newform__p_coeff_list[j][k] := Append(newform__p_coeff_list[j][k], p_r__coeff); 
        end for;
    
      end for;

    end if;


    // TO DO: Then check if the prime extends the current coefficient field at all.


    // Make the additional set of coefficients we need to test  (They're m1, not d*m1) 
    possible_prime_powers := [p^i : i in [1 .. Floor(Log(coeff_bound/d)/Log(p))]];   // Make p-powers up to coeff_bound 
    new_coeff_seq := [i*j : i in old_eligible_coeff_seq, j in possible_prime_powers | i*j lt coeff_bound/d ];
    if old_eligible_coeff_seq eq [] then
      new_coeff_seq := [ 1 ] cat new_coeff_seq;
      old_eligible_coeff_seq := new_coeff_seq;
    else 
      old_eligible_coeff_seq := Sort(old_eligible_coeff_seq cat new_coeff_seq);
    end if;



    // Add these new rows to the Big_Eigenrows_Matrix:
    // -----------------------------------------------
    print "About to compute the rows at the coefficients: ", new_coeff_seq;

    // Loop through all of our possible coefficients to make the Big_Eigenrows_Matrix
    for m1 in new_coeff_seq do
      mm := d * m1;

      // Make the mm-th row vector
      s := Cputime();
      tmp_coeff_seq := [];
      for j in allowed_j_list do

        // DIAGNOSTIC
        print "d = ", d, "   mm = ", mm, "    m1 = ", m1, "    possible_primes = ", possible_primes_500[1..k];

        // Compute the m1 coefficient of the j-th newform
        if m1 eq 1 then
          new_tmp_coeff :=  CoefficientField(form_list[j]) ! 1;
        else    
          new_tmp_coeff := &*[newform__p_coeff_list[j][Index(possible_primes_500, q)][Valuation(m1, q)] : q in possible_primes_500 | Valuation(m1, q) ge 1];
        end if;

        /*
        // SANITY CHECK -- Check that this agrees with MAGMA -- THIS TAKES TOO LONG! =(
        if (m1 lt 100) then
          print "Doing sanity check for the j = ", j, " form at the coefficient ", m1;
          assert(new_tmp_coeff eq Coefficient(form_list[j], m1));
        end if;
        */

        tmp_coeff_seq := tmp_coeff_seq cat ElementToSequence(new_tmp_coeff);
      end for;

      tmp_rowvec := RMatrixSpace(RationalField(), 1, #tmp_coeff_seq) ! tmp_coeff_seq;
      print "   TIMING: Took ", Cputime(s), " seconds to compute the eigenvalues at the coefficient ", mm;


      // Add it to the big row matrix, and to the list of (good_mm) row indices
      if m1 eq 1 then
        Big_EigenRows_Matrix := tmp_rowvec;
      else
        Big_EigenRows_Matrix := VerticalJoin(Big_EigenRows_Matrix, tmp_rowvec);
      end if;
      good_mm_list := good_mm_list cat [mm];


    end for;




    // Compute the rank modulo three random primes p > 10000, and remove rows which don't contribute to these ranks.
    // -------------------------------------------------------------------------------------------------------------
    s := Cputime();
    old_size := NumberOfRows(Big_EigenRows_Matrix);

    // Make a list of random primes
    rand_prime_list := [];
    while #rand_prime_list lt 3 do
    rp := RandomPrime(20);
      if rp gt 10000 then
        rand_prime_list := Append(rand_prime_list, rp);
      end if;
    end while;


    // Run through our random primes to do the rank checking
    for rp in rand_prime_list do

      // Find the bad rows in the mod p matrix
      BEMp := ChangeRing(Big_EigenRows_Matrix, FiniteField(rp));
      tmp_echelon := EchelonForm(Transpose(BEMp));  // The columns with pivots index the rows to keep! =)
      good_indices := [];
      for i in [1..NumberOfRows(tmp_echelon)] do
        pivot_flag := false;
        for j in [1..NumberOfColumns(tmp_echelon)] do
          if (tmp_echelon[i,j] ne 0) and (pivot_flag eq false) then
            pivot_flag := true;
            good_indices := Append(good_indices, j);
          end if;
        end for;
      end for;

      // Make a new Big_EigenRows_Matrix and good_mm_list using the good indices above
      BEM2 := Matrix([Big_EigenRows_Matrix[i] : i in good_indices]);
      good_mm2 := [good_mm_list[i] : i in good_indices];
      Big_EigenRows_Matrix := BEM2;
      good_mm_list := good_mm2;
      delete BEM2, good_mm2;

      new_size := NumberOfRows(Big_EigenRows_Matrix);
      print "   TIMING: Took ", Cputime(s), " seconds to pull out", new_size, " rows from the original ", old_size, "ones using the random prime", rp;

    end for;


    // Check if we're done! =)
    if Nrows(Big_EigenRows_Matrix) eq Ncols(Big_EigenRows_Matrix) then
      rank_is_full_flag := true;
    end if;


    // Status report and increment
    print " After using ", k, " primes, we have rank >= ", Nrows(Big_EigenRows_Matrix), 
          " with a goal of ", Ncols(Big_EigenRows_Matrix);
    k := k + 1;
    print "";
    print "";


  end while;
  k:= k - 1;



  // Print the list of coefficients for this d
  print "  For d = ", d, " using the ", #good_mm_list, "Fourier coefficients:";
  print "      ", good_mm_list;

  // Append these coefficients to the list (for each d)
  good_mm_list_vec := Append(good_mm_list_vec, good_mm_list);
  print "good_mm_list_vec = ", good_mm_list_vec;

  // Append these matrices to the list (for each d)
  Big_EigenRows_Matrix_list := Append(Big_EigenRows_Matrix_list, Big_EigenRows_Matrix);
  print "Big_EigenRows_Matrix_list has ", #Big_EigenRows_Matrix_list, " elements.";

  // TIMING: Check the time taken up to this point, and reset the time.
  print "TIMING: Time to compute the d = ", d, " Big_EigenRows_matrix and good_mm_list = ", Cputime(r), " seconds.";
  print "";
  print "";


end for;



// Make a list of primes we actually use  (THIS IS SUPERFLUOUS, BUT USED LATER!)
prime_seq := possible_primes_500[1..#newform__p_coeff_list[1]];
print " prime_seq = ", prime_seq;


//  coeff_seq := Sort(old_eligible_coeff_seq);
//  print "We have used ", #prime_seq, "primes to find ", #coeff_seq, "(possible) Fourier coefficients, ";
//  print "which we hope is enough to determine the ", Dim_S2, "dimensional cuspidal space! =)";
//  print " Eligible coeff_seq = ", coeff_seq;



// TIMING: Check the time taken up to this point, and reset the time.
print "TIMING: Time to compute ALL Big_EigenRows_matrix and good_mm_list = ", Cputime(t), " seconds.";
t := Cputime();



print "";
print "";
print "     =====================================================================";
print "     ===== Computing the prime-power coefficients dividing the level =====";
print "     =====================================================================";
print "";
print "";


// Compute the prime-power coefficients at p|N  (where N is the big level):
// ------------------------------------------------------------------------
D_LCM := LCM(oldform_d_list);
N_prime_seq := PrimeDivisors(D_LCM);
N_newform__p_coeff_list := [ [* *] : j in [1.. #form_list] ];   // one list for each newform, then the primes for each are stored according to prime_seq. =)
for p in N_prime_seq do
  k := Index(N_prime_seq, p);     // Index of p in prime_seq
  for j in [1.. #form_list] do

    // Useful definitions
    N_j := level_list[j];
    bad_primes_j := PrimeDivisors(N_j);


    // Compute the p-th eigenvalue
    print "Computing the coefficient at ", p, "for the newform at j = ", j;
    time N_newform__p_coeff_list[j] := Append(N_newform__p_coeff_list[j], [ Coefficient(form_list[j], p) ]); // Keeps a list of p-powers, starting at p.  

    // Compute the p-power coefficients
    for r in [2 .. Valuation(D_LCM,p)] do
      //print " Computing the coefficient for newform ", j, " at p^r = ", p, "^", r, " = ", p^r;
      if p in bad_primes_j then
        p_r__coeff := N_newform__p_coeff_list[j][k][1] ^ r;
      else
        if r eq 2 then
          p_r__coeff := N_newform__p_coeff_list[j][k][1]^2 - eps(p) * p;
        else
          p_r__coeff := N_newform__p_coeff_list[j][k][1] * N_newform__p_coeff_list[j][k][r-1] - eps(p) * p * N_newform__p_coeff_list[j][k][r-2];
        end if;
      end if;

      // SANITY CHECK -- Check that our computation agrees with MAGMA. =)
      /*
      print " The relevant prime-power vector is: ", newform__p_coeff_list[j][k];
      print "a(p)^2 = ", newform__p_coeff_list[j][k][1]^2;
      print "eps(p) = ", eps(p);
      print " I computed ", p_r__coeff;
      print " Magma got  ", Coefficient(form_list[j], p^r);
      if p^r lt 100 then
        print "Doing sanity check for the j = ", j, " form at the coefficient ", p^r;
        assert( p_r__coeff eq Coefficient(form_list[j], p^r) );
      end if;
     */


      // Append it to the list
      N_newform__p_coeff_list[j][k] := Append(N_newform__p_coeff_list[j][k], p_r__coeff); 
    end for;

  end for;
end for;
//print "The table of newform coefficients at p|N is: ", newform__p_coeff_list;



print "";
print "";
print "     ====================================================================";
print "     ===== Computing the Theta series and its cuspidal coefficients =====";
print "     ====================================================================";
print "";
print "";



// Make the theta series and its Eisenstein and cuspidal parts
// -----------------------------------------------------------

// Compute the precision from the good_mm_list_vec
precision := 0;
for good_mm_list in good_mm_list_vec do
  if not IsEmpty(good_mm_list) then
    precision := Max(precision, Max(good_mm_list));
  end if;
end for;
precision := precision + 1;
print " Using precision ", precision, "(computed from the good_mm_lists)";

// Compute the set of allowed Fourier coefficients
Allowed_coefficients_set := {};
for good_mm_list in good_mm_list_vec do
  Allowed_coefficients_set := Allowed_coefficients_set join Seqset(good_mm_list);
end for;
Allowed_coefficients_list := Setseq(Allowed_coefficients_set);


// Make the theta function
ThetaRing<q>:=PowerSeriesRing(Rationals());
Q2 := Z4 ! QQ;
L2 := LatticeWithGram(Q2);
Theta2 := ThetaRing ! ThetaSeries(L2, 2*precision+1);
Theta_vec := [Coefficient(Theta2, 2*m) : m in [0..precision]];

// Make the Eisenstein series
Gen := GenusRepresentatives(L2);
AvgTheta := &+[ThetaRing ! ThetaSeries(L, 2*precision+1) / #AutomorphismGroup(L) : L in Gen];
AvgAuto := &+[ 1 / #AutomorphismGroup(L) : L in Gen];
Eis2 := AvgTheta / AvgAuto;
Eis_vec := [Coefficient(Eis2, 2*m) : m in [0..precision]];

// Make the Cusp form, and set it as the cuspidal remainder.
Cusp_vec := [Coefficient(Theta2, 2*m) - Coefficient(Eis2, 2*m) : m in [0..precision]];
PS<q> := PowerSeriesRing(Rationals());
Cuspidal_Theta := (PS ! Cusp_vec) + O(q^#Cusp_vec);
remaining_cuspidal_theta := Cuspidal_Theta;

// DIAGNOSTIC -- Print the first 20 terms of Theta(z), E(z), and f(z).
print "";
print "Theta(z) =", (PS ! Theta_vec) + O(q^20);
print "E(z) =", (PS ! Eis_vec) + O(q^20);
print "f(z) =", (PS ! Cusp_vec) + O(q^20);
print "";

// Kill all coefficients of Theta(z) not in the allowed set (so the sum is zero at the end, since we only update these!).
for mm in (Seqset([1..precision]) diff Allowed_coefficients_set) do
  remaining_cuspidal_theta := remaining_cuspidal_theta - Coefficient(remaining_cuspidal_theta,mm) * q^mm;
end for;
print "The truncated theta function (at allowed coefficients) is: ", remaining_cuspidal_theta;


// TIMING: Check the time taken up to this point, and reset the time.
print "TIMING: Time to compute the Eis + Cusp decomposition of Theta(z) = ", Cputime(t), " seconds.";
t := Cputime();



print "";
print "";
print "     ====================================================================";
print "     ===== Computing the new/oldform decomposition of the cusp form =====";
print "     ====================================================================";



// Find the new/oldform decomposition for each possible oldform lift f(z) -> f(dz):
// --------------------------------------------------------------------------------
oldform_coefficient_list := [* *];

//for d in oldform_d_list do
for ind in [1..#oldform_d_list] do

  // Some precomputed info =)
  d := oldform_d_list[ind];
  allowed_j_list := allowed_j_list_vec[ind];
  Big_EigenRows_Matrix := Big_EigenRows_Matrix_list[ind];   // NOTE: THIS MAY SLOW IS DOWN SINCE WE NEED TO COPY THE MATRIX AGAIN!?!
  m_range := good_mm_list_vec[ind];
  m_range_set := Seqset(m_range);

  // TIMING: Reset the time.
  t := Cputime();


  print "";
  print "";
  print "Computing the oldforms for d =", d;
  print "=================================";



  // Compute the dimension of the new/oldform space for this d
  d_dim := &+[Degree(form_list[j]) : j in allowed_j_list];
  print "Using d =", d, "generates a", d_dim, "cuspidal subspace.";


  // Make the direct sum of the trace matrices (only for allowed indices j) 
  trace_matrix_list := [* *];
  for j in allowed_j_list do
    f := form_list[j];

    // Use the basis [1] for the rationals, and the powers of a generator otherwise...
    deg := Degree(CoefficientField(f));
    if (deg eq 1) then  
      temp_trace_matrix := MatrixRing(RationalField(), 1) ! [1];
    else
      temp_trace_matrix := ZeroMatrix(RationalField(), deg, deg);  
      K<a> := CoefficientField(f);

      // Make the trace matrix efficiently. =)
      tmp_a_pow := K ! 1;
      for i_plus_j in [2..2*deg] do

        // Compute the trace
        tmp_a_trace := Trace(tmp_a_pow);

        // Fill in the reverse diagonals
        for i in [Max(1, i_plus_j - deg) .. Min(deg, i_plus_j - 1)] do
          j1 := i_plus_j - i;
          temp_trace_matrix[i,j1] := tmp_a_trace;
        end for;

        // Increment the element
        tmp_a_pow := tmp_a_pow * a;

      end for;
//Old Way:   temp_trace_matrix := MatrixRing(RationalField(), deg) ! [Trace(a^((i-1)+(j-1))): i, j in [1..deg] | i le j];
    end if;


    trace_matrix_list := Append(trace_matrix_list, temp_trace_matrix);
  end for;

  Big_Trace_Matrix := trace_matrix_list[1];
  for j in [2..#trace_matrix_list] do
    Big_Trace_Matrix := DirectSum(Big_Trace_Matrix ,trace_matrix_list[j]);
  end for;  
  print "  The trace matrix has", Nrows(Big_Trace_Matrix), "rows and ", Ncols(Big_Trace_Matrix), "columns.";

  // TIMING: Check the time taken up to this point, and reset the time.
  print "TIMING: Time to compute the trace matrix = ", Cputime(t), " seconds.";
  t := Cputime();



/*
  // Define the m-coefficients to consider (here we only care about newforms)
  N_min := Conductor(eps);
  further_lift_num := N div (d * N_min);
  m_range := [d*i : i in [1..(precision div d)] | GCD(i, further_lift_num) eq 1];      // TO DO: Need to adjust this, and check the rank is maximal!
//  m_range := [d*i : i in [1..4*d_dim] | GCD(i, further_lift_num) eq 1];      // TO DO: Need to adjust this, and check the rank is maximal!
  print "  Using", #m_range, "Fourier coefficients to deduce the new/oldforms!";


  // Make a list of vectors for the m-th eigenform coefficients (as we vary m).
  for mm in m_range  do

    // Make the m-th row vector
    tmp_coeff_seq := [];
    for j in allowed_j_list do
        tmp_coeff_seq := tmp_coeff_seq cat ElementToSequence(Coefficient(form_list[j], mm div d));
    end for;
    tmp_rowvec := RMatrixSpace(RationalField(), 1, #tmp_coeff_seq) ! tmp_coeff_seq;

    // Add it to the big row matrix
    if mm eq d then
      Big_EigenRows_Matrix := tmp_rowvec;
    else
      Big_EigenRows_Matrix := VerticalJoin(Big_EigenRows_Matrix, tmp_rowvec);
    end if;

  end for;

  // TIMING: Check the time taken up to this point, and reset the time.
  print "TIMING: Time to compute the eigenform matrix = ", Cputime(t), " seconds.";
  t := Cputime();
*/

  


  // Make the matrix to solve for the column vector a_j
  A := Big_EigenRows_Matrix * Big_Trace_Matrix;
  print "  The total matrix A is", Nrows(A), "x", Ncols(A), ", and has rank", Rank(A);  

  // DIAGNOSTIC
  if not (Ncols(A) eq Rank(A)) then
    print "  ERROR:  THE RANK SEEMS TO BE TOO SMALL..." ;
    print "  Big_Trace_Matrix =", Big_Trace_Matrix;
    print "  Big_EigenRows_Matrix =", Big_EigenRows_Matrix;
    print "  A = ", A;
  end if;

  assert(Ncols(A) eq Rank(A));


  // TIMING: Check the time taken up to this point, and reset the time.
  print "TIMING: Time to compute the product matrix (eigenform * trace) = ", Cputime(t), " seconds.";
  t := Cputime();



  // Make the cusp form vector for the Fourier coefficients in m_range:
  // ------------------------------------------------------------------
  //new_coeff_vec := RMatrixSpace(RationalField(), #m_range, 1) ! [Cusp_vec[m+1] : m in m_range];
  /*
  N_min := Conductor(eps);
  further_lift_num := N div (d * N_min);
  m_range := [d*i : i in [1..(precision div d)] | GCD(i, further_lift_num) eq 1];      // TO DO: Need to adjust this, and check the rank is maximal!
  */
  new_coeff_vec := RMatrixSpace(RationalField(), #m_range, 1) ! [Coefficient(remaining_cuspidal_theta, m) : m in m_range];
  b := new_coeff_vec;


  // Solve the system!
  // -----------------
  ok, x := IsConsistent(Transpose(A), Transpose(b));
  if not ok then
    error "The system is not consistent... =(";
  end if;
  x_seq := ElementToSequence(x);
  print "  Is the system consistent?", ok;
  assert(ok eq true);

  // TIMING: Check the time taken up to this point, and reset the time.
  print "TIMING: Time to solve the system = ", Cputime(t), " seconds.";
  t := Cputime();



  // Break x up into coefficients for each Galois conjugacy class
  // -------------------------------------------------------------
  other_x := x_seq;
  coefficient_list := [* *];
  for j in [1..#form_list] do
    f := form_list[j];
    K := CoefficientField(f);
    deg := Degree(f);

    // Write the coefficients we found, with zeros elsewhere.
    if j in allowed_j_list then
      first_x := [other_x[i] : i in [1..deg]];
      other_x := [other_x[i] : i in [deg+1..#other_x]];
      coefficient_list := Append(coefficient_list, K ! first_x);
    else
      coefficient_list := Append(coefficient_list, K ! [0 : i in [1..deg]]);
    end if;

  end for;

  // Append the current coefficient_list to the cumulative oldform_coefficient_list
  oldform_coefficient_list := Append(oldform_coefficient_list, coefficient_list);
  print "  We have the coefficient vector", coefficient_list;


  // TIMING: Check the time taken up to this point, and reset the time.
  print "TIMING: Time to organize the coefficients of the cuspidal part = ", Cputime(t), " seconds.";
  t := Cputime();


  // Find what's left over when we subtract the old/newform part for this d-value
  // -----------------------------------------------------------------------------

/*
  // Make our newform component as a power series
  new_component_vec := [ 0 ];
  for mm in [1..precision] do
    tmp_num := 0;
    if (mm mod d eq 0) then
      for j in [1..#form_list] do

        // DIAGNOSTIC
        //print "Precision = ", precision;
        //time CoefficientField(form_list[j]);
    
        // We need to be careful since Magma doesn't deal well with AbsoluteTrace of a rational number! =(
        if CoefficientField(form_list[j]) eq RationalField() then 
          tmp_num := tmp_num + coefficient_list[j] * Coefficient(form_list[j], mm div d);
        else 
          tmp_num := tmp_num + AbsoluteTrace(coefficient_list[j] * Coefficient(form_list[j], mm div d));
        end if;
    
      end for;
    end if;
    new_component_vec := new_component_vec cat [tmp_num];
  end for;
*/


  // Make the big matrix used to find any rational coefficient
  BB := x * Big_Trace_Matrix;

  // Make our newform component as a power series
  new_component_vec := [ 0 ];

print "";
print "Subtracting the d =", d, "Fourier contribution";
print "-------------------------------------------";
print "m_range_set = ", m_range_set;

  // MAIN LOOP for subtracting the newform components
  for mm in [1..precision] do
    tmp_num_seq := [0];

    // Only compute coefficients for allowed entries! =)
    if (mm mod d eq 0) and (mm in Allowed_coefficients_set) then

      // Don't recompute the answer for elements of M_d, instead steal them to make the answer zero. =)
      if mm in m_range_set then
        tmp_num_seq := [ Coefficient(remaining_cuspidal_theta, mm) ];   // Using mm+1 since the 1st entry is the 0-th Fourier coefficient!
      else
 
        // TIMING: Set the time
        s := Cputime();

        // Make the rational eigenform coefficient (column) vector for each coefficient mm
        tmp_coeff_seq := [];
        for j in allowed_j_list do

          // Compute the m1 coefficient of the j-th newform -- NOW INCLUDING PRIMES AT EACH LEVEL! =0
	  // ----------------------------------------------------------------------------------------
  	  m1 := mm div d;
          m1_bad := GCD(m1, N);
          m1_good := m1 div m1_bad;
          assert (GCD(N, m1_good) eq 1);   // Sanity check! =)

          // Compute the part away from N
          if m1_good eq 1 then
            new_tmp_coeff :=  CoefficientField(form_list[j]) ! 1;
          else    
            new_tmp_coeff := &*[newform__p_coeff_list[j][Index(prime_seq, p)][Valuation(m1_good, p)] : p in prime_seq | Valuation(m1_good, p) ge 1];
          end if;

          // Compute the part at N
          if m1_bad eq 1 then
            new_tmp_coeff := new_tmp_coeff * CoefficientField(form_list[j]) ! 1;
          else    
            new_tmp_coeff := new_tmp_coeff * &*[N_newform__p_coeff_list[j][Index(N_prime_seq, p)][Valuation(m1_bad, p)] : p in N_prime_seq | Valuation(m1_bad, p) ge 1];
          end if;


          /*
          // SANITY CHECK -- See if it agrees with MAGMA''s computation! =)
          if p^r lt 100 then
            print "Doing sanity check for the j = ", j, " form at the coefficient ", p^r;
            assert( p_r__coeff eq Coefficient(form_list[j], p^r) );
          end if;
	  */


          // Append it to the others.
          tmp_coeff_seq := tmp_coeff_seq cat ElementToSequence(new_tmp_coeff);

        end for;
        tmp_colvec := RMatrixSpace(RationalField(), #tmp_coeff_seq, 1) ! tmp_coeff_seq;

        // Make the new_coefficient at mm
        tmp_num_matrix := BB * tmp_colvec;
        print "New coefficient at ", mm, " = ", tmp_num_matrix;
        tmp_num_seq := ElementToSequence(tmp_num_matrix);

        // TIMING: Print the time used
        print "   TIMING: Took ", Cputime(s), " seconds to compute the 'newform' coefficient at ", mm;

      end if;

    end if;
    new_component_vec := new_component_vec cat tmp_num_seq;
  end for;

  // TIMING: Check the time taken up to this point, and reset the time.
  print "TIMING: Time to create the newly found 'newform' piece as a power series = ", Cputime(t), " seconds.";
  t := Cputime();


  // Subtract it from the remaining cuspidal theta
  newform_component := PS ! new_component_vec;
  remaining_cuspidal_theta := remaining_cuspidal_theta - newform_component;
  print "  The remaining cuspidal component after removing the oldforms at", d, " is:";
  print remaining_cuspidal_theta;

  // TIMING: Check the time taken up to this point, and reset the time.
  print "TIMING: Time to subtract off this cuspidal part = ", Cputime(t), " seconds.";
  t := Cputime();

end for;




print "";
print "";
print "     ===============================================";
print "     ===== Printing the cuspidal coefficients  =====";
print "     ===============================================";
print "";
print "";



// Print the newforms, their levels, and their coefficients
print "The level list is:", level_list;
print "The newform list is:", form_list;
print "";
print "The old/newform lift are",  oldform_d_list;
print "The old/newform coefficients are given by the list:", oldform_coefficient_list;



print "";
print "";
print "     ===============================================";
print "     ===== Computing the overall cusp constant =====";
print "     ===============================================";
print "";


// Compute the overall cuspidal constant! =)
total_cusp_constant := 0;

for ind in [1..#oldform_d_list] do

  oldform_coefficient_list_by_d := oldform_coefficient_list[ind];
  d := oldform_d_list[ind];

  print "";
  print "Coefficients for d =", d;
  print "=======================";


  for c in oldform_coefficient_list_by_d do


    r1, r2 := Signature(Parent(c));    // Finds the number of real and complex embeddings
                                       // Note: This should be the same as the embeddings of the newform, 
                                       //       even if the coefficient is in a smaller number field, 
                                       //       since we coerced it into the same ring.  EXPLICITLY CHECK THIS!!!

    // Make the absolute values
    if (r1 eq 1) and (r2 eq 0) then      // Deal with the rational field separately since Magma is unhappy... =(
      c_embedding_sizes := [Abs(c)];
    else
      c_embedding_sizes := AbsoluteValues(c);
    end if;
    c_real := c_embedding_sizes[1..r1];
    c_complex := c_embedding_sizes[r1+1..r2];

    // Compute the sum of all absolute values (with multiplicity)
    c_cusp_const := &+c_real;
    c_cusp_const := c_cusp_const + 2 * &+c_complex;
 
    // Print the computed value, and add it to the overall cuspidal constant
    print "We computed the constant for c =", c, "as", c_cusp_const + 0.00000;
    print "";
    total_cusp_constant := total_cusp_constant + c_cusp_const;

  end for;

end for;

// TIMING: Check the time taken up to this point, and reset the time.
print "TIMING: Time to compute the real (embedded) cuspidal constant = ", Cputime(t), " seconds.";
t := Cputime();


print "\n This gives the overall cuspidal constant of:", total_cusp_constant + 0.00000;



// Print MAGMA Info 
print "";
print "==============================================";
print "MAGMA Info (in cusp_const() routine):";
print "-------------------------------------";
print "MAGMA Version:", GetVersion();
print "MAGMA Current Memory Usage:", GetMemoryUsage();
print "";






return total_cusp_constant;

end function;


/////////////////////////////////////////////////////////////////////////////////
// TEST CODE FOR UNDERSTANDING THE OUTPUT OF AbsoluteValues():
// -----------------------------------------------------------
// NN<a> := NumberField(Polynomial([5,0,1]));    // This makes Q(sqrt(-5))
// a^2;
// AbsoluteValues(a);
/////////////////////////////////////////////////////////////////////////////////




/*
// This finds the first n decimal digits of the number x
function chop(x, n)
     return Truncate(x * 10^n) / (10^n);
end function;
*/


