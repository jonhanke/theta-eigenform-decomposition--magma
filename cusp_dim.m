// ===================================================================
// Cuspidal dimension routine for theta-eigenform-decomposition--magma 
// (c) 2011 Jonathan Hanke
// Released under GNU Public License (GPL) version 2.0
// ===================================================================





// ----------------------------------------------------------
// Computes the dimension of the full space of cusp forms  
// containing the theta series of the given quadratic form.
// ----------------------------------------------------------

function compute_cusp_dimension(QQ)

     // TIMING: Set the current time
     t := Cputime();

     // DIAGNOSTIC
     print "Using the matrix 2*Q =",QQ;

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
print " ";
print "The form Q has level", N, "and determinant", D/16;


M2 := ModularForms(DirichletGroup(N) ! eps);
S2 := CuspidalSubspace(M2);
print "The dimension of the cuspidal subspace is ", Dimension(S2);
precision := PrecisionBound(S2);
print "Need the precision: ", precision, "to uniquely identify a cusp form.";
print " ";


// TIMING: Check the time taken up to this point, and reset the time.
print "TIMING: Time to setup the spaces of modular forms = ", Cputime(t), " seconds.";
t := Cputime();
print " ";

return Dimension(S2);
end function;
