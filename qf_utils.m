// ======================================================
// Utilities for theta-eigenform-decomposition--magma 
// (c) 2011 Jonathan Hanke, William Stein
// Released under GNU Public License (GPL) version 2.0
// ======================================================





// -------------------------------------------------------------------------
// Quadratic Form Utilities -- Level, Theta series, and Eisenstein Series
// -------------------------------------------------------------------------




// Level function -- kindly provided by W. Stein 
function LevelOfQuadForm(B)
     //intrinsic LevelOfQuadForm(A::AlgMatElt) -> RngIntElt
     //{}
     /*
The level of Q is the minimal number N such that N * B^(-1) is a
matrix with integer entries and even diagonal entries; here B is twice
the matrix associated to Q.  So for example, the four-squares form has
corresponding matrix B=diag[2,2,2,2]; its inverse is then [1/2,1/2,1/2,1/2],
so N=4.
     */

//     B := 2*A;  // Oops! We assume that we are given B to begin with! =O
Binv := B^(-1);
N := Determinant(B);
N := N/(2^Valuation(N,2));
l := LCM([Denominator(x) : x in Eltseq(N*Binv)]);
N := N*l;
m := GCD([Numerator(N*Binv[i,i]) : i in [1..4]]);
if IsOdd(m) then
N := N*2;
end if;
return Integers()!N;
//end intrinsic;
end function;




// =========================================================================================
// =========================================================================================



// Makes the Theta series for the form Q (where Q2 = 2*Q) up to precision prec.
// ---------------------------------------------------------------------------- 
function Make_Theta_Series(Q2, prec)

  // Declare some variables
  Z := Integers();
  Q := Rationals();
  n := Nrows(Q2);
  ZZ := RMatrixSpace(Z, n, n);
  QQ := RMatrixSpace(Q, n, n);
  ThetaRing<q>:=PowerSeriesRing(Q);

  // Make the Theta series
  L2 := LatticeWithGram(Q2);
  Theta2 := ThetaRing ! ThetaSeries(L2, 2*prec);

  // Substitute x for x^2 everywhere
  Theta_half := ThetaRing ! [ Coefficient(Theta2, 2*i) : i in [0..Floor((AbsolutePrecision(Theta2) - 1)/2)]];

  // Return the Theta Series
  return Theta_half;

end function;



// =========================================================================================
// =========================================================================================


// Makes the Eisinstein series for the form Q (where Q2 = 2*Q) up to precision prec.
// ---------------------------------------------------------------------------------
function Make_Eisenstein_Series(Q2, prec)

  // Declare some variables
  Z := Integers();
  Q := Rationals();
  n := Nrows(Q2);
  ZZ := RMatrixSpace(Z, n, n);
  QQ := RMatrixSpace(Q, n, n);
  ThetaRing<q>:=PowerSeriesRing(Q);

  // Make the Eisenstein series
  L2 := LatticeWithGram(Q2);
  Gen := GenusRepresentatives(L2);
  AvgTheta := &+[ThetaRing ! ThetaSeries(L, 2*prec) / #AutomorphismGroup(L) : L in Gen];
  AvgAuto := &+[ 1 / #AutomorphismGroup(L) : L in Gen];
  Eis2 := AvgTheta / AvgAuto;

  // Substitute x for x^2 everywhere
  Eis_half := ThetaRing ! [ Coefficient(Eis2, 2*i) : i in [0..Floor((AbsolutePrecision(Eis2) - 1)/2)]];

  // Return the Eisenstein Series
  return Eis_half;

end function;









