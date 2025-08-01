(* ::Package:: *)

(* Mathematica Package *)
(* :Title: Quantum Ising model Exact Diagonalization *)
(* :Context: Quantum Ising model Exact Diagonalization *)
(* :Author: Ning Sun *)
(* :Date:2024-11-20 *)
(* :Package Version:1.0 *)
(* :Mathematica Version:14.0 *)
(* :Copyright:(c) 2024 atom-sun *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["QIsingED`"];
(* Exported symbols added here with SymbolName::usage *)

qIsingED::usage = "
Quantum (Transverse-field) Ising model Exact Diagonalization overall.
Return ground state energy and with its state vector in the tensor space. 
";

htisingHamilt::usage = "
Quantum Ising model Hamiltonian. 
"


Begin["`Private`"];

On[Assert];

<<"SpinOneHalfChain.wl";

(* define transverse-field Ising model Hamiltonian. *)
htisingHamilt[nSites_Integer, g_, bc_String : "obc"] := 
  Module[{hzz, hx, hb}, 
   hzz = Sum[
     anyOp[nSites, 
      Sequence @@ 
       Join[Table["I", i - 1], Table["Z", 2], 
        Table["I", nSites - 1 - i]]], {i, 1, nSites - 1}];
   If[bc == "pbc", 
    hb = TensorProduct[sZ, IdentityMatrix[2^(nSites - 2)], sZ] // 
       ArrayFlatten // ArrayFlatten;
    hzz += hb];
   hx = Sum[
     anyOp[nSites, 
      Sequence @@ Join[Table["I", i - 1], {"X"}, Table["I", nSites - i]]], 
       {i, 1, nSites}];
   -1.0*(hzz + g*hx)];


qIsingED[nSites_Integer, g_, bc_String : "obc"] := 
  Module[{erg, vec, where, ham, e0, v0}, 

   ham = htisingHamilt[nSites, g, bc];
   {erg, vec} = Eigensystem[N[ham]];
   where = 
    DeleteDuplicates[
     Flatten[Map[Position[erg, #] &, Sort[erg, Less]]]];
   erg = erg[[where]];
   vec = vec[[where]];

   (* ground state energy. *)
   e0 = erg[[1]];
   (* ground state. *)
   v0 = vec[[1]];

   {e0, v0}];

End[];
EndPackage[];