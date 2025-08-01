(* ::Package:: *)

(* Mathematica Package *)
(* :Title: General Implementation on a Spin-1/2 Chain. *)
(* :Context: General Implementation on a Spin-1/2 Chain. *)
(* :Author: Ning Sun *)
(* :Date: 2024-11-24 *)
(* :Package Version: 1.0 *)
(* :Mathematica Version: 14.0 *)
(* :Copyright:(c) 2024 atom-sun *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["SpinOneHalfChain`"];
(* Exported symbols added here with SymbolName::usage *)

sI::usage = "Identity matrix of dimension 2: {{1, 0}, {0, 1}}. "
sX::usage = "Pauli X matrix: {{0, 1}, {1, 0}}. "
sY::usage = "Pauli Y matrix: {{0, -I}, {I, 0}}. "
sZ::usage = "Pauli Z matrix: {{1, 0}, {0, -1}}. "

anyOp::usage = "
Any operator defined on a spin-1/2 chain.
";

generalOp::usage = "
General operator defined on a spin-1/2 chain.
";

zizjOp::usage = "ZiZj Operator. ";
zjOp::usage = "Zj Operator. ";

reduce1dm::usage = "Reduce density matrix reduced over last one site. ";
reduceDensmat::usage = "Reduce density matrix over last j sites. ";

svonnEntropy::usage = "Calculate von Neumann Entropy. ";
srenyi2Entropy::usage = "Calculate Renyi2 Entropy. ";


Begin["`Private`"];

On[Assert];

(* constants. *)
sI = IdentityMatrix[2];
sX = PauliMatrix[1];
sY = PauliMatrix[2];
sZ = PauliMatrix[3];

(* define any operators on the spin-1/2 chain. *)
anyOp[nSites_Integer, ops__String] := 
  Module[{mapOps, op}, 
   Assert[Length[{ops}] == nSites];
   mapOps = <|"I" -> sI, "X" -> sX, "Y" -> sY, "Z" -> sZ|>;
   Assert[AllTrue[{ops}, MemberQ[Keys[mapOps], #] &]];
   op = Nest[ArrayFlatten, 
     TensorProduct @@ Table[mapOps[j], {j, {ops}}], nSites - 1];
   op];

(* define general operators on a spin-1/2 chain. *)
generalOp[nSites_Integer, ops__] := Module[
   {l, nSymbs, sites, symbs, where, mapOps, symb, oos, idims, imats, 
    op},
   
   l = Length[{ops}];
   Assert[Mod[l, 2] == 0];
   nSymbs = l/2;
   Assert[nSymbs <= nSites];
   
   sites = {ops}[[1 ;; l ;; 2]];
   Assert[Max[Last[Transpose[Tally[sites]]]] == 1];
   symbs = {ops}[[2 ;; l ;; 2]];
   where = 
    DeleteDuplicates[
     Flatten[Map[Position[sites, #] &, Sort[sites, Less]]]];
   sites = sites[[where]];
   symbs = symbs[[where]];
   
   mapOps = <|"I" -> sI, "X" -> sX, "Y" -> sY, "Z" -> sZ|>;
   
   oos = Table[If[symb = symbs[[i]]; StringQ[symb],
      Assert[MemberQ[Keys[mapOps], symb]]; mapOps[symb],
      Assert[Dimensions[symb] == {2, 2}]; symb
      ]
     , {i, nSymbs}];
   idims = 
    Join[{First[sites] - 1}, 
     Differences[sites] - 1, {nSites - Last[sites]}];
   imats = Table[IdentityMatrix[2^i], {i, idims}];
   oos = Riffle[imats, oos];
   op = Nest[ArrayFlatten, TensorProduct @@ oos, Length[oos] - 1];
   op
   ];

(* test: compare with anyOp. *)
generalOp[4, 1, "I", 2, "X", 3, "Y", 4, "Z"] == 
 anyOp[4, "I", "X", "Y", "Z"]


(* define correlation functions.  <ZiZj> *)
zizjOp[nSites_Integer, i_Integer, j_Integer] /; i < j := 
  TensorProduct[IdentityMatrix[2^(i - 1)], sZ, 
       IdentityMatrix[2^(j - i - 1)], sZ, 
       IdentityMatrix[2^(nSites - j)]] // ArrayFlatten // 
     ArrayFlatten // ArrayFlatten // ArrayFlatten;
zizjOp[nSites_Integer, i_Integer, j_Integer] /; i > j := 
  zizjOp[nSites, j, i];
zizjOp[nSites_Integer, i_Integer, j_Integer] /; i == j := 
  TensorProduct[IdentityMatrix[2^(j - 1)], sZ.sZ, 
     IdentityMatrix[2^(nSites - j)]] // ArrayFlatten // ArrayFlatten;

(* define order parameter. <Zj> *)
zjOp[nSites_Integer, j_Integer] := 
  TensorProduct[IdentityMatrix[2^(j - 1)], sZ, 
     IdentityMatrix[2^(nSites - j)]] // ArrayFlatten // ArrayFlatten;

(* define reduced density matrix *)
(* ResourceFunction["MatrixPartialTrace"] *)
(* partial trace from right side one-by-one. *)
reduce1dm[rho_] := Module[{l, dm1, dm2}, l = Length[rho];
   dm1 = rho[[1 ;; Length[rho] ;; 2, 1 ;; Length[rho] ;; 2]];
   dm2 = rho[[2 ;; Length[rho] ;; 2, 2 ;; Length[rho] ;; 2]];
   dm1 + dm2];
reduceDensmat[rho_, j_Integer] := Nest[reduce1dm, rho, j];
densMat[psi_] := KroneckerProduct[psi, psi];

(* calculate entanglement entropies. *)
(* svonnEntropy[rho_]:=-Tr[rho.MatrixLog[rho]]; *)
(* svonnEntropy[rho_]:=-Tr[rho.Log[rho]]; *)
(* svonnEntropy[rho_]:=-Tr[rho.MatrixLog[rho+10^(-4)*Identity[Length[rho]]]]; *)
svonnEntropy[rho_] := Module[{r}, r = Eigenvalues[rho];
   Chop[-r*Log[r]] // Total];
srenyi2Entropy[rho_] := -Log[Tr[rho . rho]];


End[];

EndPackage[];
