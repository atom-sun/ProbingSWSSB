(* ::Package:: *)

(* Mathematica Package *)
(* :Title: Simulation of Swap *)
(* :Context: Simulation of Swap to measure Renyi2 entropy and correlators. *)
(* :Author: Ning Sun *)
(* :Date: 2024-12-13 *)
(* :Package Version: 1.0 *)
(* :Mathematica Version: 14.0 *)
(* :Copyright:(c) 2024 atom-sun *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["SimSwap`"];
(* Exported symbols added here with SymbolName::usage *)

projSwapOp::usage = "
Define swap operator on iSite and jSite of nSites qubit system.
";

measureParity::usage = "
One single measurement of whole parity of state and state2.
"


Begin["`Private`"];

On[Assert];

<< "SpinOneHalfChain.wl";
<< "QIsingED.wl";
<< "OpenSystem.wl";
<< "Measurement.wl";


(* define swap operator on iSite and jSite of nSites qubit system. *)
projSwapOp[nSites_Integer, iSite_Integer, jSite_Integer] /; 
   iSite <= jSite :=
  Module[
   {genOp, projuu, projdd, projss, projaa, i, x, y, z},
   
   Assert[iSite < jSite <= nSites];
   Assert[iSite != jSite];
   
   genOp = 
    Nest[ArrayFlatten, 
      TensorProduct[IdentityMatrix[2^(iSite - 1)], #1, 
       IdentityMatrix[2^(jSite - iSite - 1)], #2, 
       IdentityMatrix[2^(nSites - jSite)]], 4] &;
   (*projuu=Nest[ArrayFlatten,TensorProduct[IdentityMatrix[
   2^(iSite-1)],(sI+sZ)/2,IdentityMatrix[2^(jSite-iSite-1)],(sI+sZ)/2,
   IdentityMatrix[2^(nSites-jSite)]],4];
   projdd=Nest[ArrayFlatten,TensorProduct[IdentityMatrix[
   2^(iSite-1)],(sI-sZ)/2,IdentityMatrix[2^(jSite-iSite-1)],(sI-sZ)/2,
   IdentityMatrix[2^(nSites-jSite)]],4];*)
   projuu = genOp @@ {(sI + sZ)/2, (sI + sZ)/2};
   projdd = genOp @@ {(sI - sZ)/2, (sI - sZ)/2};
   {i, x, y, z} = genOp @@@ {{sI, sI}, {sX, sX}, {sY, sY}, {sZ, sZ}};
   projss = i/2 + (i + x + y + z)/4 - projuu - projdd;
   projaa = (i - x - y - z)/4;
   
   (* ev: 1, 1, 1, -1 *)
   {projuu, projdd, projss, projaa}
   
   ];
projSwapOp[nSites_Integer, iSite_Integer, jSite_Integer] /; 
   iSite > jSite := projSwapOp[nSites, jSite, iSite];
  

(* one single measurement of whole parity of state and state2. *)
measureParity[state_, state2_] := Module[
   {doublestate, nSites, nSites2, uu, dd, ss, aa, puu, pdd, pss, paa, 
    s, s0, pcum, r, pr, parity},
   
   Assert[Length[state] == Length[state2]];
   doublestate = TensorProduct[state, state2];
   If[Dimensions[Dimensions[state]] == {1},
    doublestate = TensorProduct[state, state2] // Flatten;,
    Assert[Dimensions[Dimensions[state]] == {2}];
    doublestate = TensorProduct[state, state2] // ArrayFlatten;
    ];
   nSites = Log2[Length[state]];
   nSites2 = nSites*2;
   
   s = Table[
     {uu, dd, ss, aa} = projSwapOp[nSites2, jStep, jStep + nSites];
     {puu, pdd, pss, paa} = 
      calcProb[doublestate, #] & /@ {uu, dd, ss, aa};
     Assert[puu + pdd + pss + paa == 1];
     
     pcum = Accumulate[{puu, pdd, pss, paa}];
     r = RandomReal[];
     s0 = Which[r <= pcum[[1]], pr = uu; 1,
       pcum[[1]] < r <= pcum[[2]], pr = dd; 1,
       pcum[[2]] < r <= pcum[[3]], pr = ss; 1,
       pcum[[3]] < r, pr = aa; -1
       ];
     doublestate = updateState[doublestate, pr];
     
     s0
     , {jStep, nSites}];
   
   parity = FoldList[Times, 1, s][[-1]];
   parity
   ];


End[];
EndPackage[];
