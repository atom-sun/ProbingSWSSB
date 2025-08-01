(* ::Package:: *)

(* Mathematica Package *)
(* :Title: Simulation of Randomized Measurement *)
(* :Context: Simulation of Randomized Measurement to measure Renyi2. *)
(* :Author: Ning Sun *)
(* :Date: 2024-12-10 *)
(* :Package Version: 1.0 *)
(* :Mathematica Version: 14.0 *)
(* :Copyright:(c) 2024 atom-sun *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["SimRandomized`"];
(* Exported symbols added here with SymbolName::usage *)

measureState::usage = "
Measure a state with respect to specified uchoices (x/y/z).
";

measureStatesX::usage = "
Randomized measurements of X.
"


Begin["`Private`"];

On[Assert];

<< "SpinOneHalfChain.wl";
<< "OpenSystem.wl";
<< "Measurement.wl";

(* define: Project Operators on a single site. *)
xProjUp = (sI + sX)/2;
xProjDown = (sI - sX)/2;
yProjUp = (sI + sY)/2;
yProjDown = (sI - sY)/2;
zProjUp = (sI + sZ)/2;
zProjDown = (sI - sZ)/2;

(* XYZ projection operators on a nSites spin-1/2 chain of the jSite. *)
projXYZOps[nSites_Integer, jSite_Integer] := Module[
   {xup, xdown, yup, ydown, zup, zdown, op},
   
   Assert[jSite <= nSites];
   
   op = (TensorProduct[IdentityMatrix[2^(jSite - 1)], #, 
         IdentityMatrix[2^(nSites - jSite)]] // ArrayFlatten // 
       ArrayFlatten) &;
   {xup, xdown, yup, ydown, zup, zdown} =
    op /@ {xProjUp, xProjDown, yProjUp, yProjDown, zProjUp, zProjDown};
   
   (*{{xup,xdown},{yup,ydown},{zup,zdown}}*)
   <|"X" -> {xup, xdown}, "Y" -> {yup, ydown}, "Z" -> {zup, zdown}|>
   ];

(* measure state under a set of uchoices. 
state can be a vector of a pure state or a density matrix. 
*)
measureState[state_, uchoices_] := Module[
   {nSites, newstate, s, u, pr0, pr1, p0, p1, s0, r, pr},

   << "Measurement.wl";
   
   nSites = Length[uchoices];
   Assert[Length[state] == 2^nSites];
   newstate = state;
   
   s = Table[
     u = uchoices[[jStep]];
     {pr0, pr1} = projXYZOps[nSites, jStep][u];
     
     p0 = calcProb[newstate, pr0];
     p1 = calcProb[newstate, pr1];
     Assert[Abs[p0 + p1 - 1] < 10.^(-8)];
     
     r = RandomReal[];
     s0 = Which[r <= p0, pr = pr0; 0,
       r > p0, pr = pr1; 1];
     newstate = updateState[newstate, pr];
     
     s0
     , {jStep, nSites}
     ];
   
   StringRiffle[s, ""]
   ];

(* measure numerator of Renyi2 correlator of state and state2 by \
randomized measurements. 
state can be a vector of a pure state or a density matrix. 
*)
(* Tr(rho ZiZj rho ZiZj) *)
(* return all measurement results xs. *)
measureStatesX[state_, state2_, nSamplesU_Integer, 
   mSamplesS_Integer] := Module[
   {dHilbert, nSites, uchoices, s, s1, ds, ss, xs},
   
   Assert[Length[state] == Length[state2]];
   
   (* constants: *)
   dHilbert = 2;
   nSites = Log2[Length[state]];
   
   ss = Table[
     uchoices = Table[RandomChoice[{"X", "Y", "Z"}], nSites];
     Table[
      s = measureState[state, uchoices];
      s1 = measureState[state2, uchoices];
      ds = HammingDistance[s, s1] // N;
      ds
      , {jSampleS, mSamplesS}]
     , {jSampleU, nSamplesU}];
   
   xs = (-dHilbert)^(-ss);
   xs *= dHilbert^nSites;
   xs
   ];


End[];
EndPackage[];
