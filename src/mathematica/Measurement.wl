(* ::Package:: *)

(* Mathematica Package *)
(* :Title: Measurement Calculation *)
(* :Context: Measurement calculation. *)
(* :Author: Ning Sun *)
(* :Date:2024-11-30 *)
(* :Package Version:1.0 *)
(* :Mathematica Version:14.0 *)
(* :Copyright:(c) 2024 atom-sun *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["Measurement`"];
(* Exported symbols added here with SymbolName::usage *)

calcProb::usage = "
Compute the probability under projection operator with respect to certain state.
state can be a vector of a pure state or a density matrix. 

calcProb[state, projop]

";

updateState::usage = "
Update state with respect to projection operator. 
state can be a vector of a pure state or a density matrix. 

updateState[state, projop]

"


Begin["`Private`"];

On[Assert];

calcProb[state_, projop_] /; Dimensions[Dimensions[state]] == {1} :=
  Chop[Conjugate[state] . projop . state];
calcProb[state_, projop_] /; Dimensions[Dimensions[state]] == {2} :=
  ((* Assert[Transpose[Conjugate[state]] == state]; *)
   Chop[Tr[state . projop]]);

updateState[state_, projop_] :=
  Module[{newstate},
   If[Dimensions[Dimensions[state]] == {1},
    
    newstate = projop . state;
    newstate /= Sqrt[(Conjugate[newstate] . newstate)];,
    
    Assert[Dimensions[Dimensions[state]] == {2}];
    newstate = Transpose[Conjugate[projop]] . state . projop;
    newstate /= Tr[newstate];
    ];
   
   newstate
   ];

End[];
EndPackage[];
