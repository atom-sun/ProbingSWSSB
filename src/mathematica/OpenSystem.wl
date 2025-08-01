(* ::Package:: *)

(* Mathematica Package *)
(* :Title: Open System *)
(* :Context: Open System related implementations. *)
(* :Author: Ning Sun *)
(* :Date: 2024-12-02 *)
(* :Package Version: 1.0 *)
(* :Mathematica Version: 14.0 *)
(* :Copyright:(c) 2024 atom-sun *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["OpenSystem`"];
(* Exported symbols added here with SymbolName::usage *)

kkrausChannelZZ::usage = "Apply ZZ-type Kraus quantum channel. ";

zzChannelA2Adm::usage = "
Apply all-to-all ZZ quantum channel to a density matrix.
";

zzChannelA2Astate::usage = "
Probabilistically (classically) apply all-to-all ZZ gate to a state.
";


Begin["`Private`"];

On[Assert];

<<"SpinOneHalfChain.wl";


(* kkraus zz quantum channel. *)
(* KKrause quantum channel ZZ is defined as below.
   KK_{i,i+1}=(1-p) \rho+p Z_iZ_j \rho Z_iZ_j 
   kraus1=sqrt(1-p) I 
   kraus2=sqrt(p)ZiZj 
   Apply on (i,i+1): kraus1 \rho kraus1' + kraus2 \rho kraus2' 
   sum over all (i,i+1).
*)

kkrausChannelZZ[rho_, p_, bc_String : "obc"] := 
  Module[{ijloops, nSites, dm, j1, j2, i}, 
   Assert[StringMatchQ[bc, {"obc", "pbc"}]];
   nSites = Log2[Length[rho]];
   Assert[IntegerQ[nSites]];
   ijloops = 
    Join[Table[j1 = i; j2 = i + 1; {j1, j2}, {i, 1, nSites - 1, 2}], 
     Table[j1 = i; j2 = i + 1; {j1, j2}, {i, 2, nSites - 1, 2}]];
   If[bc == "pbc", ijloops = Join[ijloops, {{nSites, 1}}];];
   (* be cautious with below code.One should make sure dm is changed 
   with Table expanding,and equivalent to the For loop. *)
   dm = rho;
   Table[{j1, j2} = j1j2;
    dm = (1 - p)*dm + 
      p*zizjOp[nSites, j1, j2] . dm . zizjOp[nSites, j1, j2];
    dm, {j1j2, ijloops}];
   (*dm=rho;
   For[i=1,i<=Length[ijloops],i++,{j1,j2}=ijloops[[i]];
   dm=(1-p)*dm+p*zizjOp[nSites,j1,j2].dm.zizjOp[nSites,j1,j2];];*)
   dm];


(* all-to-all ZZ two-qubit quantum channel. *)
(* apply all-to-all ZZ quantum channel to a density matrix. *)
zzChannelA2Adm[nSites_Integer, p_, rho_] := Module[
   {p1, newrho, zz},
   Assert[Length[rho] == 2^nSites];
   
   p1 = p/nSites // N;
   newrho = rho;
   Table[
    zz = zizjOp[nSites, j1, j2];
    newrho = (1 - p1)*newrho + p1*zz . newrho . zz;
    , {j1, 1, nSites - 1}, {j2, j1 + 1, nSites}];
   newrho
   ];

(* probabilistically (classically) apply all-to-all ZZ gate to a \
state. *)
zzChannelA2Astate[nSites_Integer, p_, psi_] := Module[
   {p1, ijs, opz, i, j, ops, zz},
   Assert[Length[psi] == 2^nSites];
   
   p1 = p/nSites // N;
   ijs = Flatten[Table[
      {j1, j2}
      , {j1, 1, nSites - 1}, {j2, j1 + 1, nSites}], 1];
   
   opz = Association[Table[i -> 0, {i, nSites}]];
   Table[
    {i, j} = ij;
    r = RandomReal[];
    If[r <= p1, opz[i] += 1; opz[j] += 1;]
    , {ij, ijs}];
   opz = Association[Table[i -> Mod[opz[i], 2], {i, nSites}]];
   ops = Values[opz /. {0 -> "I", 1 -> "Z"}];
   zz = anyOp[nSites, Sequence @@ ops];
   
   zz . psi
   ];


End[];

EndPackage[];
