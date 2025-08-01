(* this script is obsolete. *)


(* Mathematica Source File *)
(* Transverse field Ising model *)
(* :Author: Ning Sun *)
(* :Date: 2024-10-30 *)
(* :Mathematica Version: 14.0 *)
(* :Copyright: (c) 2024 atom-sun *)

Clear["Global`*"];
Remove["Global`*"];
SetDirectory[NotebookDirectory[]];

(* constants. *)
sI = IdentityMatrix[2];
sX = PauliMatrix[1];
sY = PauliMatrix[2];
sZ = PauliMatrix[3];

(* define any operators on the spin-1/2 chain. *)
anyOp[nSites_Integer, ops__String] := Module[
   {mapOps, op},
   
   On[Assert];
   Assert[Length[{ops}] == nSites];
   mapOps = <|"I" -> sI, "X" -> sX, "Y" -> sY, "Z" -> sZ|>;
   Assert[AllTrue[{ops}, MemberQ[Keys[mapOps], #] &]];
   
   op = Nest[ArrayFlatten, 
     TensorProduct @@ Table[mapOps[j], {j, {ops}}], nSites - 1];
   op
   ];

(* define transverse-field Ising model Hamiltonian. *)
htisingHamilt[nSites_Integer, g_, bc_String : "obc"] := Module[
   {hzz, hx, hb},
   hzz = 
    Sum[anyOp[nSites, 
      Sequence @@ 
       Join[Table["I", i - 1], Table["Z", 2], 
        Table["I", nSites - 1 - i]]], {i, 1, nSites - 1}];
   If[bc == "pbc",
    hb = 
     TensorProduct[sZ, IdentityMatrix[2^(nSites - 2)], sZ] // 
       ArrayFlatten // ArrayFlatten;
    hzz += hb];
   hx = Sum[
     anyOp[nSites, 
      Sequence @@ 
       Join[Table["I", i - 1], {"X"}, Table["I", nSites - i]]], {i, 1,
       nSites}];
   -1.0*(hzz + g*hx)
   ];

(* define transverse field Ising model Hamiltonian obsolete. *)
hzzHamilt[nSites_Integer, bc_String : "obc"] :=
  Module[{h, hb},
   On[Assert];
   Assert[StringMatchQ[bc, {"obc", "pbc"}]];
   h = Sum[
     TensorProduct[IdentityMatrix[2^(j - 1)], sZ, sZ, 
         IdentityMatrix[2^(nSites - 1 - j)]] // ArrayFlatten // 
       ArrayFlatten // ArrayFlatten, {j, 1, nSites - 1}];
   If[bc == "pbc",
    hb = 
     TensorProduct[sZ, IdentityMatrix[2^(nSites - 2)], sZ] // 
       ArrayFlatten // ArrayFlatten;
    h += hb];
   h *= -1.0;
   h
   ];
hxHamilt[nSites_Integer, g_] := -1.0*g*
   Sum[TensorProduct[IdentityMatrix[2^(j - 1)], sX, 
       IdentityMatrix[2^(nSites - j)]] // ArrayFlatten // 
     ArrayFlatten, {j, 1, nSites}];

(* test: check Hamiltonian correctness. *)
htisingHamilt[6, 1.5] == hzzHamilt[6] + hxHamilt[6, 1.5]
htisingHamilt[6, 1.5, "obc"] == 
 hzzHamilt[6, "obc"] + hxHamilt[6, 1.5]
htisingHamilt[6, 1.5, "pbc"] == hzzHamilt[6, "pbc"] + hxHamilt[6, 1.5]

(* define correlation functions.  <ZiZj> *)
zizjOp[nSites_Integer, i_Integer, j_Integer] /; i < j := 
  TensorProduct[IdentityMatrix[2^(i - 1)], sZ, 
       IdentityMatrix[2^(j - i - 1)], sZ, 
       IdentityMatrix[2^(nSites - j)]] // ArrayFlatten // 
     ArrayFlatten // ArrayFlatten // ArrayFlatten;
zizjOp[nSites_Integer, i_Integer, j_Integer] /; i > j := 
  zizjOp[nSites, j, i];
zizjOp[nSites_Integer, i_Integer, j_Integer] /; i == j := 
  TensorProduct[IdentityMatrix[2^(j - 1)], sZ . sZ, 
     IdentityMatrix[2^(nSites - j)]] // ArrayFlatten // ArrayFlatten;

(* define order parameter. <Zj> *)
zjOp[nSites_Integer, j_Integer] := 
  TensorProduct[IdentityMatrix[2^(j - 1)], sZ, 
     IdentityMatrix[2^(nSites - j)]] // ArrayFlatten // ArrayFlatten;

(* test: anyOp with zjOp and zizjOp. *)
zjOp[4, 2] == anyOp[4, "I", "Z", "I", "I"]
zizjOp[4, 2, 3] == anyOp[4, "I", "Z", "Z", "I"]

(* define entanglement entropies. *)
(*svonnEntropy[rho_]:=-Tr[rho.MatrixLog[rho]];*)
(*svonnEntropy[rho_]:=-Tr[rho.Log[rho]];*)
(*svonnEntropy[rho_]:=-Tr[rho.MatrixLog[rho+10^(-4)*Identity[Length[\
rho]]]];*)
svonnEntropy[rho_] :=
  Module[{r},
   r = Eigenvalues[rho];
   Chop[-r*Log[r]] // Total
   ];
srenyi2Entropy[rho_] := -Log[Tr[rho . rho]];

(* define reduced density matrix *)
(*ResourceFunction["MatrixPartialTrace"]*)
reduce1dm[rho_] :=
  Module[{l, dm1, dm2},
   l = Length[rho];
   dm1 = rho[[1 ;; Length[rho] ;; 2, 1 ;; Length[rho] ;; 2]];
   dm2 = rho[[2 ;; Length[rho] ;; 2, 2 ;; Length[rho] ;; 2]];
   dm1 + dm2
   ];
reduceDensmat[rho_, j_Integer] := Nest[reduce1dm, rho, j];
densMat[psi_] := KroneckerProduct[psi, psi];

nSites = 4;
hx = 1.5;

(* obc *)
ham = hzzHamilt[nSites] + hxHamilt[nSites, hx];
ham = htisingHamilt[nSites, hx];

(* pbc *)
ham = hzzHamilt[nSites, "pbc"] + hxHamilt[nSites, hx];
ham = htisingHamilt[nSites, hx, "pbc"];

{erg, vec} = Eigensystem[N[ham]];
where = DeleteDuplicates[
   Flatten[Map[Position[erg, #] &, Sort[erg, Less]]]];
erg = erg[[where]];
vec = vec[[where]];

(* ground state energy. *)
erg[[1]]

(* ground state: *)
vec[[1]];

(* compute correlation functions. *)
Table[vec[[1]] . zizjOp[nSites, i, j] . vec[[1]], {i, 1, nSites}, {j, 
   1, nSites}] // MatrixForm

(* compute order parameter. *)
Table[vec[[1]] . zjOp[nSites, j] . vec[[1]], {j, 1, 
   nSites}] // MatrixForm

(*rho=KroneckerProduct[vec[[1]],vec[[1]]];*)
rho = densMat[vec[[1]]];
Tr[rho]

rhoa = reduceDensmat[rho, 4];
Tr[rhoa]

srenyi2Entropy[rhoa]

svonnEntropy[rhoa]

Table[rhoa = reduceDensmat[rho, j];
  {svonnEntropy[rhoa], srenyi2Entropy[rhoa]},
  {j, 0, nSites}] // MatrixForm

ListLinePlot[Table[rhoa = reduceDensmat[rho, j];
   {srenyi2Entropy[rhoa], svonnEntropy[rhoa]},
   {j, 0, nSites}] // Transpose, PlotLabel -> "g=0.5"]

(* kkraus zz quantum channel. *)
(*
KKrause quantum channel ZZ is defined as below.

KK_{i,i+1}=(1-p) \rho+p Z_iZ_j \rho Z_iZ_j
kraus1=sqrt(1-p) I
kraus2=sqrt(p)ZiZj
Apply on (i,i+1):kraus1 \rho kraus1'+kraus2 \rho kraus2'
sum over all (i,i+1).
*)

kkrausChannelZZ[rho_, p_, bc_String : "obc"] := Module[
   {ijloops, nSites, dm, j1, j2, i},
   
   On[Assert];
   Assert[StringMatchQ[bc, {"obc", "pbc"}]];
   nSites = Log2[Length[rho]];
   Assert[IntegerQ[nSites]];
   
   ijloops = Join[
     Table[j1 = i; j2 = i + 1; {j1, j2}, {i, 1, nSites - 1, 2}],
     Table[j1 = i; j2 = i + 1; {j1, j2}, {i, 2, nSites - 1, 2}]
     ];
   If[bc == "pbc", ijloops = Join[ijloops, {{nSites, 1}}];];
   
   (*be cautious with below code. 
   One should make sure dm is changed with Table expanding,
   and equivalent to the For loop.
   *)
   dm = rho;
   Table[{j1, j2} = j1j2;
    dm = (1 - p)*dm + 
      p*zizjOp[nSites, j1, j2] . dm . zizjOp[nSites, j1, j2];
    dm, {j1j2, ijloops}];
   (*dm=rho;
   For[
   i=1,
   i<=Length[ijloops],
   i++,
   {j1,j2}=ijloops[[i]];
   dm=(1-p)*dm+p*zizjOp[nSites,j1,j2].dm.zizjOp[nSites,j1,j2];
   ];*)
   dm
   ];

rho1 = kkrausChannelZZ[rho, 0.1, "pbc"];
rho1 // MatrixForm

(* check all expectation values after kraus channel. *)
rho1 = kkrausChannelZZ[rho, 0.1, "pbc"];
opset = {"I", "X", "Y", "Z"};
(*opstr={opset};*)
(* TODO: automat below. *)
opstr = Table[{j1, j2, j3, j4}, {j1, opset}, {j2, opset}, {j3, 
    opset}, {j4, opset}];
opstr = Flatten[opstr, nSites - 1];
opss = anyOp[nSites, Sequence @@ #] & /@ opstr;
(*opss//Dimensions*)

Tr[rho1 . #] & /@ opss // Chop

Tr[rho1]

Tr[rho . #] & /@ opss // Chop

ox = anyOp[nSites, Sequence @@ Table["X", nSites]];
rho2 = ox . rho1;
rho3 = rho2 . ox;

{Tr[rho1 . rho1], Tr[rho2 . rho2], Tr[rho3 . rho3], Tr[rho1 . rho2], 
 Tr[rho1 . rho3]}
