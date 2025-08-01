Clear["Global`*"];
Remove["Global`*"];
SetDirectory[NotebookDirectory[]];
On[Assert];

<< "SpinOneHalfChain.wl";
<< "QIsingED.wl";
<< "OpenSystem.wl";

(* define transverse field Ising model Hamiltonian. obsolete. *)
hzzHamilt[nSites_Integer, bc_String : "obc"] := 
  Module[{h, hb},
   Assert[StringMatchQ[bc, {"obc", "pbc"}]];
   h = Sum[
     TensorProduct[IdentityMatrix[2^(j - 1)], sZ, sZ, 
         IdentityMatrix[2^(nSites - 1 - j)]] // ArrayFlatten // 
       ArrayFlatten // ArrayFlatten, {j, 1, nSites - 1}];
   If[bc == "pbc", 
    hb = TensorProduct[sZ, IdentityMatrix[2^(nSites - 2)], sZ] // 
       ArrayFlatten // ArrayFlatten;
    h += hb];
   h *= -1.0;
   h];
hxHamilt[nSites_Integer, g_] := -1.0*g*
   Sum[TensorProduct[IdentityMatrix[2^(j - 1)], sX, 
       IdentityMatrix[2^(nSites - j)]] // ArrayFlatten // 
     ArrayFlatten, {j, 1, nSites}];

(* test: check Hamiltonian correctness. *)
htisingHamilt[6, 1.5] == hzzHamilt[6] + hxHamilt[6, 1.5]
htisingHamilt[6, 1.5, "obc"] == 
 hzzHamilt[6, "obc"] + hxHamilt[6, 1.5]
htisingHamilt[6, 1.5, "pbc"] == hzzHamilt[6, "pbc"] + hxHamilt[6, 1.5]

(* test: anyOp with zjOp and zizjOp. *)
zjOp[4, 2] == anyOp[4, "I", "Z", "I", "I"]
zizjOp[4, 2, 3] == anyOp[4, "I", "Z", "Z", "I"]

nSites = 4;
hx = 1.5;
{e0, v0} = qIsingED[nSites, hx, "pbc"]

(* compute correlation functions. *)
Table[v0 . zizjOp[nSites, i, j] . v0, {i, 1, nSites}, {j, 1, 
   nSites}] // MatrixForm

(* compute order parameter. *)
Table[v0 . zjOp[nSites, j] . v0, {j, 1, nSites}] // MatrixForm

(*rho=KroneckerProduct[vec[[1]],vec[[1]]];*)
rho = KroneckerProduct[v0, v0];
Tr[rho]

rhoa = reduceDensmat[rho, 4];
Tr[rhoa]

srenyi2Entropy[rhoa]
svonnEntropy[rhoa]

Table[rhoa = reduceDensmat[rho, j];
  {svonnEntropy[rhoa], srenyi2Entropy[rhoa]}, {j, 0, 
   nSites}] // MatrixForm

ListLinePlot[Table[rhoa = reduceDensmat[rho, j];
   {srenyi2Entropy[rhoa], svonnEntropy[rhoa]}, {j, 0, nSites}] // 
  Transpose, PlotLabel -> "g=0.5"]

rho1 = kkrausChannelZZ[rho, 0.1, "pbc"];
rho1 // MatrixForm

(* check all expectation values after kraus channel. *)
rho1 = kkrausChannelZZ[rho, 0.1, "pbc"];
opset = {"I", "X", "Y", "Z"};
(*opstr={opset};*)
(* TODO: automat below. *)
opstr = Table[
  {j1, j2, j3, j4}, {j1, opset}, {j2, opset}, {j3, opset}, {j4, opset}];
opstr = Flatten[opstr, nSites - 1];
opss = anyOp[nSites, Sequence @@ #] & /@ opstr;
(*opss//Dimensions*)

Tr[rho1.#] & /@ opss // Chop
Tr[rho1]
Tr[rho.#] & /@ opss // Chop

ox = anyOp[nSites, Sequence @@ Table["X", nSites]];
rho2 = ox.rho1;
rho3 = rho2.ox; 
{Tr[rho1.rho1], Tr[rho2.rho2], Tr[rho3.rho3], Tr[rho1.rho2], Tr[rho1.rho3]}


(* test reduce density matrix partial trace. *)
mat1 = TensorProduct[{{a, b}, {c, d}}, {{x, y}, {z, w}}] // ArrayFlatten;
mat1 // MatrixForm
reduce1dm[mat1] // MatrixForm

mat2 = TensorProduct[{{a, b}, {c, d}}, {{x, y}, {z, w}}, {{k, l}, {m, n}}] //
    ArrayFlatten // ArrayFlatten;
mat2 // MatrixForm
reduce1dm[mat2] // MatrixForm
Nest[reduce1dm, mat2, 2] // MatrixForm
