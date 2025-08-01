Clear["Global`*"];
Remove["Global`*"];
SetDirectory[NotebookDirectory[]];
On[Assert];

<< "SpinOneHalfChain.wl";
<< "QIsingED.wl";
<< "SimSwap.wl";


(* solve coefficients for swap operators. *)
sWap = {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}};
pig1 = i0*TensorProduct[sI, sI] + a*TensorProduct[sX, sX] + 
    b*TensorProduct[sY, sY] + c*TensorProduct[sZ, sZ] // 
   ArrayFlatten;
pig1 // MatrixForm;
Solve[pig1 == sWap, {a, b, c, i0}]

(* define single-site swap operators on two copies. *)
swapOp[nSites_Integer, jSite_Integer] :=
  Module[{imats, xmats, ymats, zmats, op},
   
   Assert[jSite <= nSites];
   
   imats = Table["I", nSites];
   xmats = ReplaceAt[imats, "I" -> "X", jSite];
   ymats = ReplaceAt[imats, "I" -> "Y", jSite];
   zmats = ReplaceAt[imats, "I" -> "Z", jSite];
   
   op = anyOp[2*nSites, Sequence @@ imats, Sequence @@ imats] +
     anyOp[2*nSites, Sequence @@ xmats, Sequence @@ xmats] +
     anyOp[2*nSites, Sequence @@ ymats, Sequence @@ ymats] +
     anyOp[2*nSites, Sequence @@ zmats, Sequence @@ zmats];
   
   (1/2)*op
   ];

(* test: N=4, g=1.5 ground state. ij=2,3. theoretical value: 0.120696. *)
nSites = 4;
g = 1.5;
{e0, v0} = qIsingED[nSites, g, "obc"];
{e0, v0}

(* whole parity: swap (on two copies) of total system. *)
wholeParity = Dot[Sequence @@ Table[swapOp[nSites, i], {i, nSites}]];
wholeParity // Dimensions

rho0 = KroneckerProduct[v0, v0];
rho0z = zizjOp[nSites, 2, 3] . rho0 . zizjOp[nSites, 3, 2];
Tr[wholeParity . (TensorProduct[rho0, rho0z] // ArrayFlatten)]


(* test Swap projection. *)
nSites = 4;
iSite = 2;
jSite = 3;
{uu, dd, ss, aa} = projSwapOp[nSites, iSite, jSite];
generalOp[nSites, iSite, (sI + sZ)/2, jSite, (sI + sZ)/2] == uu
generalOp[nSites, iSite, (sI - sZ)/2, jSite, (sI - sZ)/2] == dd
generalOp[nSites, iSite, "I", jSite, "I"]*3/4 + 
  generalOp[nSites, iSite, "X", jSite, "X"]/4 + 
  generalOp[nSites, iSite, "Y", jSite, "Y"]/4 + 
  generalOp[nSites, iSite, "Z", jSite, "Z"]/4 - uu - dd == ss
generalOp[nSites, iSite, "I", jSite, "I"] - 
  generalOp[nSites, iSite, "X", jSite, "X"] - 
  generalOp[nSites, iSite, "Y", jSite, "Y"] - 
  generalOp[nSites, iSite, "Z", jSite, "Z"] == aa*4


(* example. vector of state. *)
nSites = 4;
g = 1.5;
{e0, v0} = qIsingED[nSites, g, "obc"];
v0z = zizjOp[nSites, 2, 3] . v0;

(* for pure state, purity should be 1. *)
{t, p2s} = 
  ParallelTable[
    << "SimSwap.wl";
    measureParity[v0, v0], {i, 1, 1000}] // AbsoluteTiming;
{t, Total[p2s]/1000., Variance[p2s] // N}

{t, p2s} = 
  ParallelTable[
    << "SimSwap.wl";
    measureParity[v0, v0z], {i, 1, 1000}] // AbsoluteTiming;
{t, Total[p2s]/1000., Variance[p2s] // N}

(* example. density matrix. *)
nSites = 4;
g = 1.5;
{e0, v0} = qIsingED[nSites, g, "obc"];
v0z = zizjOp[nSites, 2, 3] . v0;
rho0 = KroneckerProduct[v0, v0];
rho0z = KroneckerProduct[v0z, v0z];

(* for pure state, purity should be 1.*)
{t, p2s} = 
  ParallelTable[
    << "SimSwap.wl";
    measureParity[rho0, rho0], {i, 1, 1000}] // AbsoluteTiming;
{t, Total[p2s]/1000., Variance[p2s] // N}

{t, p2s} = 
  ParallelTable[
    << "SimSwap.wl";
    measureParity[rho0, rho0z], {i, 1, 1000}] // AbsoluteTiming;
{t, Total[p2s]/1000., Variance[p2s] // N}
