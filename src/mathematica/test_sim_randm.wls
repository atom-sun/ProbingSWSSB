Clear["Global`*"];
Remove["Global`*"];
SetDirectory[NotebookDirectory[]];
On[Assert];

<< "SpinOneHalfChain.wl";
<< "OpenSystem.wl";
<< "QIsingED.wl";
<< "SimRandomized.wl";


(* test pure state. *)
(* N=4, g=1.5 ground state. ij=2,3. theoretical value: 0.120696. *)
nSites = 4;
g = 1.5;
{e0, v0} = qIsingED[nSites, g, "obc"];
v0z = zizjOp[nSites, 2, 3] . v0;
nm = Conjugate[v0] . v0;

nSamplesU = 20;
mSamplesS = 500;
ntot = nSamplesU*mSamplesS;
{t, xs} = 
  measureStatesX[v0, v0z, nSamplesU, mSamplesS] // AbsoluteTiming;
xs1 = xs/nm // Flatten;
{t, Total[xs1]/ntot, Variance[xs1]/ntot}


(* test density matrix. *)
(* N=4, g=1.5 ground state. ij=2,3. theoretical value: 0.120696. *)
nSites = 4;
g = 1.5;
{e0, v0} = qIsingED[nSites, g, "obc"];
rho0 = KroneckerProduct[v0, v0];
nm = Tr[rho0];
v0z = zizjOp[nSites, 2, 3] . v0;
rho0z = KroneckerProduct[v0z, v0z];

nSamplesU = 20;
mSamplesS = 500;
ntot = nSamplesU*mSamplesS;
{t, xs} = 
  measureStatesX[rho0, rho0z, nSamplesU, mSamplesS] // AbsoluteTiming;
xs1 = xs/nm // Flatten;
{t, Total[xs1]/ntot, Variance[xs1]/ntot}


(* demo: run for ground state. *)
nSites = 4;
g = 1.5;
{e0, v0} = qIsingED[nSites, g, "obc"];
{e0, v0}

LaunchKernels[4]
(* Science 364, 260-263 (2019). Experiment params: N=500, M=150. *)
{t, xs} = Parallelize[(<< "SimRandomized.wl";
   measureStatesX[v0, v0, #, 150]) & /@ Table[10, 50]] // AbsoluteTiming


(* all-to-all zz channel. theoretical value: 0.4780983362502319` *)
nSites = 4;
g = 1.5;
{e0, v0} = qIsingED[nSites, g, "obc"];
rho0 = KroneckerProduct[v0, v0];

p = 0.4;
rho1 = zzChannelA2Adm[nSites, p, rho0];
nm = Tr[rho1 . rho1];

{j1, j2} = {2, 3};
rho1z = zizjOp[nSites, j1, j2] . rho1 . zizjOp[nSites, j1, j2];

(* theoretical value: 0.4780983362502319` *)
theoval = 
  Tr[rho1 . zizjOp[nSites, j1, j2] . rho1 . zizjOp[nSites, j1, j2]]/
   Tr[rho1 . rho1];
Print["theoretical value: ", theoval]

nSamplesU = 1000;
mSamplesS = 1000;
{t, xs} = 
  Parallelize[(<< "SimRandomized.wl";
    measureStatesX[rho1, rho1z, #, mSamplesS]) & /@ 
     Table[1, nSamplesU]] // AbsoluteTiming;
xs1 = (Total /@ Flatten[xs, 1])/mSamplesS/nm;
{t, Total[xs1]/nSamplesU, Variance[xs1]/nSamplesU // Sqrt}

