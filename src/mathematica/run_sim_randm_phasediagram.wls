(*run on cluster using wolfram script.*)
(*math-script script.wls*)
On[Assert];

path0=DirectoryName[$InputFileName];
(* path0 = NotebookDirectory[]; *)
Print[path0];

path1 = Directory[];
Print[path1];
dumpPath = FileNameJoin[{path1, "dump"}];
If[! DirectoryQ[dumpPath], 
  CreateDirectory[dumpPath, CreateIntermediateDirectories -> True]];

SetDirectory[path0];
Print[Directory[]];

{pkgspinchain, pkgqising, pkgopensystem, pkgrandm} = 
  pkgs = FileNameJoin[{path0, #}] & /@ {"SpinOneHalfChain.wl", 
     "QIsingED.wl", "OpenSystem.wl", "SimRandomized.wl"};

Get /@ pkgs

(*check pkg import.*)
{e0, v0} = qIsingED[4, 1.5, "obc"];
Print[e0];
Print[v0 // Dimensions];
Print[v0];


pps = Range[0., 0.75, 0.075]
ggs = Tan[Range[Pi/4., Pi/2., Pi/41.]]

nSites = 7;
nSamples = 10000;


(* run phase diagram. *)
{t, xsall} = Table[
  {e0, v0} = qIsingED[nSites, g, "obc"];
  {j1, j2} = {1, nSites};
  xs0 = ParallelTable[
    Get[pkgopensystem];
    Get[pkgrandm];
    v11 = zzChannelA2Astate[nSites, p, v0];
    v12 = zzChannelA2Astate[nSites, p, v0];
    x0 = measureStatesX[v11, v12, 1, 1];
    x0, {j, nSamples}];
  DumpSave[
   FileNameJoin[{dumpPath, 
     StringTemplate["xs0-g``-p``-nSites``-nSamples``.mx"][g, p, 
      nSites, nSamples]}], xs0];
  xs1 = ParallelTable[
    Get[pkgopensystem];
    Get[pkgrandm];
    v1 = zzChannelA2Astate[nSites, p, v0];
    v1z = zizjOp[nSites, j1, j2] . zzChannelA2Astate[nSites, p, v0];
    x1 = measureStatesX[v1, v1z, 1, 1];
    x1, {j, nSamples}];
  DumpSave[
   FileNameJoin[{dumpPath, 
     StringTemplate["xs1-g``-p``-nSites``-nSamples``.mx"][g, p, 
      nSites, nSamples]}], xs1];
  {xs0, xs1}
  , {g, ggs}, {p, pps}] // AbsoluteTiming;


DumpSave[
  FileNameJoin[{dumpPath, 
    StringTemplate["xsall-nSites``-nSamples``.mx"][nSites, 
     nSamples]}], xsall];
DumpSave[
  FileNameJoin[{dumpPath, 
    StringTemplate["xsall-nSites``-nSamples``-params.mx"][nSites, 
     nSamples]}], {ggs, pps}];

Print[t];
Print[xsall // Dimensions];

