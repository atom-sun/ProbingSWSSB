(* run on cluster using wolfram script. *)
(* math -script script.wls *)

On[Assert];

path0 = DirectoryName[$InputFileName];
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

(* check pkg import. *)
{e0, v0} = qIsingED[4, 1.5, "obc"];
Print[e0];
Print[v0 // Dimensions];
Print[v0];


(*all-to-all zz channel.probabilistically.*)
g = 10;
p = 0.3;
nSamples = 1000;
nSitesChoices = {4, 5, 6};

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
  , {nSites, nSitesChoices}] // AbsoluteTiming;


DumpSave[
  FileNameJoin[{dumpPath, 
    StringTemplate["xsall-g``-p``-nSamples``.mx"][g, p, nSamples]}], 
  xsall];
DumpSave[
  FileNameJoin[{dumpPath, 
    StringTemplate["xsall-g``-p``-nSamples``-params.mx"][g, p, 
     nSamples]}], {g, p, nSamples, nSitesChoices}];

Print[t];
Print[xsall // Dimensions];

