(* run on my own computer using notebook frontend. *)

Clear["Global`*"];
Remove["Global`*"];
SetDirectory[NotebookDirectory[]];
On[Assert];

<< "SpinOneHalfChain.wl";
<< "OpenSystem.wl";
<< "QIsingED.wl";
<< "SimRandomized.wl";


(* all-to-all zz channel. *)
(* analytical solution. *)
r2Analytical[nSites_Integer, p_] := Module[
   {p1, mu, r2},
   p1 = p/nSites // N;
   mu = ArcTanh[p1/(1 - p1)];

   r2 = (1./
        nSites/(nSites - 
         1))*(-nSites + (Sum[
          Binomial[nSites, m]*(nSites - 2*m)^2*
           Exp[mu*(nSites - 2 m)^2], {m, 0, nSites}])/(Sum[
          Binomial[nSites, m]*Exp[mu*(nSites - 2 m)^2], {m, 0, 
           nSites}]));

   r2
   ];

(* ED computing. *)
g = 10;
nns = {4, 5, 6, 7, 8, 9, 10};
pps = Join[{0}, Range[10]/10.];

res = Table[
   {e0, v0} = qIsingED[nSites, g, "obc"];
   rho0 = KroneckerProduct[v0, v0];
   
   {j1, j2} = {1, nSites};
   
   Table[
    rho1 = zzChannelA2Adm[nSites, p, rho0];
    nm = Tr[rho1 . rho1];
    edval = 
     Tr[rho1 . zizjOp[nSites, j1, j2] . rho1 . 
        zizjOp[nSites, j1, j2]]/Tr[rho1 . rho1];
    (*Print["p: "<>ToString[p]<>"  theoretical value: "<>ToString[
    theoval]]*)
    (*Print[StringTemplate["N: `` p: `` theoretical value: ``"][
    nSites,p,theoval]];*)
    alval = r2Analytical[nSites, p];
    {edval, alval}
    , {p, pps}], {nSites, nns}];

res // Dimensions
Transpose[res, {2, 3, 1}] // Dimensions
{edvals, alvals} = Transpose[res, {2, 3, 1}];


extNinfdata1 = 
  Table[Interpolation[Transpose[{1./nns, Transpose[edvals][[i]]}]][0]
    , {i, Length[pps]}];

extNinfdata2 = Table[
   data = Transpose[{1./nns, Transpose[edvals][[i]]}];
   ff = Fit[#, {1, x, x^2}, x] &;
   ff[data] /. {x -> 0}
   , {i, Length[pps]}];

extNinfdatasqrt = 
  Table[Interpolation[
     Transpose[{1./Sqrt[nns], Transpose[edvals][[i]]}]][0], {i, Length[pps]}];


(* plotting. *)
(* plot utilities. *)
styletext := Style[#, FontSize -> 18, Black] &;

(* Theoretical results: Analytical solution. *)
fig0 = ListLinePlot[
   Thread[{pps, #}] & /@ alvals,
   PlotRange -> {-0.01, 1},
   Frame -> True,
   FrameTicksStyle -> Directive[FontSize -> 14, Black],
   FrameLabel -> (styletext /@ {"p", "R2"}),
   (*PlotLabel->StringTemplate["g ``"][g];*)
   PlotLegends -> 
    Placed[StringTemplate["al N=``"] /@ nns, {0.1, 0.75}],
   AspectRatio -> .8
   (*ImageSize->400*)
   ];

(* ED results. *)
fig1 = ListPlot[
   Thread[{pps, #}] & /@ edvals,
   PlotRange -> {-0.01, 1},
   Frame -> True,
   FrameTicksStyle -> Directive[FontSize -> 14, Black],
   FrameLabel -> (styletext /@ {"p", "R2"}),
   PlotLegends -> 
    Placed[StringTemplate["ed N=``"] /@ nns, {0.85, 0.4}],
   AspectRatio -> .8
   ];

(* ED finite-size scaling. *)
fig2 = ListLinePlot[
   Thread[{pps, #}] & /@ {extNinfdata1, extNinfdata2, extNinfdatasqrt},
   Frame -> True,
   FrameTicksStyle -> Directive[FontSize -> 14, Black],
   FrameLabel -> (styletext /@ {"p", "R2"}),
   PlotStyle -> {{Red, Dashed}, {Orange, Dashed}, {Blue, Dashed}},
   ImageSize -> 400,
   PlotLegends -> 
    Placed[StringTemplate["inf N: ``"] /@ {"linear", "quadratic", 
       "sqrt"}, {0.15, 0.4}]
   ];

(* Theoretical results: mean field calculation. *)
mfres = Table[{u, 
    x^2 /. FindMaximum[
       Log[Cosh[2  u  x]] - u  x^2, {x, 0.5}][[2]]}, {u, 0.001, 1, 
    0.05}];
fig3 = ListLinePlot[mfres,
   PlotStyle -> Directive[Gray, Dashed],
   PlotLegends -> Placed[{"inf N: mean field"}, {0.15, 0.25}]
   ];

fig = Show[fig1, fig2, fig3, fig0, ImageSize -> 600, 
  PlotLabel -> StringTemplate["g=``"][g], 
  LabelStyle -> Directive[Bold, FontSize -> 20]]
Export[StringTemplate["figs/r2-p-g``-ed-al-mf.png"][g], fig]


(* show finite size scaling r2 vs. 1/N for each p. *)
pps
ListLinePlot[#] & /@ (Thread[{1./nns, #}] & /@ Transpose[edvals])

