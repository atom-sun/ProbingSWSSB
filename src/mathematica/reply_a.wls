SetDirectory[NotebookDirectory[]];

<< "SpinOneHalfChain.wl";
<< "OpenSystem.wl";
<< "QIsingED.wl";


Plot[Sin[x],{x,-10,10}]


colors = System`PlotThemeDump`$ThemeDefaultGradient
color1 = colors[[-1]];
color2 = colors[[2]];
{color1, color2}


frame1[legend_] := Framed[legend,
   Background -> White,
   FrameStyle -> Black,
   FrameMargins -> 0,
   RoundingRadius -> 10
   ];
{circle, uptriangle, diamond, square, downtriangle} =
  Charting`CommonDump`GraphicsOpenPlotMarkers[][[;; 5]];
leg1 = SwatchLegend[
   {Style["\!\(\*SubscriptBox[\(P\), \(I\)]\)", color1], 
    Style["\!\(\*SubscriptBox[\(P\), \(ZZ\)]\)", color2]},
   LabelStyle -> {FontSize -> 18},
   LegendFunction -> frame1,
   LegendMarkers -> {circle, uptriangle}];

aspratio = 4/5;
imgsize = 400;


ArcTan[10.]/Pi

ArcTan[2.]/Pi


(* Q1.1  eps vs. Na *)
nSites = 7;
nSamplesPw = {3, 4, 5, 6, 7};
nSamplesLst = 10^nSamplesPw;

Table[
 
 {e0, v0} = qIsingED[nSites, g, "obc"];
 rho0 = KroneckerProduct[v0, v0];
 {j1, j2} = {1, nSites};
 rho1 = zzChannelA2Adm[nSites, p, rho0];
 nm = Tr[rho1 . rho1];
 x0 = Tr[rho1 . rho1];
 x1 = Tr[
   rho1 . zizjOp[nSites, j1, j2] . rho1 . zizjOp[nSites, j1, j2]];
 resED = {x0, x1};
 Print[resED];
 
 Import[StringTemplate[
    "data-reply/xs0-g``-p``-nSites``-nSamples10000000.mx"][g, p, 
   nSites, nSamples]];
 Import[StringTemplate[
    "data-reply/xs1-g``-p``-nSites``-nSamples10000000.mx"][g, p, 
   nSites, nSamples]];
 xs0 = xs0 // Flatten;
 xs1 = xs1 // Flatten;
 resRM = Table[
   xx0 = RandomChoice[xs0, nSamples] - x0;
   xx1 = RandomChoice[xs1, nSamples] - x1;
   v0 = StandardDeviation[xx0]/Sqrt[nSamples];
   v1 = StandardDeviation[xx1]/Sqrt[nSamples];
   {v0, v1}
   , {nSamples, nSamplesLst}];
 Print[resRM];
 
 data = Thread[{nSamplesLst, #}] & /@ (resRM // Transpose);
 model1 = LinearModelFit[{Log[#1], Log[#2]} & @@@ N[data[[1]]], {x}, x];
 model2 = LinearModelFit[{Log[#1], Log[#2]} & @@@ N[data[[2]]], {x}, x];
 fig1 = ListLinePlot[
   data,
   AspectRatio -> aspratio,
   Axes -> False,
   Frame -> True,
   FrameLabel -> (Style[#, 
        FontSize -> 24] & /@ {"\!\(\*SubscriptBox[\(N\), \(a\)]\)", 
       "\[Epsilon]"}),
   FrameStyle -> Thickness[0.004],
   ImageSize -> imgsize,
   LabelStyle -> {18, Black, FontFamily -> "Times"},
   PlotMarkers -> "OpenMarkers",
   PlotStyle -> {{Thickness[0.008], color1}, {Thickness[0.008], 
      color2}},
   PlotLegends -> Placed[leg1, {0.8, 0.8}],
   PlotLabel -> StringTemplate["g=``  \[Mu]=``   ``sites"][g, p, nSites],
   
   ScalingFunctions -> {"Log", "Log"}
   ];
 fig2 = Plot[(model1[Log[x]] + model2[Log[x]])/2, {x, 10^2, 10^8}, 
   PlotStyle -> {Dashed, Gray, Thickness[0.005]}, 
   ScalingFunctions -> {"Log", Automatic}];
 fig3 = Show[fig1, fig2];
 Export[StringTemplate[
    "figs/reply/sample_complexity-g``-\[Mu]``-nSites``.pdf"][g, p, 
   nSites], fig3];
 Show[fig3]

 , {g, {10, 2}}, {p, {0.3, 0.7}}]


(* Q1.2  Na vs. nSites, fix eps=0.05 *)
eps = 0.05;
nSitesLst = {4, 5, 6, 7, 8};
nSamplesPw = {3, 4, 5, 6, 7};
nSamplesLst = 10^nSamplesPw;

Table[
 
 nc = Table[
   
   {e0, v0} = qIsingED[nSites, g, "obc"];
   rho0 = KroneckerProduct[v0, v0];
   {j1, j2} = {1, nSites};
   rho1 = zzChannelA2Adm[nSites, p, rho0];
   nm = Tr[rho1 . rho1];
   x0 = Tr[rho1 . rho1];
   x1 = Tr[
     rho1 . zizjOp[nSites, j1, j2] . rho1 . zizjOp[nSites, j1, j2]];
   resED = {x0, x1};
   
   Import[
    StringTemplate[
      "data-reply/xs0-g``-p``-nSites``-nSamples10000000.mx"][g, p, 
     nSites, nSamples]];
   Import[
    StringTemplate[
      "data-reply/xs1-g``-p``-nSites``-nSamples10000000.mx"][g, p, 
     nSites, nSamples]];
   xs0 = xs0 // Flatten;
   xs1 = xs1 // Flatten;
   
   resRM = Table[
     xx0 = RandomChoice[xs0, nSamples] - x0;
     xx1 = RandomChoice[xs1, nSamples] - x1;
     v0 = StandardDeviation[xx0]/Sqrt[nSamples];
     v1 = StandardDeviation[xx1]/Sqrt[nSamples];
     {v0, v1}
     , {nSamples, nSamplesLst}];
   
   data = Thread[{nSamplesPw, #}] & /@ (resRM // Transpose);
   n0 = FindRoot[
     Interpolation[data[[1]], InterpolationOrder -> 1][x] == eps, {x, 
      3}];
   n1 = FindRoot[
     Interpolation[data[[2]], InterpolationOrder -> 1][x] == eps, {x, 
      3}];
   {n0, n1}
   
   , {nSites, nSitesLst}];
 
 data2 = Thread[{nSitesLst, #}] & /@ Transpose[10^x /. nc];
 
 model1 = LinearModelFit[{#1, Log10[#2]} & @@@ N[data2[[1]]], {x}, x];
 model2 = LinearModelFit[{#1, Log10[#2]} & @@@ N[data2[[2]]], {x}, x];
 slp1 = model1["BestFitParameters"][[2]];
 slp2 = model2["BestFitParameters"][[2]];
 leg2 = SwatchLegend[
   {Style["\!\(\*SubscriptBox[\(P\), \(I\)]\):  " Superscript[
   ToString[NumberForm[10^slp1, {\[Infinity], 3}]], "N"], color1],
   Style["\!\(\*SubscriptBox[\(P\), \(ZZ\)]\):  " Superscript[
   ToString[NumberForm[10^slp2, {\[Infinity], 3}]], "N"], color2]},
   LabelStyle -> {FontSize -> 18},
   LegendFunction -> frame1,
   LegendMarkers -> {circle, uptriangle}];
 fig1 = ListLinePlot[
   data2,
   AspectRatio -> aspratio,
   Axes -> False,
   Frame -> True,
   FrameLabel -> (Style[#, 
        FontSize -> 
         24] & /@ {"\!\(\*SubscriptBox[\(N\), \(sites\)]\)", 
       "\!\(\*SubscriptBox[\(N\), \(a\)]\)"}),
   FrameStyle -> Thickness[0.004],
   ImageSize -> imgsize,
   LabelStyle -> {18, Black, FontFamily -> "Times"},
   PlotMarkers -> "OpenMarkers",
   PlotStyle -> {{Thickness[0.008], color1}, {Thickness[0.008], 
      color2}},
   PlotLegends -> Placed[leg2, {0.75, 0.2}],
   PlotLabel -> 
    StringTemplate["g=``  \[Mu]=``   fix eps=``"][g, p, eps],
   
   ScalingFunctions -> "Log"
   ];
 fig2 = Plot[{10^model1[x], 10^model2[x]}, {x, 3, 9},
   PlotStyle -> {{Dashed, Thickness[0.005], color1}, {Dashed, 
      Thickness[0.005], color2}},
   ScalingFunctions -> "Log"];
 fig3 = Show[fig1, fig2];
 Export[StringTemplate["figs/reply/sample_sites-g``-\[Mu]``-eps``.pdf"][g,
    p, eps], fig3];
 Show[fig3]
 
 , {g, {10, 2}}, {p, {0.3, 0.7}}]


(* Q1.3  Na vs g *)
(* Q1.3  Na vs. g, fix eps=0.05 *)
eps = 0.05;
ggs0 = Range[1/4., 1/2., 1/41.];
ggs = Tan[Pi*ggs0];
g0Lst = ggs0[[2 ;; -2]];
gLst = ggs[[2 ;; -2]];
nSamplesPw = {3, 4, 5, 6};
nSamplesLst = 10^nSamplesPw;

nSites = 7;
p = 0.3;

nc = Table[
  
  {e0, v0} = qIsingED[nSites, g, "obc"];
  rho0 = KroneckerProduct[v0, v0];
  {j1, j2} = {1, nSites};
  rho1 = zzChannelA2Adm[nSites, p, rho0];
  nm = Tr[rho1 . rho1];
  x0 = Tr[rho1 . rho1];
  x1 = Tr[
    rho1 . zizjOp[nSites, j1, j2] . rho1 . zizjOp[nSites, j1, j2]];
  resED = {x0, x1};
  
  Import[
   StringTemplate[
     "data-reply/xs0-g``-p``-nSites``-nSamples1000000.mx"][g, p, 
    nSites, nSamples]];
  Import[
   StringTemplate[
     "data-reply/xs1-g``-p``-nSites``-nSamples1000000.mx"][g, p, 
    nSites, nSamples]];
  xs0 = xs0 // Flatten;
  xs1 = xs1 // Flatten;
  
  resRM = Table[
    xx0 = RandomChoice[xs0, nSamples] - x0;
    xx1 = RandomChoice[xs1, nSamples] - x1;
    v0 = StandardDeviation[xx0]/Sqrt[nSamples];
    v1 = StandardDeviation[xx1]/Sqrt[nSamples];
    {v0, v1}
    , {nSamples, nSamplesLst}];
  
  data = Thread[{nSamplesPw, #}] & /@ (resRM // Transpose);
  n0 = FindRoot[
    Interpolation[data[[1]], InterpolationOrder -> 1][x] == eps, {x, 
     3}];
  n1 = FindRoot[
    Interpolation[data[[2]], InterpolationOrder -> 1][x] == eps, {x, 
     3}];
  {n0, n1}
  
  , {g, gLst}]

data2 = Thread[{g0Lst, #}] & /@ Transpose[10^x /. nc]

fig = ListLinePlot[
   data2,
   AspectRatio -> aspratio,
   Axes -> False,
   Frame -> True,
   FrameLabel -> (Style[#, FontSize -> 24] & /@ {"arctan(g)/\[Pi]", 
       "\!\(\*SubscriptBox[\(N\), \(a\)]\)"}),
   FrameStyle -> Thickness[0.004],
   ImageSize -> imgsize,
   LabelStyle -> {18, Black, FontFamily -> "Times"},
   PlotMarkers -> "OpenMarkers",
   PlotStyle -> {{Thickness[0.008], color1}, {Thickness[0.008], 
      color2}},
   PlotLegends -> Placed[leg1, {0.8, 0.5}],
   PlotLabel -> 
    StringTemplate["Nsites=``  \[Mu]=``   fix eps=``"][nSites, p, eps],
   
   ScalingFunctions -> "Log"
   ];
Export[StringTemplate["figs/reply/sample_g-nSites``-\[Mu]``-eps``.pdf"][
   nSites, p, eps], fig];
Show[fig]


(* Q1.4  Na vs mu *)
(* Q1.4  Na vs. mu, fix eps=0.05 *)
eps = 0.05;
pLst = pps = Range[0., 0.75, 0.075];
ggs0 = Range[1/4., 1/2., 1/41.];
ggs = Tan[Pi*ggs0];

nSamplesPw = {3, 4, 5, 6};
nSamplesLst = 10^nSamplesPw;

nSites = 7;
g = ggs[[-2]];

nc = Table[
  
  {e0, v0} = qIsingED[nSites, g, "obc"];
  rho0 = KroneckerProduct[v0, v0];
  {j1, j2} = {1, nSites};
  rho1 = zzChannelA2Adm[nSites, p, rho0];
  nm = Tr[rho1 . rho1];
  x0 = Tr[rho1 . rho1];
  x1 = Tr[
    rho1 . zizjOp[nSites, j1, j2] . rho1 . zizjOp[nSites, j1, j2]];
  resED = {x0, x1};
  
  Import[
   StringTemplate[
     "data-reply/xs0-g``-p``-nSites``-nSamples1000000.mx"][g, p, 
    nSites, nSamples]];
  Import[
   StringTemplate[
     "data-reply/xs1-g``-p``-nSites``-nSamples1000000.mx"][g, p, 
    nSites, nSamples]];
  xs0 = xs0 // Flatten;
  xs1 = xs1 // Flatten;
  
  resRM = Table[
    xx0 = RandomChoice[xs0, nSamples] - x0;
    xx1 = RandomChoice[xs1, nSamples] - x1;
    v0 = StandardDeviation[xx0]/Sqrt[nSamples];
    v1 = StandardDeviation[xx1]/Sqrt[nSamples];
    {v0, v1}
    , {nSamples, nSamplesLst}];
  
  data = Thread[{nSamplesPw, #}] & /@ (resRM // Transpose);
  n0 = FindRoot[
    Interpolation[data[[1]], InterpolationOrder -> 1][x] == eps, {x, 
     3}];
  n1 = FindRoot[
    Interpolation[data[[2]], InterpolationOrder -> 1][x] == eps, {x, 
     3}];
  {n0, n1}
  
  , {p, pLst}]

data2 = Thread[{pLst, #}] & /@ Transpose[10^x /. nc]

fig = ListLinePlot[
   data2,
   AspectRatio -> aspratio,
   Axes -> False,
   Frame -> True,
   FrameLabel -> (Style[#, FontSize -> 24] & /@ {"\[Mu]", 
       "\!\(\*SubscriptBox[\(N\), \(a\)]\)"}),
   FrameStyle -> Thickness[0.004],
   ImageSize -> imgsize,
   LabelStyle -> {18, Black, FontFamily -> "Times"},
   PlotMarkers -> "OpenMarkers",
   PlotStyle -> {{Thickness[0.008], color1}, {Thickness[0.008], 
      color2}},
   PlotLegends -> Placed[leg1, {0.8, 0.5}],
   PlotLabel -> 
    StringTemplate["Nsites=``  g=``   fix eps=``"][nSites, g, eps],
   
   ScalingFunctions -> "Log"
   ];
Export[StringTemplate["figs/reply/sample_mu-nSites``-g``-eps``.pdf"][
   nSites, g, eps], fig];
Show[fig]


(* Q5. fit exponent. |mu-muc|^beta, fit the exponent beta, for fix g. \
*)
SetDirectory[NotebookDirectory[]];
FileNames["*", Directory[]];
SetEnvironment[
  "PATH" -> "/opt/homebrew/bin:" <> Environment["PATH"]];
Environment["PATH"]

data = RunProcess[{"julia", "-e", 
    "using JLD2; data=load(\"data/gs-avgzizjs2.jld2\"); \
println(data)"}];
(*Print[data["StandardOutput"]];*)

out = data["StandardOutput"];
StringPosition[out, "\"sumzizjs2\""];
StringPosition[out, "\"params\""];
StringPosition[out, "\"sumzizjs\""];

{starts, tmp} = StringPosition[out, "["] // Transpose;
{ends, tmp} = StringPosition[out, "]"] // Transpose;

sumzizjs2 = 
  ToExpression /@ (StringTrim /@ 
     StringSplit[StringTake[out, {starts[[1]] + 1, ends[[1]] - 1}], 
      ","]);
paramgs = 
  ToExpression /@ (StringTrim /@ 
     StringSplit[StringTake[out, {starts[[2]] + 1, ends[[2]] - 1}], 
      ","]);

ddata0 = Thread[{paramgs, sumzizjs2/100}];
fun0 = Interpolation[ddata0];


(* g=10 or 2 around. *)
pps = Range[0., 0.75, 0.075];
ggs0 = Range[1/4., 1/2., 1/41.];
ggs = Tan[Pi*ggs0]

nSites = 7;
nSamples = 10000;

g = ggs[[-2]]
(*g=ggs[[5]]*)

ddata0 = Thread[{paramgs, sumzizjs2/100}];
fun0 = Interpolation[ddata0];
mu0 = 1./(2. + 2*fun0[g]);


data = Table[
   Import[
    StringTemplate["dump/xs0-g``-p``-nSites``-nSamples``.mx"][g, p, 
     nSites, nSamples]];
   Import[
    StringTemplate["dump/xs1-g``-p``-nSites``-nSamples``.mx"][g, p, 
     nSites, nSamples]];
   
   xs0 = xs0 // Flatten;
   xs1 = xs1 // Flatten;
   
   keys = xs0 // Counts // Keys // Sort;
   p0 = Table[Counts[xs0][k], {k, keys}] // N;
   p0 /= Total[p0];
   p1 = Table[Counts[xs1][k], {k, keys}] // N;
   p1 /= Total[p1];
   relentrp = p0*Log[p0/p1] // Total;
   
   {p - mu0, relentrp}
   
   , {p, pps}];

data
ListLinePlot[data]

fitdata = {Abs[#1], #2} & @@@ 
  Select[data, (#[[1]] < 0 && #[[1]] > -0.35) &]
model = NonlinearModelFit[
  fitdata, {cf*x^bf + df}, {{cf, 0.1}, {bf, 1.1}, {df, 0.1}}, x]
fitdata2 = {Abs[#1], #2} & @@@ Select[data, (#[[1]] < 0) &]
slope = (bf /. model["BestFitParameters"])


fig1 = ListLinePlot[data,
   
   AspectRatio -> aspratio,
   AxesLabel -> (Style[#, Black, FontFamily -> "Times", 
        FontSize -> 
         14] & /@ {"\[Mu]-\!\(\*SubscriptBox[\(\[Mu]\), \(c\)]\)", 
       "KL"}),
   AxesStyle -> Thickness[0.004],
   ImageSize -> 320,
   PlotMarkers -> "OpenMarkers",
   PlotStyle -> {Blue, Thickness[0.008]}
   ];

fig2 =
  ListLinePlot[fitdata2,
   Axes -> False,
   AspectRatio -> aspratio,
   Frame -> True,
   FrameLabel -> (Style[#, Black, FontFamily -> "Times", 
        FontSize -> 
         16] & /@ {"\!\(\*SubscriptBox[\(\[Mu]\), \(c\)]\)-\[Mu]", 
       "KL"}),
   FrameStyle -> Thickness[0.004],
   ImageSize -> 400,
   LabelStyle -> {16, Black, FontFamily -> "Times"},
   PlotMarkers -> "OpenMarkers",
   PlotStyle -> {Blue, Thickness[0.01]}
   ];
leg2 = SwatchLegend[
   {Style[StringTemplate["\[Beta]=``"][slope], Red]},
   LabelStyle -> {FontSize -> 18},
   LegendFunction -> frame1,
   LegendMarkers -> {"--"}];
fig3 = Plot[
   model[x], {x, 0, 1},
   PlotStyle -> {Dashed, Red, Thickness[0.01]},
   PlotLegends -> Placed[leg2, {0.75, 0.8}]
   ];
fig4 = Show[fig2, fig3];

fig5 =
  Labeled[GraphicsRow[{fig1, fig4}, ImageSize -> 800], 
   Style[StringTemplate["g=``"][g], FontFamily -> "Times", 
    FontSize -> 24, Black], Top];
Export[StringTemplate["figs/reply/beta-g``.pdf"][g], fig5]
Show[fig5]


(* Q6. KL vs Nsites, @SWSSB and @Symmetric respectively. *)
(* Symmetric. g=10 p=0.3 *)
nSitesLst = {4, 5, 6, 7, 8};

g = 10;
p = 0.3;

{t1, data1} = Table[
    Import[
     StringTemplate[
       "data-reply/xs0-g``-p``-nSites``-nSamples10000000.mx"][g, p, 
      nSites, nSamples]];
    Import[
     StringTemplate[
       "data-reply/xs1-g``-p``-nSites``-nSamples10000000.mx"][g, p, 
      nSites, nSamples]];
    
    xs0 = xs0 // Flatten;
    xs1 = xs1 // Flatten;
    
    keys = xs0 // Counts // Keys // Sort;
    p0 = Table[Counts[xs0][k], {k, keys}] // N;
    p0 /= Total[p0];
    p1 = Table[Counts[xs1][k], {k, keys}] // N;
    p1 /= Total[p1];
    relentrp = p0*Log[p0/p1] // Total;
    
    {nSites, relentrp}
    
    , {nSites, nSitesLst}] // AbsoluteTiming;

{t1, data1}

(* SWSSB. g=10 p=0.7 *)
nSitesLst = {4, 5, 6, 7, 8};

g = 10;
p = 0.7;

{t2, data2} = Table[
    Import[
     StringTemplate[
       "data-reply/xs0-g``-p``-nSites``-nSamples10000000.mx"][g, p, 
      nSites, nSamples]];
    Import[
     StringTemplate[
       "data-reply/xs1-g``-p``-nSites``-nSamples10000000.mx"][g, p, 
      nSites, nSamples]];
    
    xs0 = xs0 // Flatten;
    xs1 = xs1 // Flatten;
    
    keys = xs0 // Counts // Keys // Sort;
    p0 = Table[Counts[xs0][k], {k, keys}] // N;
    p0 /= Total[p0];
    p1 = Table[Counts[xs1][k], {k, keys}] // N;
    p1 /= Total[p1];
    relentrp = p0*Log[p0/p1] // Total;
    
    {nSites, relentrp}
    
    , {nSites, nSitesLst}] // AbsoluteTiming;

{t2, data2}

data = Transpose[data2][[2]]/Transpose[data2][[1]];
data = Transpose[{nSitesLst, data}];

model = LinearModelFit[{#1, Log[#2]} & @@@ data, {x}, x]
fig1 = ListLinePlot[data,
   
   Axes -> False,
   AspectRatio -> aspratio,
   Frame -> True,
   FrameLabel -> (Style[#, Black, FontFamily -> "Times", 
        FontSize -> 
         16] & /@ {"\!\(\*SubscriptBox[\(N\), \(sites\)]\)", 
       "KL ratio"}),
   FrameStyle -> Thickness[0.004],
   ImageSize -> 400,
   LabelStyle -> {16, Black, FontFamily -> "Times"},
   PlotMarkers -> "OpenMarkers",
   PlotStyle -> {Blue, Thickness[0.01]},
   PlotLabel -> StringTemplate["SWSSB/Symmetric KL ratio "][g, p],
   
   ScalingFunctions -> "Log"
   ];
fig2 = Plot[model[x], {x, 3, 9}, 
   PlotStyle -> {Dashed, Red, Thickness[0.01]}];
fig3 = Show[fig1];
Export["figs/reply/Q6-SSratio.pdf", fig3]
Show[fig3]


