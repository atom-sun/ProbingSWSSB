SetDirectory[NotebookDirectory[]];

<< "SpinOneHalfChain.wl";
<< "OpenSystem.wl";
<< "QIsingED.wl";


ColorData[97, "ColorList"]
{color1, color2} = ColorData[97, "ColorList"][[;; 2]]

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
imgsize = 360;


nSamples = 1000000;


g = 10;
p = 0.3;
mark = "\[FivePointedStar]";

resED = Table[
      {e0, v0} = qIsingED[nSites, g, "obc"];
      rho0 = KroneckerProduct[v0, v0];
      
      {j1, j2} = {1, nSites};
   
       rho1 = zzChannelA2Adm[nSites, p, rho0];
       nm = Tr[rho1 . rho1];
   x0 = Tr[rho1 . rho1];
   x1 = Tr[
     rho1 . zizjOp[nSites, j1, j2] . rho1 . zizjOp[nSites, j1, j2]];
       {x0, x1}
      , {nSites, {4, 5, 6, 7}}];

resED


data = Table[
   Import[
    StringTemplate["dump/xs0-g``-p``-nSites``-nSamples``.mx"][g, p, 
     nSites, nSamples]];
   Import[
    StringTemplate["dump/xs1-g``-p``-nSites``-nSamples``.mx"][g, p, 
     nSites, nSamples]];
   
   Assert[{nSamples, 1, 1} == Dimensions[xs0]];
   Assert[Dimensions[xs1] == Dimensions[xs0]];
   
   xs0 = xs0 // Flatten;
   xs1 = xs1 // Flatten;
   
   meanxs0 = Mean[xs0];
   stdxs0 = Sqrt[Variance[xs0]/Length[xs0]];
   
   meanxs1 = Mean[xs1];
   stdxs1 = Sqrt[Variance[xs1]/Length[xs1]];
   
   meanres = meanxs1/meanxs0;
   stdres = stdxs1/meanxs0 + stdxs0*meanxs1/meanxs0^2;
   
   {{nSites, Around[meanxs0, stdxs0]}, {nSites, 
     Around[meanxs1, stdxs1]}}
   , {nSites, 4, 7}];

{data0, data1} = data // Transpose


fig1 = ListPlot[
   {data0, data1},
   
   AspectRatio -> aspratio,
   Axes -> False,
   Frame -> True,
   FrameLabel -> {"\!\(\*SubscriptBox[\(N\), \(site\)]\)", 
     "\!\(\*SubscriptBox[\(P\), \(\\\ \)]\) "},
   FrameStyle -> Thickness[0.004], 
   FrameTicks -> {{{0, 0.1, 0.2, 0.3, 0.4}, None}, {{4, 5, 6, 7}, 
      None}},
   GridLinesStyle -> Directive[Gray, Dashed],
   GridLines -> {{4, 5, 6, 7}, {0, 0.1, 0.2, 0.3, 0.4}},
   ImageSize -> imgsize,
   LabelStyle -> {24, Black, FontFamily -> "Times"},
   PlotMarkers -> "OpenMarkers",
   PlotRange -> {{3.5, 7.5}, {-0.03, 0.48}},
   PlotStyle -> {{Thickness[0.008], color1}, {Thickness[0.008], color2}},

   PlotLegends -> Placed[
     leg1,
     {0.85, 0.8}]
   ];

dataED = Thread[{Range[4, 7], #}] & /@ Transpose[resED];
fig2 = ListLinePlot[dataED,
   PlotStyle -> {{color1, Dashed, Thickness[0.006]}, {color2, Dashed, 
      Thickness[0.006]}}
   ];
fig0 = ListPlot[{{{3.7, 0.45}}}, PlotMarkers -> {{mark, 20}}, 
   PlotStyle -> Red];
fig3 = Show[fig1, fig2, fig0]
Export[StringTemplate["figs/measurebar_randm-g``-p``-nSamples``.pdf"][
   g, p, nSamples], fig3];


g = 10;
p = 0.7;
mark = "\[SixPointedStar]";

resED = Table[
      {e0, v0} = qIsingED[nSites, g, "obc"];
      rho0 = KroneckerProduct[v0, v0];
      
      {j1, j2} = {1, nSites};
   
       rho1 = zzChannelA2Adm[nSites, p, rho0];
       nm = Tr[rho1 . rho1];
   x0 = Tr[rho1 . rho1];
   x1 = Tr[
     rho1 . zizjOp[nSites, j1, j2] . rho1 . zizjOp[nSites, j1, j2]];
       {x0, x1}
      , {nSites, {4, 5, 6, 7}}];

resED


data = Table[
   Import[
    StringTemplate["dump/xs0-g``-p``-nSites``-nSamples``.mx"][g, p, 
     nSites, nSamples]];
   Import[
    StringTemplate["dump/xs1-g``-p``-nSites``-nSamples``.mx"][g, p, 
     nSites, nSamples]];
   
   Assert[{nSamples, 1, 1} == Dimensions[xs0]];
   Assert[Dimensions[xs1] == Dimensions[xs0]];
   
   xs0 = xs0 // Flatten;
   xs1 = xs1 // Flatten;
   
   meanxs0 = Mean[xs0];
   stdxs0 = Sqrt[Variance[xs0]/Length[xs0]];
   
   meanxs1 = Mean[xs1];
   stdxs1 = Sqrt[Variance[xs1]/Length[xs1]];
   
   meanres = meanxs1/meanxs0;
   stdres = stdxs1/meanxs0 + stdxs0*meanxs1/meanxs0^2;
   
   {{nSites, Around[meanxs0, stdxs0]}, {nSites, 
     Around[meanxs1, stdxs1]}}
   , {nSites, 4, 7}];

{data0, data1} = data // Transpose


fig1 = ListPlot[
   {data0, data1},
   
   AspectRatio -> aspratio,
   Axes -> False,
   Frame -> True,
   FrameLabel -> {"\!\(\*SubscriptBox[\(N\), \(site\)]\)", 
     "\!\(\*SubscriptBox[\(P\), \(\\\ \)]\) "},
   FrameStyle -> Thickness[0.004], 
   FrameTicks -> {{{0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4}, 
      None}, {{4, 5, 6, 7}, None}},
   GridLinesStyle -> Directive[Gray, Dashed],
   GridLines -> {{4, 5, 6, 7}, {0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4}},
   ImageSize -> imgsize,
   LabelStyle -> {24, Black, FontFamily -> "Times"},
   PlotMarkers -> "OpenMarkers",
   PlotRange -> {{3.5, 7.5}, {-0.03, 0.22}},
   PlotStyle -> {{Thickness[0.008], color1}, {Thickness[0.008], color2}},

   PlotLegends -> Placed[
     leg1,
     {0.85, 0.8}]
   ];

dataED = Thread[{Range[4, 7], #}] & /@ Transpose[resED];
fig2 = ListLinePlot[dataED,
   PlotStyle -> {{color1, Dashed, Thickness[0.006]}, {color2, Dashed, 
      Thickness[0.006]}}
   ];
fig0 = ListPlot[{{{3.7, 0.205}}}, PlotMarkers -> {{mark, 20}}, 
   PlotStyle -> Red];
fig3 = Show[fig1, fig2, fig0]
Export[StringTemplate["figs/measurebar_randm-g``-p``-nSamples``.pdf"][
   g, p, nSamples], fig3];


nSites = 7;


g = 10;
p = 0.3;
mark = "\[FivePointedStar]";

Import[StringTemplate["dump/xs0-g``-p``-nSites``-nSamples``.mx"][g, p,
    nSites, nSamples]];
Import[StringTemplate["dump/xs1-g``-p``-nSites``-nSamples``.mx"][g, p,
    nSites, nSamples]];

xs0 = xs0 // Flatten;
xs1 = xs1 // Flatten;

keys = Table[2.^nSites*(-2)^(-i), {i, 0, 7}];
p0 = Table[Counts[xs0][k], {k, keys}] // N;
p0 /= Total[p0]
p1 = Table[Counts[xs1][k], {k, keys}] // N;
p1 /= Total[p1]

p0*Log[p0/p1] // Total


fig1 = BarChart[Transpose[{p0, p1}],
   AspectRatio -> aspratio,
   BarSpacing -> {Large, 2},
   Frame -> True,
   FrameLabel -> {Column[{Style["D(s, s\[Prime])", 24, Italic, 
     FontFamily -> "Times"]}, BaselinePosition -> Top], 
     Style["p(D)", 24, Italic, FontFamily -> "Times"]},
   FrameStyle -> Thickness[0.004],
   GridLines -> Automatic,
   GridLinesStyle -> Directive[Gray, Dashed],
   ImageSize -> imgsize,
   LabelStyle -> {24, Black, FontFamily -> "Times"},
   PlotRange -> {Automatic, 0.3},
   
   ChartLabels -> {Placed[{"0", "1", "2", "3", "4", "5", "6", "7"}, 
      Below], Automatic},
   ChartLegends -> Placed[
     leg1,
     {0.85, 0.8}],
   ChartStyle -> {color1, color2}
   ];
fig0 = ListPlot[{{{1.5, 0.28}}}, PlotMarkers -> {{mark, 20}}, 
   PlotStyle -> Red];
fig3 = Show[fig1, fig0]
Export[StringTemplate["figs/hist_randm-g``-p``-nSamples``.pdf"][g, p, 
   nSamples], fig3];


g = 10;
p = 0.7;
mark = "\[SixPointedStar]";

Import[StringTemplate["dump/xs0-g``-p``-nSites``-nSamples``.mx"][g, p,
    nSites, nSamples]];
Import[StringTemplate["dump/xs1-g``-p``-nSites``-nSamples``.mx"][g, p,
    nSites, nSamples]];

xs0 = xs0 // Flatten;
xs1 = xs1 // Flatten;

keys = Table[2.^nSites*(-2)^(-i), {i, 0, 7}];
p0 = Table[Counts[xs0][k], {k, keys}] // N;
p0 /= Total[p0]
p1 = Table[Counts[xs1][k], {k, keys}] // N;
p1 /= Total[p1]

p0*Log[p0/p1] // Total


fig1 = BarChart[Transpose[{p0, p1}],
   AspectRatio -> aspratio,
   BarSpacing -> {Large, 2},
   Frame -> True,
   FrameLabel -> {Column[{Style["D(s, s\[Prime])", 24, Italic, 
     FontFamily -> "Times"]}, BaselinePosition -> Top], 
     Style["p(D)", 24, Italic, FontFamily -> "Times"]},
   FrameStyle -> Thickness[0.004],
   GridLines -> Automatic,
   GridLinesStyle -> Directive[Gray, Dashed],
   ImageSize -> imgsize,
   LabelStyle -> {24, Black, FontFamily -> "Times"},
   PlotRange -> {Automatic, 0.3},
   
   ChartLabels -> {Placed[{"0", "1", "2", "3", "4", "5", "6", "7"}, 
      Below], Automatic},
   ChartLegends -> Placed[
     leg1,
     {0.85, 0.8}],
   ChartStyle -> {color1, color2}
   ];
fig0 = ListPlot[{{{1.5, 0.28}}}, PlotMarkers -> {{mark, 20}}, 
   PlotStyle -> Red];
fig3 = Show[fig1, fig0]
Export[StringTemplate["figs/hist_randm-g``-p``-nSamples``.pdf"][g, p, 
   nSamples], fig3];

