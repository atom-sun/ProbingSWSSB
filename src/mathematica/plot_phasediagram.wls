SetDirectory[NotebookDirectory[]];
FileNames["*", Directory[]];

SetEnvironment[
  "PATH" -> "/opt/homebrew/bin:" <> Environment["PATH"]];
Environment["PATH"]


(* anal fit with julia MPS data of N=200 system size. *)
data = RunProcess[{"julia", "-e", 
    "using JLD2; data=load(\"data/gs-avgzizjs2.jld2\"); \
println(data)"}];
(*Print[data["StandardOutput"]];*)

out = data["StandardOutput"];
StringPosition[out, "\"sumzizjs2\""]
StringPosition[out, "\"params\""]
StringPosition[out, "\"sumzizjs\""]

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

{sumzizjs2, paramgs}

ddata0 = Thread[{paramgs, sumzizjs2/100}];
fun0 = Interpolation[ddata0];


(* plot phase diagram. sim_randm density plot vs. anal fit boundary. *)
pps = Range[0., 0.75, 0.075]
ggs0 = Range[1/4., 1/2., 1/41.]
ggs = Tan[Pi*ggs0]

nSites = 7;
nSamples = 10000;

data = Table[
   g = Tan[Pi*g0];
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
   
   keys = xs0 // Counts // Keys // Sort;
   p0 = Table[Counts[xs0][k], {k, keys}] // N;
   p0 /= Total[p0];
   p1 = Table[Counts[xs1][k], {k, keys}] // N;
   p1 /= Total[p1];
   relentrp = p0*Log[p0/p1] // Total;
   
   {g0, p, relentrp}
   , {g0, ggs0}, {p, pps}];
data1 = Flatten[data, 1];


fig1 = ListDensityPlot[data1,
   ColorFunction -> "BlueGreenYellow",
   Frame -> True,
   FrameLabel -> {"arctan(g)/\[Pi]", 
     "\!\(\*SubscriptBox[\(\[Mu]\), \(\\\ \)]\)"},
   FrameStyle -> Thickness[0.004],
   InterpolationOrder -> 0,
   LabelStyle -> {20, Black, FontFamily -> "Times"},
   MeshStyle -> Opacity[0.5, Black],
   PlotLegends -> 
    Placed[BarLegend[Automatic, LegendMarkerSize -> 320, 
      LegendLabel -> 
       Placed["KL divergence", Right, 
        Rotate[#, 90 Degree] &]], {After, Top}],
   ImageSize -> 360
   ];
fig2 = Plot[
   dataal = 1./(2. + fun0[Tan[Pi*x]]), {x, 1/4 + 0.003, 1/2 - 0.005},
   PlotRange -> All,
   PlotStyle -> {Thick, LightGray, Dashed}
   ];
fig0 = ListPlot[{{{ArcTan[10]/Pi, 0.3}}, {{ArcTan[10]/Pi, 0.7}}}, 
   PlotMarkers -> {{"\[FivePointedStar]", 20}, {"\[SixPointedStar]", 
      20}}, PlotStyle -> Red];
txt1 = Graphics[
   Text[Style["SW-SSB", 20, Bold, White, 
     FontFamily -> "Times"], {0.32, 0.57}]];
txt2 = Graphics[
   Text[Style["Symmetric", 20, Bold, White, 
     FontFamily -> "Times"], {0.38, 0.25}]];
fig3 = Show[fig1, fig2, fig0, txt1, txt2]
Export["figs/phasediagram.png", fig3, "PNG"]


(* choose color function below. *)
cols = ColorData["Gradients"];

Table[
 fig1 = ListDensityPlot[data1,
   ColorFunction -> c,
   Frame -> True,
   FrameLabel -> {"arctan(g)/\[Pi]", 
     "\!\(\*SubscriptBox[\(p\), \(\\\ \)]\)"},
   FrameStyle -> Thickness[0.004],
   InterpolationOrder -> 0,
   LabelStyle -> {20, Black, FontFamily -> "Times"},
   MeshStyle -> Opacity[0.5, Black],
   PlotLegends -> Automatic,
   
   PlotLabel -> c,
   ImageSize -> 360
   ];
 fig2 = Plot[
   1./(2. + fun0[Tan[Pi*x]]), {x, 1/4 + 0.003, 1/2 - 0.005},
   PlotRange -> All,
   PlotStyle -> {Thick, Black, Dashed}
   ];
 fig3 = Show[fig1, fig2];
 Export["figs/colorfunc/" <> c <> ".png", fig3];
 fig3
 , {c, cols}
 ]
