SetDirectory[NotebookDirectory[]];
Directory[]

SetEnvironment[
  "PATH" -> "/opt/homebrew/bin:" <> Environment["PATH"]];
Environment["PATH"]

data = RunProcess[{"julia", "-e", 
    "using JLD2; data=load(\"data/gs-avgzizjs2.jld2\"); \
println(data)"}];
Print[data["StandardOutput"]];


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
