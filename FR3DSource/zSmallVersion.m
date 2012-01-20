% zSmallVersion(File) removes some fields from File 

function [File] = zSmallVersion(File)

for f=1:length(File),
  File(f).Pair = [];
  File(f).CI   = [];
  File(f).SizeCode = 1;

  for n=1:length(File(f).NT),
    File.NT(n).Loc = [];
  end
end
