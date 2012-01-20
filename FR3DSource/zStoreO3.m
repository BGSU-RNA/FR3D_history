
function [File] = zStoreO3(File)

for f = 1:length(File),
  File(f).NT(1).Sugar(13,:) = [Inf Inf Inf];
  for i = 2:length(File(f).NT),
    File(f).NT(i).Sugar(13,:) = File(f).NT(i-1).Sugar(5,:);
  end
end