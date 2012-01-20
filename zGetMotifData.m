
function [M] = zGetMotifData(M)

File = zGetNTData(M.Filename,0);

Numbers = cat(1,{File.NT(:).Number});
for i=1:length(M.Bases)
    M.Indices(i)= find(ismember(Numbers,M.Bases(i)));
    M.NT(i)     = File.NT(M.Indices(i));
end
for i=1:length(M.AllBases)
    M.AllIndices(i) = find(ismember(Numbers,M.AllBases(i)));
    M.AllNT(i)      = File.NT(M.AllIndices(i));
end
M.NumNT=length(M.NT);

