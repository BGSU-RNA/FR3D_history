%mGetBaseFromNTnumberAndChain.m
%File = zGetNTData('5S_Hm_Ec',0);
%Bas=mGetBaseFromNTnumberAndChain(File,'2','A')

% should be replaced by zIndexLookup

function Bas = mGetBaseFromNTnumberAndChain(File,NTnum,cha)

N=[];
for i=1:length(File.NT)
    if strcmp(File.NT(i).Number,NTnum)
        N=[N i];
    end
end

for i=1:length(N)
    if strcmp(File.NT(N(i)).Chain,cha)
        ind=N(i);
        break
    end
end

% if isempty(ind)
%     fprintf('%s\n','No match for NT number and NT chain requested');
% end

Bas=File.NT(ind).Base;
