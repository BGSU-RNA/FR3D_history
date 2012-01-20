%mGetChainsOfNTnumber.m
%File = zGetNTData('5S_Hm_Ec',0);
%Ch=mGetChainsOfNTnumber(File,'1')
%Ch=mGetChainsOfNTnumber(File,'3')
function Ch = mGetChainsOfNTnumber(File,NTnum) %File here is the Query PDB  %%num must be a character not number

ct=1;
for i=1:length(File.NT)
    if strcmp(File.NT(i).Number,NTnum)
        Ch{ct}=File.NT(i).Chain;
        ct=ct+1;
    end
end
