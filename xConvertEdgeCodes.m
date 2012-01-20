% xConvertEdgeCodes.m

function Inter = mGetReqInterFromGUI(str)

if strmatch(str,'any')%just to save time, because most Matrix boxes in the GUI will be like this
    Inter=[];
else
    BPcodes={'cWW','tWW','cWH','tWH','cWS','tWS','cHH','tHH','cHS','tHS','cSS','tSS'};
    BPequiv1=[1,2,3,4,5,6,7,8,9,10,11,12];
    BPequiv2=[-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12];

    str=regexprep(str,';| ',',');
    commas=strfind(str,',');
    lim=[0 commas length(str)+1];

    count=0;
    for i=1:length(lim)-1
        InterStr=str(lim(i)+1:lim(i+1)-1);
        InterNum=str2num(InterStr);
        if isempty(InterNum)
            %         Inter{i}=InterStr;
            %%%make BPcode into number
            ind=find(strcmpi(BPcodes,InterStr));
            if ~isempty(ind)
                count=count+1;
                Inter(count)=BPequiv1(ind);
                count=count+1;
                Inter(count)=BPequiv2(ind);
            else Inter=[]; %This happens if a starnge word is entered, this will be equivalet to "any" interaction
            end
        else
            count=count+1;
            Inter(count)=InterNum;
        end
    end
end
