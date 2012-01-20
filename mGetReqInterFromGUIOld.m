%mGetReqInterFromGUI.m

function Inter = mGetReqInterFromGUI(str)

if strmatch(str,'any')%just to save time, because most Matrix boxes in the GUI will be like this
    Inter=[];
else
    BPcodes= {'cWW','tWW','cWH','cHW','tWH','tHW','cWS','cSW','tWS','tSW','cHH','tHH','cHS','cSH','tHS','tSH','cSS','cSs','csS','tSS','tSs','tsS','sP','sA'};
    BPequiv1=[1    ,2    ,3    ,-3   ,4    ,-4   ,5    ,-5   ,6    ,-6   ,7    ,8    ,9    ,-9   ,10   ,-10  ,11   ,11   ,-11  ,12   ,12   ,-12  ,15, 16   ];
    BPequiv2=[-1   ,-2   ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,-7   ,-8   ,0    ,0    ,0    ,0    ,-11  ,0    ,0    ,-12  ,0    ,0    ,17, 18   ];
    BPequiv3=[0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0,   0     ];
    BPequiv4=[0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0    ,0,   0     ];
    
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
            if (strcmp('cSs',InterStr))||(strcmp('csS',InterStr))||(strcmp('tSs',InterStr))||(strcmp('tsS',InterStr)) %only these are case sensitive, and tss will not give an error either
            ind=find(strcmp(BPcodes,InterStr));
            else
            ind=find(strcmpi(BPcodes,InterStr));
            if ~isempty(ind)
            ind=ind(1);%in case there are several hits (as in cSS and tSS)
            end
            end
            if ~isempty(ind)
                count=count+1;
                Inter(count)=BPequiv1(ind);
                if BPequiv2(ind)~=0
                    count=count+1;
                    Inter(count)=BPequiv2(ind);
                end
                if BPequiv3(ind)~=0
                    count=count+1;
                    Inter(count)=BPequiv3(ind);
                end
                if BPequiv4(ind)~=0
                    count=count+1;
                    Inter(count)=BPequiv4(ind);
                end
            else Inter=[]; %This happens if a starnge word is entered, this will be equivalet to "any" interaction
            end
        else
            count=count+1;
            Inter(count)=InterNum;
        end
    end
end
