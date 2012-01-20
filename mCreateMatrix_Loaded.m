%mCreateMatrix_Loaded

%%%%%%Create the basepair matrix (and delete extra ones) for determining
%%%%%%Query.Diff, Query.Edges, Query.Diagonal
%%%%%%Query.Config added
for i=1:12
    h=findobj('Tag',strcat('Config',num2str(i)));
    delete(h);
end
for i=1:length(NT)
    str={'','anti','syn'};
    if isfield(Search.Query,'Config'),
      switch Search.Query.Config{i}
        case 'anti' 
          PreviousConfig(i) = 2;
        case 'syn'  
          PreviousConfig(i) = 3;
        otherwise   
          PreviousConfig(i) = 1;
      end
    else
      PreviousConfig(i) = 1;
    end
    handles.Config(i) = uicontrol('Tag',strcat('Config',num2str(i)),'Style','popupmenu','Units','normalized','Position',[(0.25+0.057*i) (0.755) .054 .04],'Background',[1 1 1],'String',str,'Value',PreviousConfig(i));
end

for i=1:12
    hh=findobj('Tag',strcat('BPtexth',num2str(i)));
    delete(hh);
    hv=findobj('Tag',strcat('BPtextv',num2str(i)));
    delete(hv);
    hd=findobj('Tag',strcat('Diagonal',num2str(i)));
    delete(hd);
    for j=1:12
        hh=findobj('Tag',strcat('BPType',num2str(i),num2str(j)));
        delete(hh);
        hv=findobj('Tag',strcat('Diff',num2str(i),num2str(j)));
        delete(hv);
    end
end
for i=1:length(NT)
    if isfield(Search.Query,'NT')
        Bas=Search.Query.NT(i).Base;
        Num=Search.Query.NT(i).Number;
        str=strcat(Bas,Num);
    else
        Bas='NT';
        str=strcat(Bas,num2str(i));
    end
    handles.BPtexth(i) = uicontrol('Tag',strcat('BPtexth',num2str(i)),'Style','text','Units','normalized','Position',[(0.25+0.057*i) (0.71) .054 .04],'String',str);
    handles.BPtextv(i) = uicontrol('Tag',strcat('BPtextv',num2str(i)),'Style','text','Units','normalized','Position',[(0.26) (0.715-0.0445*i) .054 .04],'String',str);
end

%%%Now the matrix:
%MaskList={'N (ACGU)','A','C','G','U','R (AG)','Y (CU)','M (AC)','W (AU)','S (GC)','K (GU)','V (ACG)','H (ACU)','D (AGU)','B (CGU)'};
PreviousMasks=ones(1,12); %declared
% PreviousDiff={''};
% PreviousBPType={''};

for i=1:length(NT)
    for j=1:length(NT)
        if j==i
            if isfield(Search.Query,'Diagonal')
                PreviousDiagonal{i}=Search.Query.Diagonal{i};
            else PreviousDiagonal{i}='N';
            end
            handles.Diagonal(i) = uicontrol('Tag',strcat('Diagonal',num2str(i)),'Style','edit','Units','normalized','Position',[(0.25+0.057*j) (0.73-0.045*i) .054 .04],'Background',[1 1 1],'String',PreviousDiagonal{i});
        end
        if i<j
            if isfield(Search.Query,'Edges')
                if ~isempty(Search.Query.Edges{i,j}),
                    PreviousBPType{i,j}=Search.Query.Edges{i,j};
                else
                    PreviousBPType{i,j}='';
                end
            else PreviousBPType{i,j}='';
            end
            handles.BPType(i,j) = uicontrol('Tag',strcat('BPType',num2str(i),num2str(j)),'Style','edit','Units','normalized','Position',[(0.25+0.057*j) (0.73-0.045*i) .054 .04],'String',PreviousBPType{i,j},'Background',[1 1 0]);
        elseif i>j
            if isfield(Search.Query,'Diff')
                if ~isempty(Search.Query.Diff(i,j)),
                    PreviousDiff{i,j}=Search.Query.Diff{i,j};
                else
                    PreviousDiff{i,j}='';
                end
            else PreviousDiff{i,j}='';
            end
            handles.Diff(i,j) = uicontrol('Tag',strcat('Diff',num2str(i),num2str(j)),'Style','edit','Units','normalized','Position',[(0.25+0.057*j) (0.73-0.045*i) .054 .04],'String',PreviousDiff{i,j},'Background',[0 1 1]);
        end
    end
end

set(handles.ConfigText,'Visible','on');
set(handles.NTmaskText,'Visible','on');
set(handles.MaxDistText,'Visible','on');
set(handles.InteractionText,'Visible','on');

% Search.Query.Edges
% PreviousBPType
% Search.Query.Diff
% PreviousDiff

