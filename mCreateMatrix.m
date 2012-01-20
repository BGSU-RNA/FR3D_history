% mCreateMatrix is a script that sets up the interaction matrix for FR3D_GUI

% modify the nucleotide index lookup to be aware of the given chain

%%%%%%Create the basepair matrix (and delete extra ones) for determining
%%%%%%Query.Diff, Query.ReqInter, Query.Diagonal
%%%%%%Query.Config added

PreviousConfig=ones(1,12);
for i=1:12
    h=findobj('Tag',strcat('Config',num2str(i)));
    try,PreviousConfig(i)=get(h,'Value');end
    delete(h);
end

for i=1:length(NT)
    str={'','anti','syn'};
    handles.Config(i) = uicontrol('Tag',strcat('Config',num2str(i)),'Style','popupmenu','Units','normalized','Position',[(0.25+0.057*i) (0.755) .054 .04],'Background',[1 1 1],'String',str,'Value',PreviousConfig(i));
end

for i=1:12
    hh=findobj('Tag',strcat('BPtexth',num2str(i)));
    delete(hh);
    hv=findobj('Tag',strcat('BPtextv',num2str(i)));
    delete(hv);
end

for i=1:length(NT)
    if get(handles.Geometric,'Value') == 1
        ind = zIndexLookup(File(QIndex),NT{i},ChainList{i});
        Bas = File(QIndex).NT(ind).Base;
    else
        Bas='NT';
    end
    %%%Horizontal title line:
%     handles.BPtexth(i) = uicontrol('Tag',strcat('BPtexth',num2str(i)),'Style','text','Units','normalized','Position',[(0.25+0.057*i) (0.73) .054 .04],'String',strcat(File(QIndex).NT(str2num(NT{i})).Base,NT(i)));
    handles.BPtexth(i) = uicontrol('Tag',strcat('BPtexth',num2str(i)),'Style','text','Units','normalized','Position',[(0.25+0.057*i) (0.71) .054 .04],'String',strcat(Bas,NT(i)));
    %%%Vertical title line:
    handles.BPtextv(i) = uicontrol('Tag',strcat('BPtextv',num2str(i)),'Style','text','Units','normalized','Position',[(0.26) (0.715-0.0445*i) .054 .04],'String',strcat(Bas,NT(i)));
end

%%%Now the matrix:
for i=1:12
    hd=findobj('Tag',strcat('Diagonal',num2str(i)));
    try,PreviousDiagonal{i}=get(hd,'String');end
    delete(hd);
    for j=1:12
        hh=findobj('Tag',strcat('BPType',num2str(i),num2str(j)));
        try,PreviousBPType{i,j}=get(hh,'String');end
        delete(hh);
        hv=findobj('Tag',strcat('Diff',num2str(i),num2str(j)));
        try,PreviousDiff{i,j}=get(hv,'String');end
        delete(hv);
    end
end

% % % if ~isempty(PreviousMasks)
% % %     PreviousMasks(1)
% % %     PreviousBPType{1,2}
% % %     PreviousDiff{2,1}
% % % end

%MaskList={'N (ACGU)','A','C','G','U','R (AG)','Y (CU)','M (AC)','W (AU)','S (GC)','K (GU)','V (ACG)','H (ACU)','D (AGU)','B (CGU)'};
for i=1:length(NT)
    for j=1:length(NT)
        if j==i
            if isempty(PreviousDiagonal{i}),PreviousDiagonal{i}='N';end
            handles.Diagonal(i) = uicontrol('Tag',strcat('Diagonal',num2str(i)),'Style','edit','Units','normalized','Position',[(0.25+0.057*j) (0.73-0.045*i) .054 .04],'Background',[1 1 1],'String',PreviousDiagonal{i});
        end
        if i<j
            if isempty(PreviousBPType{i,j}),PreviousBPType{i,j}='';end
            handles.BPType(i,j) = uicontrol('Tag',strcat('BPType',num2str(i),num2str(j)),'Style','edit','Units','normalized','Position',[(0.25+0.057*j) (0.73-0.045*i) .054 .04],'String',PreviousBPType{i,j},'Background',[1 1 0]);
        elseif i>j
            if isempty(PreviousDiff{i,j}),PreviousDiff{i,j}='';end
            handles.Diff(i,j) = uicontrol('Tag',strcat('Diff',num2str(i),num2str(j)),'Style','edit','Units','normalized','Position',[(0.25+0.057*j) (0.73-0.045*i) .054 .04],'String',PreviousDiff{i,j},'Background',[0 1 1]);
        end
    end
end

set(handles.ConfigText,'Visible','on');
set(handles.NTmaskText,'Visible','on');
set(handles.MaxDistText,'Visible','on');
set(handles.InteractionText,'Visible','on');
