%mSpecifyQuery_forGUI
% The variable Query has several fields:
% xSpecifyQuery returns the description of a model motif
% The variable Query has several fields:
%   Query.Filename       a string like rr0033_23S
%   Query.NTList         a list of nucleotide numbers
%   Query.Mask           a mask for which nucleotides to allow
%   Query.AngleWeight    weights to put on the angles
%   Query.DistanceWeight weights to put on nucleotide distances
%   Query.DiscCutoff     discrepancy cutoff D_0
%   Query.RelCutoff      relaxed cutoff D_1
%   Query.Diff        maximum difference between sorted nucleotide
%     indices; use with Sequential = 1 to screen for sequential constraints
%   Query.Geometric      set to 0 to ignore geometry, only use screens. 
%                        Default is 1.
%   Query.ExcludeOverlap set to 1 to eliminate highly redundant motifs with
%     larger discrepancies; often, you get the same candidate with one or two
%     nucleotides different but much higher discrepancy.  Default is 1 when
%     Query.NumNT > 6.
%Added: Query.Config

for i=1:NTlen
    h=findobj('Tag',strcat('Config',num2str(i)));
    v=get(h,'Value');
    s=get(h,'String');
    Query.Config(i)=s(v); %use curly baracelets here {} because Query.Mask must be a character array, not a cell array. The (1) means that only first letter is taken (what is after it may be just description)
end
%%%End Read Query.ChainList from GUI

DiscCutoff=str2num(get(handles.GuarCutoff,'String')); %%%GuarCutoff name changed into DiscCutoff
if ~isempty(DiscCutoff)
   Query.DiscCutoff = DiscCutoff;
end

RelCutoff=str2num(get(handles.RelCutoff,'String'));
if ~isempty(RelCutoff)
   Query.RelCutoff = RelCutoff;
end

v=get(handles.Overlap,'Value');
if v==1
    Query.ExcludeOverlap=1;
else
    Query.ExcludeOverlap=0;
end

Query.Name=get(handles.SearchName,'String');
desc=get(handles.SearchDescription,'String');
if ~isempty(desc) && ~strcmp(desc,'-')
    Query.Description=desc;
end

if get(handles.Geometric,'Value') == 1
    Query.Geometric=1;
else
    Query.Geometric=0;
end


% % %%%Read Query.ChainList from GUI  %%%%Now this is done just once before the Matrix is created, and then exported to this function
% % % Query.ChainList      = {'9' '9' '9'};
% % for i=1:NTlen %length(Query.NTList)
% %     h=findobj('Tag',strcat('ChainList',num2str(i)));
% %     s=get(h,'String');
% %     v=get(h,'Value');
% %     Query.ChainList(i)=s(v);
% % end
% % %%%End Read Query.ChainList from GUI


%%%Read Query.Diagonal from GUI
% Query.Diagonal includes Mask information and also some other info
for i=1:NTlen
    h=findobj('Tag',strcat('Diagonal',num2str(i)));
    s=get(h,'String');
    Query.Diagonal{i}=s;
end
%%%End Read Query.Diagonal from GUI


%%%Read Query.Diff from GUI
% Query.Diff        = [4 4];
for i=1:NTlen
    for j=1:NTlen
        if i>j
            h=findobj('Tag',strcat('Diff',num2str(i),num2str(j)));
            s=get(h,'String');
            Query.Diff{i,j}=s;
%             v=str2num(s);
%             if ~isempty(v)
%                 Query.Diff(i,j)=round(v);
%             else
%                 Query.Diff(i,j)=Inf;
%             end
%             if Query.Diff(i,j)==0;
%                 Query.Diff(i,j)=Inf;
%             end
        end
    end
end

for i=1:NTlen
    for j=1:NTlen
        if i<j
            h=findobj('Tag',strcat('BPType',num2str(i),num2str(j)));
            s=get(h,'String');
            %Query.ReqInter{i,j}=mGetReqInterFromGUI(s);
            Query.Edges{i,j}= regexprep(s,',',' '); %users can enter , instead of spaces between interaction codes
            
% % %             [Query.ReqEdge{i,j},Query.ExEdge{i,j},Query.OKPairs{i,j},Query.ExPairs{i,j}]=...
% % %                 xInterpretEdgeText(Query.Edges{i,j});
        end
    end
end


