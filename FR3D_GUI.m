%<div class="moz-text-flowed" style="font-family: -moz-fixed">
function varargout = FR3D_GUI(varargin)   %By Ali Mokdad - March 2006
% FR3D_GUI M-file for FR3D_GUI.fig
%      FR3D_GUI, by itself, creates a new FR3D_GUI or raises the existing
%      singleton*.
%
%      H = FR3D_GUI returns the handle to a new FR3D_GUI or the handle to
%      the existing singleton*.
%
%      FR3D_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FR3D_GUI.M with the given input arguments.
%
%      FR3D_GUI('Property','Value',...) creates a new FR3D_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FR3D_GUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FR3D_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to runsearch (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FR3D_GUI

% Last Modified by GUIDE v2.5 17-Jul-2006 19:31:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FR3D_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @FR3D_GUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FR3D_GUI is made visible.
function FR3D_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
savedf=ls('SearchSaveFiles/*.mat');
if isempty(savedf),savedf=' ';end
set(handles.LOAD,'String',savedf);

mGetPDBfilenames %defines s
set(handles.SearchPDBs,'String',s);
set(handles.SearchPDBs,'Min',1);
set(handles.SearchPDBs,'Max',length(s)+1);
set(handles.QueryPDB,'String',s);

handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = FR3D_GUI_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;






function LOAD_Callback(hObject, eventdata, handles)
savedf=ls('SearchSaveFiles/*.mat');
if isempty(savedf),savedf=' ';end
set(handles.LOAD,'String',savedf); %Update without having to reopen GUI
v=get(handles.LOAD,'Value');
f=savedf(v,1:end);
l=strcat('SearchSaveFiles/',f);
load(l)


mSetLoadedParameters

handles.Search=Search;
guidata(hObject, handles);

function QueryPDBedit_Callback(hObject, eventdata, handles)
function QueryPDBedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function QueryNTs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function QueryNTs_Callback(hObject, eventdata, handles)


% --- Executes on button press in ReadQuery.
function ReadQuery_Callback(hObject, eventdata, handles)
%%%DetermineQuery.NTList %%%%This must be done before running mCreateDynamicGUI
NTs = get(handles.QueryNTs,'String');
NTs    = regexprep(NTs,';| ',',');
while strfind(NTs,',,'),
    NTs = regexprep(NTs,',,',',');
end
ind=findstr(',',NTs);
ind=[0 ind length(NTs)+1];
for i=1:length(ind)-1
    NT{i}=NTs(ind(i)+1:ind(i+1)-1);
end
Query.NTList=NT;
%%%End DetermineQuery.NTList

if length(NT)<=12 %this is a limitation by the size of the GUI
    %%%Read Query PDB File or just use the one in memory if it is present and same as the new Query PDB
    x=get(handles.QueryPDB,'Value');
    s=get(handles.QueryPDB,'String');
    Query.Filename = s{x};
    set(handles.Status,'String','Reading query PDB, please wait ...');
    drawnow
    %%handles.File.Filename 
    %%%%%%%%%%%%%%%ALERT: CRAIG'S HELP NEEDED: I need to use File(SIndex), but how to set SIndex???????
    
    if isfield(handles,'File')
        File = handles.File;
        [File,SIndex]=zAddNTData(Query.Filename,0,File);
    else
        [File,SIndex]=zAddNTData(Query.Filename,0);
    end
    %%%

    %     mCreateDynamicGUI %Creates all dynamic popum menus and edit boxes and hides extra ones from previous searches
    mCreateChains


    %Pass some variables created here to be used by other functions
    handles.NT=NT;
    handles.Filename = Query.Filename;
    handles.File = File;
    handles.SIndex=SIndex;
    handles.NTList=Query.NTList;
    guidata(hObject, handles);

    if get(handles.ViewQuery,'Value')==1
        figure(3)
        zDisplayNT(File(SIndex),NT);
    end
    set(handles.Status,'String','Choose "Query Chains" and click "Generate Interaction Matrix"');
    set(handles.GenerateMatrix,'Visible','on');

else
    set(handles.Status,'String','Error: It is not possible to use this GUI for a motif longer than 12 NTs');
    set(handles.GenerateMatrix,'Visible','off');
end
set(handles.RunSearch,'Visible','off');
set(handles.ListCandidates,'Visible','off');
set(handles.DisplayCandidates,'Visible','off');


function GenerateMatrix_Callback(hObject, eventdata, handles)

if get(handles.Geometric,'Value') == 1
    NT=handles.NT;
    File=handles.File;
    SIndex=handles.SIndex;

    %%%Read Query.ChainList from GUI
    % Query.ChainList      = {'9' '9' '9'};
    for i=1:length(NT)
        h=findobj('Tag',strcat('ChainPopup',num2str(i)));
        s=get(h,'String');
        v=get(h,'Value');
        Query.ChainList(i)=s(v);
    end
    %%%End Read Query.ChainList from GUI
    ChainList=Query.ChainList;
    handles.ChainList=ChainList;
    handles.NTlen=length(NT);
    set(handles.Status,'String','Ready to "Search"');
else
    NTlen=str2num(get(handles.NumberOfNTs,'String'));
    if NTlen>12
        set(handles.Status,'String','12 is the maximum length of motifs that can be addressed in this GUI. Ready to "Search"');
    else
        set(handles.Status,'String','Ready to "Search"');
    end
    NTlen=min(12,NTlen);
    for i=1:NTlen
        NT{i}=num2str(i);
    end
    handles.NTlen=NTlen;
end

mCreateMatrix
set(handles.GuarCutoffText,'visible','on');
set(handles.GuarCutoff,'visible','on');
set(handles.RelCutoffText,'visible','on');
set(handles.RelCutoff,'visible','on');
set(handles.Overlap,'visible','on');
set(handles.RunSearch,'visible','on');
set(handles.SearchNameText,'visible','on');
set(handles.SearchName,'visible','on');
set(handles.SearchDescriptionText,'visible','on');
set(handles.SearchDescription,'visible','on');


guidata(hObject, handles);


% --- Executes on button press in RunSearch.
function RunSearch_Callback(hObject, eventdata, handles)
if ~isdeployed,clc,end

% if get(handles.Geometric,'Value') == 1
if isfield(handles,'File')
    File=handles.File;
    SIndex=handles.SIndex;
end
if get(handles.Geometric,'Value') == 1
    Query.Filename=handles.Filename;
    Query.ChainList=handles.ChainList;
    Query.NTList=handles.NTList;
end
NTlen=handles.NTlen;
mSpecifyQuery_forGUI;

%%Now determine PDB "Filenames" to search in:
x=get(handles.SearchPDBs,'Value');
s=get(handles.SearchPDBs,'String');
for i=1:length(x)
    Filenames{i,1}=s{x(i)};
end
%%

GUIactive = 1;

set(handles.Status,'String','Searching... Please wait (press Ctrl+C to interrupt)');
drawnow

FR3D
% Search.Query.Inter

handles.File=File;
handles.SIndex=SIndex;

if ~isempty(Candidates),
    set(handles.ListCandidates,'Visible','on');
    set(handles.DisplayCandidates,'Visible','on');
    if length(Discrepancy)>1
        s='s';
    else s='';
    end
    str=strcat(num2str(length(Discrepancy)),' candidate',s,' found');
    set(handles.Status,'String',str);

    handles.Search=Search;

else
    set(handles.Status,'String','No candidates found');
    set(handles.ListCandidates,'Visible','off');
    set(handles.DisplayCandidates,'Visible','off');
end
guidata(hObject, handles);


function SearchName_Callback(hObject, eventdata, handles)
function SearchName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SearchPDBs_Callback(hObject, eventdata, handles)
function SearchPDBs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function QueryPDB_Callback(hObject, eventdata, handles)
function QueryPDB_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ViewQuery_Callback(hObject, eventdata, handles)


function NumberOfNTs_Callback(hObject, eventdata, handles)
function NumberOfNTs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on button press in Geometric.
function Geometric_Callback(hObject, eventdata, handles)
set(handles.Geometric,'Value',1);
set(handles.NonGeometric,'Value',0);

set(handles.NumberOfNTsTitle,'Visible','off');
set(handles.NumberOfNTs,'Visible','off');
set(handles.GenerateMatrix,'Visible','off');

set(handles.QueryPDBTitle,'Visible','on');
set(handles.QueryPDB,'Visible','on');
set(handles.ViewQuery,'Visible','on');
set(handles.QueryNTsTitle,'Visible','on');
set(handles.QueryNTs,'Visible','on');
set(handles.ReadQuery,'Visible','on');

set(handles.Status,'String','Hint: Choose "Query PDB" and "Query NTs", then click "Read Query" - Repeat anytime to restart');
Query.Geometric = 1;


% --- Executes on button press in NonGeometric.
function NonGeometric_Callback(hObject, eventdata, handles)
set(handles.NonGeometric,'Value',1);
set(handles.Geometric,'Value',0)

set(handles.QueryPDBTitle,'Visible','off');
set(handles.ViewQuery,'Visible','off');
set(handles.QueryNTsTitle,'Visible','off');
set(handles.ReadQuery,'Visible','off');

set(handles.QueryPDB,'Visible','off');
set(handles.QueryNTs,'Visible','off');
% h=findobj('Tag','QueryPDB');delete(h);
% h=findobj('Tag','QueryNTs');delete(h);
for i=1:12
    h=findobj('Tag',strcat('ChainPopup',num2str(i)));
    delete(h);
end
set(handles.QueryChains,'Visible','off')

set(handles.NumberOfNTsTitle,'Visible','on');
set(handles.NumberOfNTs,'Visible','on');
set(handles.GenerateMatrix,'Visible','on');

set(handles.Status,'String','Hint: Input "Number of NTs" (positive integer) and click "Generate Interaction Matrix" - Repeat anytime to restart');
Query.Geometric = 0;



function SearchDescription_Callback(hObject, eventdata, handles)
function SearchDescription_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function DisplayCandidates_Callback(hObject, eventdata, handles)

Search=handles.Search;
if isfield(handles,'SIndex')
    File=handles.File;
    SIndex=handles.SIndex;

    xDisplayCandidates(File(SIndex),Search);

else %If data is loaded from saved search results
    [File,SIndex]=zAddNTData(Search.Filenames,0);
    handles.File=File;
    handles.SIndex=SIndex;
    guidata(hObject, handles);

    xDisplayCandidates(File(SIndex),Search);   %%%%%ALERT: SIndex here is not correct?
end
%xGroupCandidates(File(SIndex),Search);  % doesn't work very well yet!



function ListCandidates_Callback(hObject, eventdata, handles)
Search=handles.Search;

if isfield(handles,'SIndex')
File=handles.File;
SIndex=handles.SIndex;
xListCandidates(File(SIndex),Search,Inf);
%winopen(OUT) %%OUT is not outputted by the above function
%winopen('Query_51_Results.txt')
else %If data is loaded from saved search results
    [File,SIndex]=zAddNTData(Search.Filenames,0);
    handles.File=File;
    handles.SIndex=SIndex;
    guidata(hObject, handles);
    
    xListCandidates(File(SIndex),Search,Inf);   %%%%%ALERT: SIndex here is not correct?
end
%xGroupCandidates(File(SIndex),Search);  % doesn't work very well yet!



% --- Executes on selection change in Overlap.
function Overlap_Callback(hObject, eventdata, handles)
% hObject    handle to Overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Overlap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Overlap


% --- Executes during object creation, after setting all properties.
function Overlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function GuarCutoff_Callback(hObject, eventdata, handles)
% hObject    handle to GuarCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GuarCutoff as text
%        str2double(get(hObject,'String')) returns contents of GuarCutoff as a double


% --- Executes during object creation, after setting all properties.
function GuarCutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GuarCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%</div>



function RelCutoff_Callback(hObject, eventdata, handles)
% hObject    handle to RelCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RelCutoff as text
%        str2double(get(hObject,'String')) returns contents of RelCutoff as a double


% --- Executes during object creation, after setting all properties.
function RelCutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RelCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in InteractionText.
function InteractionText_Callback(hObject, eventdata, handles)
xDescribeInteraction

% --- Executes on button press in MaxDistText.
function MaxDistText_Callback(hObject, eventdata, handles)
xDescribeDistance

% --- Executes on button press in NTmaskText.
function NTmaskText_Callback(hObject, eventdata, handles)
xDescribeMask

