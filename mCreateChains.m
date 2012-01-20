% mCreateChains
% Create popupmenus (and delete extra ones) for determining Query.ChainPopup
% mDetermineChains 
% this is a program that extracts chain information from the already read 
% query PDB file

for i=1:12
    h=findobj('Tag',strcat('ChainPopup',num2str(i)));
    delete(h);
end

for i=1:length(NT)
  [a,b,Ch] = zIndexLookup(File(QIndex),NT(i));
  handles.ChainPopup(i) = uicontrol('Tag',strcat('ChainPopup',num2str(i)),'Style','popupmenu','Units','normalized','Position',[(0.25+0.057*i) (0.795) .054 .04],'String',Ch,'Background',[1 1 1]);
end
set(handles.QueryChains,'Visible','on')
                             % this is the text just to the left of the popups
%%%End




return

for i=1:length(NT)
%    Ch=mGetChainsOfNTnumber(File(QIndex),NT(i)); % So now Ch is individual for each numcleotide number
%     if length(Ch)>1
        handles.ChainPopup(i) = uicontrol('Tag',strcat('ChainPopup',num2str(i)),'Style','popupmenu','Units','normalized','Position',[(0.25+0.057*i) (0.795) .054 .04],'String',Ch,'Background',[1 1 1]);
%     else
%         handles.ChainPopup(i) = uicontrol('Tag',strcat('ChainPopup',num2str(i)),'Style','text','Units','normalized','Position',[(0.25+0.057*i) (0.795) .054 .04],'String',Ch,'Background',[1 1 1]);
%     end
end
set(handles.QueryChains,'Visible','on')%this the text just to the left of the popups
%%%End
