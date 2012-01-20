
figure(2)



% a figure
clear sh;

figure(1)
clf
sh(1)=subplot(1,1,1);
peaks(16);

figure(2)
clf
sh(2)=subplot(1,1,1);
peaks(32);

% the engine
linkobj = linkprop(sh,...
{'cameraposition',...
 'cameraupvector',...
 'cameratarget',...
 'cameraviewangle'});

set(gcf, 'UserData', linkobj);

break

% a figure
clear sh;

sh(1)=figure(3);
peaks(16);

sh(2)=figure(4);
peaks(32);

% the engine
linkobj = linkprop(sh,...
{'cameraposition',...
 'cameraupvector',...
 'cameratarget',...
 'cameraviewangle'});

set(gcf, 'UserData', linkobj);

