%function test_zBestTransformation(option)

%X and Y are two 3D objects (molecules). %Each has a subgroup of atoms (A and B) that need to be superimposed

%on each other as best as possible. The object Y needs to be rotated and translated so that its points B fall as

%close as possible to points A from object X. Scaling of object Y is not allowed (since these are real molecules,

%so you do not want to make them bigger or smaller)


close all %close opened figures

if option ==1

X=[-.5,1,0; 0 0 0; 1,.75,0; 0.5,1,0];
Y=[-3,-1.5,0; -5.5,-.75,0; -7,-3.3,0; -5,-3,0];

Y = 1.1*X + ones(size(X));

A=X;
B=Y;

elseif option ==2

X=[-.5,1,0; 0,0,0; 1,0,0; 2,.75,0];

Y=[-4.5,-3,0; -5.5,-.75,0; -7,-0.75,0; -9,-2.3,0];

A=X(1:3,:,:); %Align only first 3 points from each object

B=Y(1:3,:,:);

elseif option ==3

X=[-.5,1,0; 0,0,0; .5,0,0; 0,-1,0; 1,-1,0; 1,0,0; 2,.75,0];

Y=[-4.5,-3,0; -5.5,-.75,0; -6,-.75,0; -7,1,0; -6,2,0; -6.5,1,0; -5.8,.5,0; -6.5,-.75,0; -7,-0.75,0; -9,-2.3,0];

A=[X(1:2,:,:); X(end-1:end,:,:)]; %Align first 2 points and last 2 points

B=[Y(1:2,:,:); Y(end-1:end,:,:)];

end



% % A=X(grp,:);

% % % mAngle(X(1,:), X(2,:), X(3,:))

% % B=Y(grp,:);

% % % mAngle(Y(1,:), Y(2,:), Y(3,:))



%Transformations:

%Which of these transformations works best?

[R,scale,shift,sshift] = zBestTransformation(B,A); %first variable moves onto second (Y's onto X's; X's stay stationary)

R
scale
shift
sshift

Y1 = (ones(length(Y(:,1)),1)*sshift' + scale*Y*R'); %This scales the Y object to fit X better, changing its initial size

Y2 = (ones(length(Y(:,1)),1)*sshift' + scale*Y*R')/scale; %This should screw up the shift, since now it is shifted by shift'/scale. For the two examples here it looks as if it working, but for other examples it does not work

Y3 = ones(length(Y(:,1)),1)*shift' + Y*R'; %This should work, but it doesn't. Shift seems to be wrong

%Y3 = (shift*ones(1,length(Y(:,1))) + R*Y')'; %best fit without scaling %%SAME AS PREVIOUS LINE

R = zBestTransformation(B,A);

Y4 = Y*R'; %This should rotate so that objects (A and B) are parallel, but it doesn't work

% CLZ: I modified the previous line to use R' instead of R

Y5 = Y4 - ones(length(Y(:,1)),1)*Y4(2,:); %This is rotation followed by translation




%Check how each transformation looks like:

%Plot before transformation:

figure(100)

set(gcf,'name','Original objects','Units','Normalized','Position',[.20 .20 .4 .4])

plot(X(:,1),X(:,2))

hold on

plot(A(:,1),A(:,2),'*') %These are the points that need to be superimposed as best as possible

plot(Y(:,1),Y(:,2),'r')

plot(B(:,1),B(:,2),'r*') % //

axis equal



%%Plot after transformation:

figure(1)

set(gcf,'name','Y1: transformed & scaled','Units','Normalized','Position',[.22 .22 .4 .4])

plot(X(:,1),X(:,2))

hold on

plot(Y1(:,1),Y1(:,2),'r')

axis equal


figure(2)

set(gcf,'name','Y2: transformed & scaled / scale','Units','Normalized','Position',[.24 .24 .4 .4])

plot(X(:,1),X(:,2))

hold on

plot(Y2(:,1),Y2(:,2),'r')

axis equal


figure(3)

set(gcf,'name','Y3: transformed & NOT scaled -- should work but does not! (shift is way off)','Units','Normalized','Position',[.26 .26 .4 .4])

plot(X(:,1),X(:,2))

hold on

plot(Y3(:,1),Y3(:,2),'r')

axis equal



figure(4)

set(gcf,'name','Y4: Rotated only','Units','Normalized','Position',[.28 .28 .4 .4])

plot(X(:,1),X(:,2))

hold on

plot(Y4(:,1),Y4(:,2),'r')

axis equal


figure(5)

set(gcf,'name','Y5: Rotated & translated','Units','Normalized','Position',[.30 .30 .4 .4])

plot(X(:,1),X(:,2))

hold on

plot(Y5(:,1),Y5(:,2),'r')

axis equal



