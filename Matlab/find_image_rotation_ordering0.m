%% read rotated image and filtering it
close all;
clear all;
clc;

%% load image
[y,map]=imread('image12.jpg');
imshow(y,map);
I=imcrop(y,map);
imshow(I,map);

%% filtering steps
%Use double precision and scale image values to the range of [0 1]
y = im2double(I); 
y=y./(max(max(y)));

%2D Gaussian filter
y=imgaussfilt(y,2);
figure;imshow(y);title('After Gaussian filtering')

% Adjust image contrast (stretch values)
yy = imadjust(y,[0.5 0.6],[]);
h=figure('Position',[100 100 1000 1000]);hold on; %Normal plot
imshow(yy);title('After image adjustment')

%Thresholding (create binary image)
yy= im2bw(yy, 0.9);
figure;imshow(yy);title('After threshold (binary)')

% Extract only circle larger than 7
mask = bwareaopen(yy,150);

% Get rid of white frame around outside border.
mask= imclearborder(mask);
imshow(mask);
title('Binary Image Mask');



% Find centroids and plot them with the binary image

s = regionprops(mask,'centroid');
centroids = cat(1, s.Centroid);
centroids = round(centroids);
figure
imshow(mask)
hold on
plot(centroids(:,1),centroids(:,2), 'b*')
hold off



%% ordering the 7 points
% calculating the euclidean distance
for i=1:7;
    for j=1:7;
  d=centroids(i,:)-centroids(j,:)
  c(i,j)=norm(d)
end
end
n=max(c)
[row,col]=find(c==max(c))
% now peida kardane do ta az noghate max k yek san nabashand
m=sort(n,'descend') % sorted from max to min
m= unique(m)        % removal of duplicate value 
m=sort(m,'descend')  % sorted from max to min

[H1,O1] = find(c==m(1))
[H2,O2] = find(c==m(2))
H=[H1(1);H1(2);H2(1);H2(2)]
H=sort(H,'ascend')

U1=centroids(H(1),:)
U3=centroids(H(2),:)
U5=centroids(H(3),:)
U7=centroids(H(4),:)

if abs(U5(1)-U7(1))<5 & U5(2)>U7(2)
    U7=centroids(H(3),:)
    U5=centroids(H(4),:)
end


u=[U1;U3;U5;U7]

% [~,idx] = sort(u(:,1)); % sort just the first column
% u = u(idx,:)   % sort the whole matrix using the sort indices

% https://stackoverflow.com/questions/27070878/get-all-the-neighborhood-of-point-x-of-distance-r-using-matlab


% [U5]=centroids(find(centroids(:,1)>u(1,1) & centroids(:,1)<u(2,1) & centroids(:,2)>u(1,2) & centroids(:,2)<u(2,2)))
% [U6]=centroids(find(centroids(:,1)>u(1,1) & centroids(:,1)<u(3,1) & centroids(:,2)>u(1,2) & centroids(:,2)<u(3,2)))

%[U5]=centroids(find(centroids(2,1)>u(1,1) & centroids(2,1)<u(2,1) & centroids(2,2)>u(1,2) & centroids(2,2)<u(2,2)))

% determie if U1 has two neighbors so two shapes of z-frame are possible

i=1
y=[0 0;0 0]
while y(1,:)==0 & i<7
    i=i+1
if centroids(i,1)>U1(1,1) & centroids(i,1)<U5(1,1)& centroids(i,2)>=U1(1,2)& centroids(i,2)<=U5(1,2)
   y(1,:)=centroids(i,:)
elseif centroids(i,1)>U1(1,1) & centroids(i,1)<U5(1,1)& centroids(i,2)>=U5(1,2)& centroids(i,2)<=U1(1,2)
  y(1,:)=centroids(i,:)
elseif centroids(i,1)<U1(1,1) & centroids(i,1)>U5(1,1)& centroids(i,2)>=U1(1,2)& centroids(i,2)<=U5(1,2)
   y(1,:)=centroids(i,:)
elseif centroids(i,1)<U1(1,1) & centroids(i,1)>U5(1,1)& centroids(i,2)>=U5(1,2)& centroids(i,2)<=U1(1,2)
   y(1,:)=centroids(i,:)
else 
     y(1,:)=[0 0]
end
end

i=1

while y(2,:)==0 & i<7
    i=i+1
if centroids(i,1)>=U1(1,1) & centroids(i,1)<=U3(1,1)& centroids(i,2)>U1(1,2)& centroids(i,2)<U3(1,2)
   y(2,:)=centroids(i,:)
elseif centroids(i,1)>=U1(1,1) & centroids(i,1)<=U3(1,1)& centroids(i,2)>U3(1,2)& centroids(i,2)<U1(1,2)
  y(2,:)=centroids(i,:)
elseif centroids(i,1)<=U1(1,1) & centroids(i,1)>=U3(1,1)& centroids(i,2)>U1(1,2)& centroids(i,2)<U3(1,2)
   y(2,:)=centroids(i,:)
elseif centroids(i,1)<=U1(1,1) & centroids(i,1)>=U3(1,1)& centroids(i,2)>U3(1,2)& centroids(i,2)<U1(1,2)
   y(2,:)=centroids(i,:)
else 
     y(2,:)=[0 0]
end
end

if (y(1,:)~=0 & y(2,:)~=0)
    t=1
else
    t=0
end  
% define a  degree counter-clockwise rotation matrix
% estimate the shape of the Z-frame if it is left_down or rleft_up
% find the rotation cloclwise or not
d35=U3-U5
normd35=d35/norm(d35) 
if normd35(2)<0
    % U3 is down and U5 is up
    downpoint=U3;
    if t==1
        innerpoint=centroids(find(centroids(:,1)>U3(1) & centroids(:,1)<U7(1) & centroids(:,2)>U3(2) & centroids(:,2)<U7(2)))
    elseif t==0
               innerpoint=centroids(find(centroids(:,1)>U1(1) & centroids(:,1)<U3(1) & centroids(:,2)>U3(2) & centroids(:,2)<U1(2)))
    end
else
    % U5 is down
      downpoint=U5;
    if t==1
        innerpoint=centroids(find(centroids(:,1)>U5(1) & centroids(:,1)<U7(1) & centroids(:,2)>U5(2) & centroids(:,2)<U7(2))) 
    elseif t==0
          innerpoint=centroids(find(centroids(:,1)>U1(1) & centroids(:,1)<U5(1) & centroids(:,2)>U5(2) & centroids(:,2)<U1(2)))     
        end
 end

%% angle calculation
j=2.9
center=U1
U7_orianted=U7

% I calculate the angle between the lines U1 and the diagonal of the square
% ( line between U1 and U7)



if t==1 
    theta=abs(atan((downpoint(1)-U1(1))/(downpoint(2)-U1(2))))
    if innerpoint~=0
    theta=theta+1.5708
    end
    
end
if t== 0
    theta=abs(atan((downpoint(1)-U1(1))/(downpoint(2)-U1(2))))
    if innerpoint~=0
    theta=theta+3.14159
    end
  if isempty(innerpoint)
    theta=theta+4.71239
  end
 end



% 
% if t==1 
% theta=abs(atan((U7(1)-U1(1))/(U7(2)-U1(2))))-0.785398
% if innerpoint~=0
%     theta=theta+1.5708
% end
% elseif t== 0
% theta=3.14159-abs(atan((U7(1)-U1(1))/(U7(2)-U1(2))))-0.785398
% theta=-theta
% if innerpoint~=0
%     theta=-(abs(theta)+1.5708)
% end
% end

orianted_centroids=centroids-center
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
orianted_centroids =orianted_centroids* R           % apply the rotation about the origin
orianted_centroids=orianted_centroids+center

plot(orianted_centroids(:,1),orianted_centroids(:,2), 'b*')


% while ( abs(atan((U7_orianted(1)-U1(1))/(U7_orianted(2)-U1(2))))<0.75 | abs(atan((U7_orianted(1)-U1(1))/(U7_orianted(2)-U1(2))))>0.80)
% j=j+0.1
% theta = pi/j  % pi/j radians 
% 
% if innerpoint~=0
% theta=1.5708+theta+1.5708
% end
% 
% if t==0
%     theta=-theta
% end
% U7_orianted=U7;
% 
% orianted_centroids=centroids-center
% R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% orianted_centroids =orianted_centroids* R           % apply the rotation about the origin
% orianted_centroids=orianted_centroids+center
% 
% U7_orianted=U7_orianted-center
% U7_orianted=U7_orianted* R           % apply the rotation about the origin for the U7
% U7_orianted=U7_orianted+center
% end

    
    
% estimate the shape of the Z-frame 


%% matrix calculation
P(1:7,1:2) = 0;
    centsort = sortrows(orianted_centroids) %Sort by first row (x-position)
    P(1:3,1:2) = sortrows(centsort(1:3,1:2),2); %First three points are most left
    P(4,1:2) = centsort(4,1:2); %Fourth point is in the middle (along the x-axis)
    P(5:7,1:2) = flipud(sortrows(centsort(5:7,1:2),2)); %First three points are most left (*right)

f(1:6)=0;
for i=[2,4,6]
f_i= norm(P(i+1,:)-P(i,:))/ norm(P(i+1,:)-P(i-1,:))
f(i)=f_i;
end
% estimation of The three corresponding points in the frame coordinate system
lx=40;
ly=40;
lz=40;
a=lx * [-1/2 1/2-f(4) 1/2];
b=ly *[1/2-f(2) 1/2 -1/2+f(6)];
c=lz *[-1/2+f(2) -1/2+f(4) -1/2+f(6)];
pf=[a;b;c]

