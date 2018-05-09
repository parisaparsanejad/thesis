function varargout = Zframe(varargin)
% ZFRAME MATLAB code for Zframe.fig
%      ZFRAME, by itself, creates a new ZFRAME or raises the existing
%      singleton*.
%
%      H = ZFRAME returns the handle to a new ZFRAME or the handle to
%      the existing singleton*.
%
%      ZFRAME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ZFRAME.M with the given input arguments.
%
%      ZFRAME('Property','Value',...) creates a new ZFRAME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Zframe_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Zframe_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Zframe

% Last Modified by GUIDE v2.5 07-May-2018 15:38:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Zframe_OpeningFcn, ...
                   'gui_OutputFcn',  @Zframe_OutputFcn, ...
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


% --- Executes just before Zframe is made visible.
function Zframe_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Zframe (see VARARGIN)

% Choose default command line output for Zframe
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Zframe wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Zframe_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% here we load the image
% global im im2 i
% [path, user_cancel]=imgetfile();
% if user_cancel
%     msgbox(sprintf('Error'),'Error','Error')
%     return
% end
% im=imread(path);
% im=im2double(im); % convert to double
% im2=im;  % for backup process
% axes(handles.axes1);
% imshow(im);



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% here we reset the image to original
global im im3
im=im3;
axes(handles.axes1);
imshow(im);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% here we extract the Z-marker from the background
global im im3
axes(handles.axes1);
im=imcrop(im);
im3=im;
axes(handles.axes1);
imshow(im)

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% here we do filtering
global im
im = imguidedfilter(im);
im = imadjust(im,[0.5 .6],[0 1]);
im = im2bw(im, 0.5);
axes(handles.axes1);
imshow(im)




% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% here we implement the masking procedure
global im 
im = im2double(im); 
im=im./(max(max(im)));
sliderVal=get(hObject,'value')
im = imgaussfilt(im,sliderVal);
axes(handles.axes1);
imshow(im)
% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im centroids im5
s = regionprops(im,'centroid');
centroids = cat(1, s.Centroid);
centroids = round(centroids);
axes(handles.axes1);
imshow(im);
hold on
plot(centroids(:,1),centroids(:,2), 'b*')
[m n]=size(centroids)
if m~=7
    msgbox(sprintf('please repeat the pre-processing step'),'Error','Error')
    return
end
im5=im;

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% here we do automatic detection of the Z-frame with neural network

global im im3

set(handles.figure1, 'pointer', 'watch')
drawnow;


% loading the detector only one time should be done
data=load('detectoroflargefasterCNN');
detector=data.detector;
[bboxes, scores] = detect(detector, im);
im = insertObjectAnnotation(im, 'rectangle', bboxes, scores);
axes(handles.axes1);
imshow(im);
set(handles.figure1, 'pointer', 'arrow')
% crop the detected Zframe area from the image
[limit2,m]=size(bboxes)
for i=1:limit2
    if bboxes(i,:)==0
       bboxes(i,:)=[ ]
    end 
   [limit2,m]=size(bboxes) 
   if limit2==1
       break
   end
end

im= imcrop(im,bboxes);
im3=im;
axes(handles.axes1);
imshow(im);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% here we trigger the webcam


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% here is the calibration process which containg transfering points from pixel to mm and ordering the points as opposit U shape
global im centroids P


[m n]=size(centroids)
if m~=7
    n=0
end
% ordering the 7 points
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
if abs(U1(1)-U3(1))<5 & U1(2)>U3(2)
   U3=centroids(H(1),:)
U1=centroids(H(2),:)

end

u=[U1;U3;U5;U7]

% [~,idx] = sort(u(:,1)); % sort just the first column
% u = u(idx,:)   % sort the whole matrix using the sort indices

% determie if U1 has two neighbors so two shapes of z-frame are possible

i=1
y=[0 0;0 0]
while y(1,:)==0 & i<7
    
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
i=i+1
end

i=1

while y(2,:)==0 & i<7
   
if centroids(i,1)>=U1(1,1) & centroids(i,1)<=U3(1,1)& centroids(i,2)>U1(1,2)& centroids(i,2)<U3(1,2)
   y(2,:)=centroids(i,:)
elseif centroids(i,1)>=U1(1,1) && centroids(i,1)<=U3(1,1)& centroids(i,2)>U3(1,2)& centroids(i,2)<U1(1,2)
  y(2,:)=centroids(i,:)
elseif centroids(i,1)<=U1(1,1) & centroids(i,1)>=U3(1,1)& centroids(i,2)>U1(1,2)& centroids(i,2)<U3(1,2)
  y(2,:)=centroids(i,:)
elseif centroids(i,1)<=U1(1,1) & centroids(i,1)>=U3(1,1)& centroids(i,2)>U3(1,2)& centroids(i,2)<U1(1,2)
   y(2,:)=centroids(i,:)
else 
     y(2,:)=[0 0]
end
 i=i+1
end

if (y(1,:)~=0 & y(2,:)~=0)
    t=1
else
    t=0
end  
% define a  degree counter-clockwise rotation matrix
% estimate the shape of the Z-frame if it is left_down or rleft_up
% find the rotation cloclwise or not
d35=U3-U5;
normd35=d35/norm(d35) 
if normd35(2)<0
    % U3 is down and U5 is up
    downpoint=U3;
    if t==1
        innerpoint=centroids(find(centroids(:,1)>=U3(1) & centroids(:,1)<=U7(1) & centroids(:,2)>U3(2) & centroids(:,2)<U7(2)))
        if isempty(innerpoint)
        innerpoint=centroids(find(centroids(:,1)<=U3(1) & centroids(:,1)>=U7(1) & centroids(:,2)>U3(2) & centroids(:,2)<U7(2)))
        end
        elseif t==0
        innerpoint=centroids(find(centroids(:,1)>U1(1) & centroids(:,1)<U3(1) & centroids(:,2)>=U3(2) & centroids(:,2)<=U1(2)))
        if isempty(innerpoint)
        innerpoint=centroids(find(centroids(:,1)>U1(1) & centroids(:,1)<U3(1) & centroids(:,2)<=U3(2) & centroids(:,2)>=U1(2)))
        end
        end
else
    % U5 is down
      downpoint=U5;
    if t==1
   
      innerpoint=centroids(centroids(:,1)>=U5(1) & centroids(:,1)<=U7(1) & centroids(:,2)>U5(2) & centroids(:,2)<U7(2))
    if isempty(innerpoint)
      innerpoint=centroids(find(centroids(:,1)<=U5(1) & centroids(:,1)>=U7(1) & centroids(:,2)>U5(2) & centroids(:,2)<U7(2)))
    end
    elseif t==0
          innerpoint=centroids(find(centroids(:,1)>U1(1) & centroids(:,1)<U5(1) & centroids(:,2)>U5(2) & centroids(:,2)<U1(2))) 
           if isempty(innerpoint)
          innerpoint=centroids(find(centroids(:,1)>U1(1) & centroids(:,1)<U5(1) & centroids(:,2)<=U5(2) & centroids(:,2)>=U1(2)))     
           end
     end
 end

% angle calculation
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

orianted_centroids=centroids-center
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
orianted_centroids =orianted_centroids* R           % apply the rotation about the origin
orianted_centroids=orianted_centroids+center
% labels = {'P1','P2','P3','P4','P5','P6','P7'};
% plot(orianted_centroids(:,1),orianted_centroids(:,2), 'b*')
% text(orianted_centroids(:,1),orianted_centroids(:,2),labels,'VerticalAlignment','top','HorizontalAlignment','left')

 
 
P(1:7,1:2) = 0;
    centsort = sortrows(orianted_centroids) %Sort by first row (x-position)
    P(1:3,1:2) = sortrows(centsort(1:3,1:2),2); %First three points are most left
    P(4,1:2) = centsort(4,1:2); %Fourth point is in the middle (along the x-axis)
    P(5:7,1:2) = flipud(sortrows(centsort(5:7,1:2),2)); %First three points are most left (*right)
    
    
labels = {'P1','P2','P3','P4','P5','P6','P7'};
plot(P(:,1),P(:,2), 'b*')
text(P(:,1),P(:,2),labels,'VerticalAlignment','top','HorizontalAlignment','left')

    
% line 1 with centers 1 and 7
 
 x1= [P(1) P(5)]
 y1= [P(1,2) P(5,2)]
 
 % line 2 with centers 3 and 6
 
 x2=[P(3) P(7)]
 y2=[P(3,2) P(7,2)]
 
%fit linear polynomial 

 p1 = polyfit(x1,y1,1);
 p2 = polyfit(x2,y2,1);
 
%calculate intersection

x_intersect = fzero(@(x) polyval(p1-p2,x),3)
y_intersect = polyval(p1,x_intersect)
% figure
% imshow(im)
hold on
line(x1,y1);
line(x2,y2);
plot(x_intersect,y_intersect,'r*')
hold off

% transfer the points into new coordinate system

P(:,1)=P(:,1)-x_intersect
P(:,2)=P(:,2)-y_intersect


 % Convert Pixel to mm
 [ps_camera_x,ps_camera_y] = camera_calibration(P) 
 P(:,1)=P(:,1)*ps_camera_x
 P(:,2)=P(:,2)*ps_camera_y

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  im P   
                                                                                             
  % calculation of fi which is a function that measures the fraction along the length of the frame 
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

% transformation matrix
start_points = pf'
goal_points = [P(2,:) 0;P(4,:) 0;P(6,:) 0] 
[d,start_points_transposed,tr] = procrustes(goal_points, start_points, 'scaling', false, 'reflection', false)
start_points_transposed = start_points * tr.T + tr.c
% transformation matrix end

%3D graph for nonuniform_data  
% https://blogs.mathworks.com/pick/2007/11/02/advanced-matlab-surface-plot-of-nonuniform-data/
x=a;
y=b;
z=c;
xlin=linspace(min(x),max(x),150);  % 150 shows the condence of the surface plot
ylin=linspace(min(y),max(y),150);
[X,Y]=meshgrid(xlin,ylin);
Z=griddata(x,y,z,X,Y,'cubic');
axes(handles.axes1)
surface=mesh(X,Y,Z); % interpolated
hold on
axes(handles.axes1);
labels = {'P2','P4','P4'};
plot3(x,y,z,'.','MarkerSize',15)  %marker size is the size of the points
text(x,y,z,labels,'VerticalAlignment','bottom','HorizontalAlignment','left')
xlabel('x')
ylabel('y')
zlabel('z')

hold on

% to calculate the center of the triangle
wx=mean(x)
wy=mean(y)
wz=mean(z)
plot3(wx,wy,wz,'*','MarkerSize',15)  %marker size is the size of the points


%calculate normal vector of the points
% https://stackoverflow.com/questions/17950002/calculating-the-normal-from-3-points
P0=pf(1:3,1)
P1=pf(1:3,2);
P2=pf(1:3,3);
normal = cross(P0-P1, P0-P2);
normal = normal / norm( normal ); % just to make it unit length
% figure 
% hold on 
axes(handles.axes1);
quiver3(wx, wy, wz*1.5, normal(1), normal(2), normal(3)*1.5)

% to confirm a vector is normal to the plane by confirming that the dot product is zero

disp(dot(P0 - P1,normal));
disp(dot(P0 - P2,normal));


% make a surface in yz plane (x=0 , rotation around x axis).
x_x=[0 0 0];
y_x=b;
z_x=c;
zlin_x=linspace(min(z_x),max(z_x),150);
ylin_x=linspace(min(y_x),max(y_x),150);
[Z_x,Y_x]=meshgrid(zlin_x,ylin_x);
X_x=griddata(z_x,y_x,x_x,Z_x,Y_x,'cubic');
% mesh(X_x,Y_x,Z_x); % interpolated
% axes(handles.axes1);
% hold on
% plot3(x_x,y_x,z_x,'.','MarkerSize',15)
% hold off

%calculate the normal vector in the yz plane
pf_x=[x_x;y_x;z_x];
P0_x=pf_x(1:3,1)
P1_x=pf_x(1:3,2);
P2_x=pf_x(1:3,3);
normal_x = cross(P0_x-P1_x, P0_x-P2_x);
normal_x = normal_x / norm( normal_x); % just to make it unit length
% quiver3(P0_x(1), P0_x(2), P0_x(3), normal_x(1), normal_x(2), normal_x(3))
% hold off

%calculate the angle between two normal vector  ( calculation of the angle
%between the marker plane and different axes
%The result of the atan2d is expressed in degrees

ax = atan2d(norm(cross(normal,normal_x)),dot(normal,normal_x))

% make a surface in xz plane (y=0 , rotation around y axis).
x_y=a;
y_y=[0 0 0];
z_y=c;
zlin_y=linspace(min(z_y),max(z_y),33);
xlin_y=linspace(min(x_y),max(x_y),33);
[Z_y,X_y]=meshgrid(zlin_y,xlin_y);
Y_y=griddata(z_y,x_y,y_y,Z_y,X_y,'cubic');
% mesh(X_y,Y_y,Z_y); % interpolated
% axes(handles.axes1);
% hold on
% plot3(x_y,y_y,z_y,'.','MarkerSize',15)
% hold off

%calculate the normal vector in the plane
pf_y=[x_y;y_y;z_y];
P0_y=pf_y(1:3,1)
P1_y=pf_y(1:3,2);
P2_y=pf_y(1:3,3);
normal_y = cross(P0_y-P1_y, P0_y-P2_y);
normal_y = normal_y / norm( normal_y); % just to make it unit length
% v= quiver3(P0_y(1), P0_y(2), P0_y(3), normal_y(1), normal_y(2), normal_y(3))
% hold off

%calculate the angle between two normal vector  ( calculation of the angle
%between the marker plane and different axes
%The result of the atan2d is expressed in degrees
ay = atan2d(norm(cross(normal,normal_y)),dot(normal,normal_y))

% make a surface in xy plane (z=0 , rotation around z axis).
x_z=a;
y_z=b;
z_z=[0 0 0];
xlin_z=linspace(min(x_z),max(x_z),33);
ylin_z=linspace(min(y_z),max(y_z),33);
[X_z,Y_z]=meshgrid(xlin_z,ylin_z);
Z_z=griddata(x_z,y_z,z_z,X_z,Y_z,'cubic');
% mesh(X_z,Y_z,Z_z); % interpolated
% axes(handles.axes1);
% hold on
% plot3(x_z,y_z,z_z,'.','MarkerSize',15)
hold off

%calculate the normal vector in the plane
pf_z=[x_z;y_z;z_z];
P0_z=pf_z(1:3,1);
P1_z=pf_z(1:3,2);
P2_z=pf_z(1:3,3);
normal_z = cross(P0_z-P1_z, P0_z-P2_z);
normal_z = normal_z / norm( normal_z ); % just to make it unit length
% v= quiver3(P0_z(1), P0_z(2), P0_z(3), normal_z(1), normal_z(2), normal_z(3))
% hold off

%calculate the angle between two normal vector  ( calculation of the angle
%between the marker plane and different axes
%The result of the atan2d is expressed in degrees
az = atan2d(norm(cross(normal,normal_z)),dot(normal,normal_z))
rotation_around_xaxis=90-ax
rotation_around_yaxis=90-ay
% fh = figure;
% prompt = {1 3};
% eh = uicontrol('Style','edit','String',num2str(90-ax));
set (handles.text3,'string',num2str(90-ax))
set (handles.text4,'string',num2str(90-ay))
% caption = sprintf('rotationx',num2str(90-ay));
caption1='rotation around x axis:'
set(handles.text5,'string',caption1);
caption2='rotation around y axis:'
set(handles.text6,'string',caption2);


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)

% here we adjust the contrast
global im 
cont=(get(hObject,'value'))
im = imadjust(im,[cont+0.2 cont+0.5],[0 1]);
axes(handles.axes1);
imshow(im)

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)


% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% here we do thresholding (binary image)
global im
binary=(get(hObject,'value'))
im= im2bw(im, binary);
axes(handles.axes1);
imshow(im)


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im im2
im=im2;

axes(handles.axes1);
imshow(im);

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% Here we do Autmatic image-preprocessing
% filtering
global im im4
im = imguidedfilter(im);
im = imadjust(im,[0.5 .6],[0 1]);
im = im2bw(im, 0.5);

im = bwareaopen(im,30);
im = imclearborder(im);
im4=im;
axes(handles.axes1);
imshow(im)

% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% here we do manual center detection
global im centerha centroids im5
size(im)
zoom(1.1);
zoomedImage = getframe(); % Get zoomed portion that is visible.
% size(zoomedImage.cdata)
axes(handles.axes1);
% imshow(zoomedImage.cdata);
s = regionprops(im,'centroid');
centroids = cat(1, s.Centroid);
centerha = round(centroids);
% axes(handles.axes1);
% imshow(im);
hold on

% plot(centroids(:,1),centroids(:,2), 'b*')

% Create draggable point
h1 = impoint(gca,centerha(1,:));
% Update position in title using newPositionCallback
addNewPositionCallback(h1,@(h1) title(sprintf('(%1.0f,%1.0f)',h1(1),h1(2))));
% Interactively place a point. Use wait to block the MATLAB command line. 
centroids(1,:)= wait(h1);
h2 = impoint(gca,centerha(2,:));
addNewPositionCallback(h2,@(h2) title(sprintf('(%1.0f,%1.0f)',h2(1),h2(2))));
centroids(2,:) = wait(h2);
h3 = impoint(gca,centerha(3,:));
addNewPositionCallback(h3,@(h3) title(sprintf('(%1.0f,%1.0f)',h3(1),h3(2))));
centroids(3,:) = wait(h3);
h4 = impoint(gca,centerha(4,:));
addNewPositionCallback(h4,@(h4) title(sprintf('(%1.0f,%1.0f)',h4(1),h4(2))));
centroids(4,:) = wait(h4);
h5 = impoint(gca,centerha(5,:));
addNewPositionCallback(h5,@(h5) title(sprintf('(%1.0f,%1.0f)',h5(1),h5(2))));
centroids(5,:)= wait(h5);
h6 = impoint(gca,centerha(7,:));
addNewPositionCallback(h6,@(h6) title(sprintf('(%1.0f,%1.0f)',h6(1),h6(2))));
centroids(7,:) = wait(h6);
h7 = impoint(gca,centerha(6,:));
addNewPositionCallback(h7,@(h7) title(sprintf('(%1.0f,%1.0f)',h7(1),h7(2))));
centroids(6,:)= wait(h7);
uiwait(msgbox('Center points selection completed'));
im5=im;
axes(handles.axes1);
imshow(im);
hold on
plot(centroids(:,1),centroids(:,2), 'r*')


% button=0;
% zlvl=1;
% xl = get(gca,'xlim');
% xlen = size(im,2);
% yl = get(gca,'ylim');
% ylen = size(im,1);
% 
% while button~=2
%     % Get the mouse position on the axes (needed for binary image editing) and button number
%     [y,x,button]=ginput(1);
% 
%     % Determine if it is a zoom-in or zoom-out
%     if button==1
%         zlvl = zlvl*2;
%         zoom(2);
%     elseif button==3
%         zlvl = zlvl/2;
%         if zlvl<1, zlvl=1; end % No zoom level smaller than 1
%         zoom(0.5);
%     end
% 
%     % Change the axis limits to where the mouse click has occurred
%     % and make sure that the display window is within the image dimensions
%     xlimit = [x-xlen/zlvl/2+0.5 x+xlen/zlvl/2+0.5];
%     if xlimit(1)<0.5, xlimit=[0.5 xlen/zlvl+0.5]; end
%     if xlimit(2)>0.5+xlen, xlimit=[xlen-xlen/zlvl+0.5 xlen+0.5]; end
%     xlim(xlimit);
% 
%     ylimit = [y-ylen/zlvl/2+0.5 y+ylen/zlvl/2+0.5];
%     if ylimit(1)<=0.5, ylimit=[0.5 ylen/zlvl+0.5]; end
%     if ylimit(2)>=0.5+ylen, ylimit=[ylen-ylen/zlvl+0.5 ylen+0.5]; end
%     ylim(ylimit);
% end
% 
% 
% 
% uiwait(msgbox('Locate the points'));
% centroids_input(1:7,1:2) = 0;
% centroids_input(1,:)=ginput(1);
% hold on; % Prevent image from being blown away.
% plot(centroids_input(1,1),centroids_input(1,2),'r+', 'MarkerSize', 10);



% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

%get value as a string from text box and convert to double
global ps 
ps = str2double(get(handles.edit1,'string'));
%check that value is a number and display error if it is not
if isnan(ps) || ~isreal(ps)
  errordlg('Unvalid Pixel Spacing Value','Try again','modal') 
    return
end
% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
global d
d = str2double(get(handles.edit2,'string'));
%check that value is a number and display error if it is not
if isnan(d) || ~isreal(d)
  errordlg('Unvalid Value for Size of Marker','Try again','modal') 
    return
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% here we Remove small objects which don't blong to the Zframe
global im im4
mask=round(get(hObject,'value'));
im = bwareaopen(im,mask);
im= imclearborder(im);
im4=im;
axes(handles.axes1);
imshow(im)


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im im4
im=im4;
axes(handles.axes1);
imshow(im);

% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im im5 centroids
im=im5;
axes(handles.axes1);
imshow(im);
hold on
plot(centroids(:,1),centroids(:,2), 'b*')
hold off
[m n]=size(centroids)
if m~=7
    msgbox(sprintf('please repeat the pre-processing step'),'Error','Error')
    return
end


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
% here we choose the input image
global im im2 I I2
contents=cellstr(get(hObject,'String'));
popChoice=contents{get(hObject,'Value')};
if(strcmp(popChoice,'Import image from folder'))
[path, user_cancel]=imgetfile();
if user_cancel
    msgbox(sprintf('Please Select an Image'),'Error','Error')
    return
end
im=imread(path);
im=im2double(im); % convert to double
im2=im;  % for backup process
axes(handles.axes1);
imshow(im);    
    
end
% reading the CORONAL image for needle pathway through target
I=imread('CORONAL.jpg');
% i=im2double(i); % convert to double
I2=I;
axes(handles.axes3);
imshow(I);
% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
% here we choose if we want to enter the value of Pixel spacing
global button_ps
if (get(hObject,'Value')==get(hObject,'Max'))
    button_ps=1;
else 
    button_ps=0;
end





% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% here we use mouse to correct circle boundary
global im im4
size(im)
zoom(1.3);
zoomedImage = getframe(); % Get zoomed portion that is visible.
size(zoomedImage.cdata)
axes(handles.axes1)
% figure; % Bring up separate figure.
imshow(zoomedImage.cdata);
hFH = imfreehand();
% Get the xy coordinates of where they drew.
xy = hFH.getPosition
% get rid of imfreehand remnant.
delete(hFH);
% Overlay what they drew onto the image.
hold on; % Keep image, and direction of y axis.
xCoordinates = xy(:, 1);
yCoordinates = xy(:, 2);
plot(xCoordinates, yCoordinates,'white', 'MarkerSize', 20);
caption = sprintf('Original Grayscale Image.\nPoints may not lie on adjacent pixels, depends on your speed of drawing!');
% title(caption, 'FontSize', fontSize);
 
% Ask user if they want to burn the line into the image.
promptMessage = sprintf('Do you want to burn the line into the image?');
titleBarCaption = 'Continue?';
button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
if strcmpi(button, 'Yes')
	cla;
	hold off;
	for k = 1 : length(xCoordinates)
		row = int32(yCoordinates(k));
		column = int32(xCoordinates(k));
		im(row, column) = 255;
	end
	imshow(im, []);
	axis on;
	caption = sprintf('Grayscale Image with Burned In Curve.\nPoints may not lie on adjacent pixels, depends on your speed of drawing!');
% 	title(caption, 'FontSize', fontSize);
end

% Ask user if they want to interpolate the line to get the "in-between" points that are missed..
% promptMessage = sprintf('Do you want to interpolate the curve into intervening pixels?');
% titleBarCaption = 'Continue?';
% button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
% if strcmpi(button, 'Cancel')
% 	return;
% end
xCoordinates = xy(:, 1);
yCoordinates = xy(:, 2);
numberOfKnots = length(xCoordinates);

% Close gaps that you get when you draw too fast.
% Use splines to interpolate a smoother curve,
% with 10 times as many points,
% that goes exactly through the same data points.
samplingRateIncrease = 5;
newXSamplePoints = linspace(1, numberOfKnots, numberOfKnots * samplingRateIncrease);
% smoothedY = spline(xCoordinates, yCoordinates, newXSamplePoints);
% Make the 2D array where the top row is the x coordinates and the bottom row is the y coordinates,
% but with the exception that the left column and right column is a vector that gives the direction of the slope.
yy = [0, xCoordinates', 0; 1, yCoordinates', 1]
pp = spline(1:numberOfKnots, yy); % Get interpolant
smoothedY = ppval(pp, newXSamplePoints); % Get smoothed y values in the "gaps".
% smoothedY is a 2D array with the x coordinates in the top row and the y coordinates in the bottom row.
smoothedXCoordinates = smoothedY(1, :);
smoothedYCoordinates = smoothedY(2, :);
% Plot smoothedY and show how the line is
% smooth, and has no sharp bends.
hold on; % Don't destroy the first curve we plotted.
hGreenCurve = plot(smoothedXCoordinates, smoothedYCoordinates, '-g');
% title('Spline Interpolation Demo', 'FontSize', 8);
% But smoothedXCoordinates and smoothedYCoordinates are not in pixel coordinates, they have fractional values.
% If you want integer pixel values, you have to round.
intSmoothedXCoordinates = int32(smoothedXCoordinates);
intSmoothedYCoordinates = int32(smoothedYCoordinates);
% But now it's possible that some coordinates will be on the same pixel if that's
% how they rounded according to how they were located to the nearest integer pixel location.
% So use diff() to remove elements that have the same x and y values.
diffX = [1, diff(intSmoothedXCoordinates)];
diffY = [1, diff(intSmoothedYCoordinates)];
% Find out where both have zero difference from the prior point.
bothZero = (diffX==0) & (diffY == 0);
% Remove those from the arrays.
finalX = intSmoothedXCoordinates(~bothZero);
finalY = intSmoothedYCoordinates(~bothZero);
% Now remove the green line.
delete(hGreenCurve);
% Plot the final coordinates.
hGreenCurve = plot(finalX, finalY, '-y');

% bw=bwareaopen(im,10)
% bw = imfill(bw,'holes');

m = imfill(im,'holes');
im4=im;
axes(handles.axes1)
imshow(im)


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% here we define the center of the middle point and then calculate the
% trajectory 
global I I2
% Filtering Steps
%Use double precision and scale image values to the range of [0 1]
% axes(handles.axes3);
% imshow(I);
% 
% I=imcrop(I);
% axes(handles.axes3);
% imshow(I)
% I = im2double(I); 
% I=I./(max(max(I)));
% 
% %2D Gaussian filter
% I=imgaussfilt(I,2);
% figure;imshow(I);title('After Gaussian filtering')
% 
% % Adjust image contrast (stretch values)
% I= imadjust(I,[0.4 0.9],[]);
% h=figure('Position',[100 100 1000 1000]);hold on; %Normal plot
% imshow(I);title('After image adjustment')
% 
% %Thresholding (create binary image)
% I= im2bw(I, 0.25);
% figure;imshow(I);title('After threshold (binary)')
axes(handles.axes3);
imshow(I);
[x,y] = ginput(1)
axes(handles.axes3);
imshow(I);
hold on
plot(x,y, 'b*')
hold off


A = [x 500]; 
B = [y y]; 

% figure
% imshow(IBi)
axes(handles.axes3);
imshow(I);
hold on
plot(A,B,'*')
line(A,B)
hold off


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% here we choose the target and calculate the trajectry from the needle
% point to te trajectory
global I2
axes(handles.axes3);
imshow(I2);
[x,y] = ginput(2)
axes(handles.axes3);
imshow(I2);
hold on
plot(x,y, 'b*')
hold off
% make a line
A = [x(1) x(2)]; 
B = [y(2) y(2)]; 
axes(handles.axes3);
imshow(I2);
hold on
plot(A,B,'*')
p1 = polyfit(A,B,1);
%line(A,B)
hold off


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function Theta_1_Callback(hObject, eventdata, handles)
% hObject    handle to Theta_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Theta_1 as text
%        str2double(get(hObject,'String')) returns contents of Theta_1 as a double


% --- Executes during object creation, after setting all properties.
function Theta_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Theta_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Theta_2_Callback(hObject, eventdata, handles)
% hObject    handle to Theta_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Theta_2 as text
%        str2double(get(hObject,'String')) returns contents of Theta_2 as a double


% --- Executes during object creation, after setting all properties.
function Theta_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Theta_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Theta_3_Callback(hObject, eventdata, handles)
% hObject    handle to Theta_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Theta_3 as text
%        str2double(get(hObject,'String')) returns contents of Theta_3 as a double


% --- Executes during object creation, after setting all properties.
function Theta_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Theta_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_forward.
function btn_forward_Callback(hObject, eventdata, handles)
% hObject    handle to btn_forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Theta_1 Theta_2 Theta_3

Th_1=str2double(handles.Theta_1.String)*pi/180;
Th_2=str2double(handles.Theta_2.String)*pi/180;
Th_3=str2double(handles.Theta_3.String)*pi/180;

L_1=20;
L_2=50;
L_3=40;

L(1)=Link([0 L_1 0 pi/2]);
L(2)=Link([0 0 L_2 0]);
L(3)=Link([0 0 L_3 0]);

Robot=SerialLink(L);
Robot.name='RRR_Robot';
axes(handles.axes3);
Robot.plot([Th_1 Th_2 Th_3]);
T=Robot.fkine([Th_1 Th_2 Th_3]);
handles.Pos_X.String=num2str(floor(T(1,4)));
handles.Pos_Y.String=num2str(floor(T(2,4)));
handles.Pos_Z.String=num2str(floor(T(3,4)));




function Pos_X_Callback(hObject, eventdata, handles)
% hObject    handle to Pos_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pos_X as text
%        str2double(get(hObject,'String')) returns contents of Pos_X as a double


% --- Executes during object creation, after setting all properties.
function Pos_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pos_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pos_Y_Callback(hObject, eventdata, handles)
% hObject    handle to Pos_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pos_Y as text
%        str2double(get(hObject,'String')) returns contents of Pos_Y as a double


% --- Executes during object creation, after setting all properties.
function Pos_Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pos_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pos_Z_Callback(hObject, eventdata, handles)
% hObject    handle to Pos_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pos_Z as text
%        str2double(get(hObject,'String')) returns contents of Pos_Z as a double


% --- Executes during object creation, after setting all properties.
function Pos_Z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pos_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_Inverse.
function btn_Inverse_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Inverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Theta_1 Theta_2 Theta_3
PX=str2double(handles.Pos_X.String);
PY=str2double(handles.Pos_Y.String);
PZ=str2double(handles.Pos_Z.String);

L_1=20;
L_2=50;
L_3=40;

L(1)=Link([0 L_1 0 pi/2]);
L(2)=Link([0 0 L_2 0]);
L(3)=Link([0 0 L_3 0]);

Robot=SerialLink(L);
Robot.name='RRR_Robot';

T=[1 0 0 PX;0 1 0 PY; 0 0 1 PZ;0 0 0 1];
J=Robot.ikine(T,[0 0 0],[1 1 1 0 0 0])*180/pi;
handles.Theta_1.String=num2str(floor(J(1)));
handles.Theta_2.String=num2str(floor(J(2)));
handles.Theta_3.String=num2str(floor(J(3)));
axes(handles.axes3);
Robot.plot(J*pi/180);
