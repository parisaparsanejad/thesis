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

% Last Modified by GUIDE v2.5 17-Apr-2018 11:40:52

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
global im im2
[path, user_cancel]=imgetfile();
if user_cancel
    msgbox(sprintf('Error'),'Error','Error')
    return
end
im=imread(path);
im=im2double(im); % convert to double
im2=im;  % for backup process
axes(handles.axes1);
imshow(im);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% here we reset the image to original
global im2
axes(handles.axes1);
imshow(im2);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% here we extract the Z-marker from the background
global im 
im=imcrop(im);
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
im = bwareaopen(im,60);
im = imclearborder(im);
axes(handles.axes1);
imshow(im)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im
s = regionprops(im,'centroid');
centroids = cat(1, s.Centroid);
centroids = round(centroids);
axes(handles.axes1);
imshow(im);
hold on
plot(centroids(:,1),centroids(:,2), 'b*')
hold off
[m n]=size(centroids)
if m~=7
    msgbox(sprintf('please press the pre-processing bottom'),'Error','Error')
    return
end