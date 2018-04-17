function varargout = first_trial(varargin)
% FIRST_TRIAL MATLAB code for first_trial.fig
%      FIRST_TRIAL, by itself, creates a new FIRST_TRIAL or raises the existing
%      singleton*.
%
%      H = FIRST_TRIAL returns the handle to a new FIRST_TRIAL or the handle to
%      the existing singleton*.
%
%      FIRST_TRIAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIRST_TRIAL.M with the given input arguments.
%
%      FIRST_TRIAL('Property','Value',...) creates a new FIRST_TRIAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before first_trial_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to first_trial_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help first_trial

% Last Modified by GUIDE v2.5 17-Apr-2018 10:02:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @first_trial_OpeningFcn, ...
                   'gui_OutputFcn',  @first_trial_OutputFcn, ...
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


% --- Executes just before first_trial is made visible.
function first_trial_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to first_trial (see VARARGIN)

% Choose default command line output for first_trial
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes first_trial wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = first_trial_OutputFcn(hObject, eventdata, handles) 
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
% here we load image
global im im2
[path,user_cance]=imgetfile();
if user_cance
    msgbox(sprintf('Error'),'Error','Error')
    return
end
im=imread(path);
im=im2double(im); % convert to double
im2=im % for backup process
axes(handles.axes1);
imshow( im);
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% here we reset the image to original
global im2 
axes(handles.axes1)
imshow(im2);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% here is the black and white image
global im
imblack=im;
imblack=1-im;
axes(handles.axes1);
imshow(imblack);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% here we convert the mage to grayscale using rgb=(r+G+b)/3
global im 
imgray=(im(:,:,1)+im(:,:,2)+im(:,:,2))/3;
axes(handles.axes1);
imshow(imgray);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% here we implemet simple brightness
global im2
val=0.5* get(hObject,'value')-0.5;
imbright=im2+val;
axes(handles.axes1);
imshow(imbright);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
