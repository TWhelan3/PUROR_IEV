function varargout = ISMRMRDViewer(varargin)
% ISMRMRDVIEWER MATLAB code for ISMRMRDViewer.fig
%      ISMRMRDVIEWER, by itself, creates a new ISMRMRDVIEWER or raises the existing
%      singleton*.
%
%      H = ISMRMRDVIEWER returns the handle to a new ISMRMRDVIEWER or the handle to
%      the existing singleton*.
%
%      ISMRMRDVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ISMRMRDVIEWER.M with the given input arguments.
%
%      ISMRMRDVIEWER('Property','Value',...) creates a new ISMRMRDVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ISMRMRDViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ISMRMRDViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ISMRMRDViewer

% Last Modified by GUIDE v2.5 22-Oct-2015 16:58:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ISMRMRDViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @ISMRMRDViewer_OutputFcn, ...
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


% --- Executes just before ISMRMRDViewer is made visible.
function ISMRMRDViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ISMRMRDViewer (see VARARGIN)

% Choose default command line output for ISMRMRDViewer
handles.output = hObject;
handles.channel = 1;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ISMRMRDViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ISMRMRDViewer_OutputFcn(hObject, eventdata, handles) 
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
handles.openFile = uigetfile('*.ismrmrd');
set(handles.text3, 'String', handles.openFile);
set(handles.popupmenu1, 'Value', 1);
set(handles.popupmenu1, 'String', '/');
guidata(hObject, handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% looking for 'data'
%probably won't be this, but will be something standard

info=h5info(handles.openFile,[handles.toRead '/data']);
sz=info.Dataspace.Size;

complex_struct = h5read(handles.openFile, [handles.toRead '/data'],[1 1 1 handles.channel 1],[sz(1) sz(2) sz(3) 1 sz(end)]);

sliderStep = [1, 1] / (sz(end));
set(handles.slider1, 'Min', 1);
set(handles.slider1, 'Max', sz(end));
set(handles.slider1, 'SliderStep', sliderStep);

if(isstruct(complex_struct))
handles.complex_data = complex(complex_struct.real, complex_struct.imag);
handles.phase=angle(handles.complex_data);
else
handles.phase=complex_struct;
end
hold on
imagesc(handles.phase(:,:,1));
guidata(hObject,handles);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.slice=ceil(get(hObject, 'Value'));

hold on
imagesc(handles.phase(:,:,handles.slice));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'SliderStep',[1,1]);
set(hObject, 'Value', 1)
handles.slice=1;

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.channel=str2double(get(hObject,'String'));
guidata(hObject,handles);
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


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

groupcell=cellstr(get(hObject, 'String'));
groupname=groupcell{get(hObject, 'Value')};

info=h5info(handles.openFile, groupname);
handles.toRead=groupname;
if ~isempty(info.Groups)
    numGroups=size(info.Groups);
    for i=1:numGroups
        if ~ismember(info.Groups(i).Name, groupcell)
            groupcell = [groupcell ; cellstr(info.Groups(i).Name)];
        end
    end
end
% if ~isempty(info.Datasets) 
%     numSets=size(info.Datasets);
%     for i=1:numSets
%          if ~ismember(info.Datasets(i).Name, groupcell)
%             groupcell = [groupcell ; cellstr(info.Datasets(i).Name)];
%          end
%     end
% end
set(hObject, 'String', groupcell);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function combinedMag = getSOSMagnitude(filename, toRead)
info=h5info(filename, toRead);

sz=info.Dataspace.Size;
combinedMag=[];

stride=sz(5)/8;

combinedMag=zeros(sz(1),sz(2),sz(5));

for i=1:stride:sz(end)
    base=h5read(filename, toRead,[1 1 1 1 i], [sz(1) sz(2) sz(3) sz(4) stride]);

    if isstruct(base)
        base = complex(base.real, base.imag);
    end

    base=squeeze(base);
    base=abs(base).*abs(base);
    basesum=sum(base,3);
    basesum=squeeze(basesum);
    basesum=abs(basesum);
    combinedMag(:,:,i:(i+stride-1))=sqrt(basesum);
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.phase=getSOSMagnitude(handles.openFile, [handles.toRead '/data']);
sliderStep = [1, 1] / (size(handles.phase,3));
set(handles.slider1, 'Max', size(handles.phase,3));
set(handles.slider1, 'SliderStep', sliderStep);
imagesc(handles.phase(:,:,1)*10000); %this isn't scaling properly
guidata(hObject,handles);
