function varargout = PVT(varargin)
% PVT MATLAB code for PVT.fig
%      PVT, by itself, creates a new PVT or raises the existing
%      singleton*.
%
%      H = PVT returns the handle to a new PVT or the handle to
%      the existing singleton*.
%
%      PVT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PVT.M with the given input arguments.
%
%      PVT('Property','Value',...) creates a new PVT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PVT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PVT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PVT

% Last Modified by GUIDE v2.5 14-Nov-2018 18:13:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PVT_OpeningFcn, ...
                   'gui_OutputFcn',  @PVT_OutputFcn, ...
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


% --- Executes just before PVT is made visible.
function PVT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PVT (see VARARGIN)

% Choose default command line output for PVT
handles.output = hObject;

addpath(genpath('./Library'));
handles.ref_pos = [4627537.2739   119698.4035  4373317.5742];
% Maximum number of epochs to process
handles.nEpoch_max = 10000;
% ENAC reference position (ECEF) - [x,y,z]
handles.ENAC_xyz = 1e6*[4.627536601003540,0.119700014080275,4.373318373560944];
% ENAC reference position (ECEF) - [latitude, longitude, heigth]
handles.ENAC_llh = [43.564758116,1.48173363,203.8171];

set(handles.Plot1,'visible','on');
set(handles.Plot2,'visible','off');
set(handles.Plot3,'visible','off');
set(handles.Plot4,'visible','off');


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PVT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PVT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in HSelection.
function HSelection_Callback(hObject, eventdata, handles)
% hObject    handle to HSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns HSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from HSelection
string = get(handles.HSelection,'String');
value = get(handles.HSelection,'Value');
switch(value)
    case 1
        handles.HMode = 'NLSE';
    case 2
        handles.HMode = 'NWLSE';
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function HSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PseudorangeModelSelection.
function PseudorangeModelSelection_Callback(hObject, eventdata, handles)
% hObject    handle to PseudorangeModelSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PseudorangeModelSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PseudorangeModelSelection


string = get(handles.PseudorangeModelSelection,'String');
value = get(handles.PseudorangeModelSelection,'Value');
switch(value)
    case 1
        handles.PseudorangeModel = 'Code';
    case 2
        handles.PseudorangeModel = 'CodeAndCarrier';
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function PseudorangeModelSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PseudorangeModelSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrowseOBS.
function BrowseOBS_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseOBS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.filename_o,handles.path_filename_o] = uigetfile('.obs','Select OBSERVATION .obs data file');
if handles.filename_o
    set(handles.OBSfile,'String',handles.filename_o);
    handles.path_filename_o = [handles.path_filename_o handles.filename_o];
else
    set(handles.OBSfile,'String','');
    msgbox('.obs file not selected!')
end
guidata(hObject,handles);


% --- Executes on button press in BrowseNAV.
function BrowseNAV_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseNAV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.filename_n,handles.path_filename_n] = uigetfile('.nav','Select NAVIGATION .nav data file');
if handles.filename_o
    set(handles.NAVfile,'String',handles.filename_n);
    handles.path_filename_n = [handles.path_filename_n handles.filename_n];
else
    set(handles.NAVfile,'String','');
    msgbox('.nav file not selected!')
end
guidata(hObject,handles);

% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Run Main
HSelection_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
PseudorangeModelSelection_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
if sum(handles.filename_n ~= 0) && sum(handles.filename_o ~= 0)
    [handles] = Main_Students_PVT(handles);
else
    msgbox('Unable to run because OBS and NAV files were not properly loaded.')
end
guidata(hObject,handles);


% --- Executes on button press in SNRButton.
function SNRButton_Callback(hObject, eventdata, handles)
% hObject    handle to SNRButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Plot1,'visible','on');
set(handles.Plot2,'visible','off');
set(handles.Plot3,'visible','off');
set(handles.Plot4,'visible','off');

axes(handles.Plot1)
hold off
plot(handles.mS1,'linewidth',1)
ylabel('dB')
xlabel('Epoch')
title('Signal to Noise Ratio','fontweight','bold')
legend(strcat('PRN # ', string(find(sum(handles.mS1) ~= 0))))


% --- Executes on button press in OrbitsButton.
function OrbitsButton_Callback(hObject, eventdata, handles)
% hObject    handle to OrbitsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Plot1,'visible','on');
set(handles.Plot2,'visible','off');
set(handles.Plot3,'visible','off');
set(handles.Plot4,'visible','off');

axes(handles.Plot1)
hold off
for index = 1 : handles.INDEX
        EpochToPlot = round(length(handles.SV(index).Result_x)/2);
        scatter3(handles.SV(index).Result_x(EpochToPlot)/1000,handles.SV(index).Result_y(EpochToPlot)/1000,handles.SV(index).Result_z(EpochToPlot)/1000,50,'d','filled','DisplayName',strcat('PRN # ', num2str(handles.SVTracked(index))));
        legend('-DynamicLegend')
        hold all
end
Legend = get(gca,'Legend');
Legend = Legend.String;
earth_sphere("km")
hold on
for index = 1 : handles.INDEX
    plot3(handles.SV(index).Result_x/1000,handles.SV(index).Result_y/1000,handles.SV(index).Result_z/1000,'k.');%,'DisplayName',sprintf("PRN %d", PRN));
    hold on
end
grid on
title('Satellites orbits during data collection','fontweight','bold')
legend(Legend);


% --- Executes on button press in PseduoButton.
function PseduoButton_Callback(hObject, eventdata, handles)
% hObject    handle to PseduoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Plot1,'visible','on');
set(handles.Plot2,'visible','off');
set(handles.Plot3,'visible','off');
set(handles.Plot4,'visible','off');

axes(handles.Plot1)
plot(handles.mC1,'linewidth',1)
ylabel('m')
xlabel('Epoch')
title('Pseudorange','fontweight','bold')
legend(strcat('PRN # ', string(find(sum(handles.mC1) ~= 0))))


% --- Executes on button press in CarrierPhaseButton.
function CarrierPhaseButton_Callback(hObject, eventdata, handles)
% hObject    handle to CarrierPhaseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Plot1,'visible','on');
set(handles.Plot2,'visible','off');
set(handles.Plot3,'visible','off');
set(handles.Plot4,'visible','off');

axes(handles.Plot1)
plot(handles.mL1,'linewidth',1)
ylabel('m')
xlabel('Epoch')
title('Carrier Phase','fontweight','bold')
legend(strcat('PRN # ', string(find(sum(handles.mL1) ~= 0))))


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Plot1,'visible','off');
set(handles.Plot2,'visible','on');
set(handles.Plot3,'visible','on');
set(handles.Plot4,'visible','on');

axes(handles.Plot2)
hold off;
plot(handles.RX_Position_ENU(:,1));
ylabel('m')
xlabel('Epoch')
title('East','fontweight','bold')

axes(handles.Plot3)
hold off;
plot(handles.RX_Position_ENU(:,1));
ylabel('m')
xlabel('Epoch')
title('North','fontweight','bold')

axes(handles.Plot4)
hold off;
plot(handles.RX_Position_ENU(:,1));
ylabel('m')
xlabel('Epoch')
title('Up','fontweight','bold')