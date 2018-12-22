% handles = guidata(hObject);  % Read the struct from the figure
% handles.abc = rand           % Modify the struct
% guidata(hObject, handles);   % Write the modified struct back to the figure

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

% Last Modified by GUIDE v2.5 22-Dec-2018 17:58:51

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
warning('off')

addpath(genpath('./Library'));
handles.ref_pos = [4627537.2739   119698.4035  4373317.5742];
% Maximum number of epochs to process
handles.nEpoch_max = 10000;
% ENAC reference position (ECEF) - [x,y,z]
handles.ENAC_xyz = 1e6*[4.627536601003540,0.119700014080275,4.373318373560944];
% ENAC reference position (ECEF) - [latitude, longitude, heigth]
handles.ENAC_llh = [43.564758116,1.48173363,203.8171];
% Weighting Function Options
handles.StringWType = {'SNR', 'SNR + SV GEO'};

global lambda;
handles.lambda = lambda;

DisplayPlot(hObject,handles,'1','PVT_OpeningFcn')

set(handles.WTypeSelection,'Enable','off')

axes(handles.ENACLogo)
[YourImage, ~, ImageAlpha] = imread('ENAC.png');
image(YourImage, 'AlphaData', ImageAlpha)
% image(imread('ENAC.png'));
axis off

% axes(handles.Liu)
% imshow('Hello.jpg');
handles.handles_BackUp = handles;
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
handles.handles_BackUp = handles;
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
handles.handles_BackUp = handles;
guidata(hObject,handles);


% --- Executes on button press in BrowseNAV.
function BrowseNAV_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseNAV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.filename_n,handles.path_filename_n] = uigetfile('.nav','Select NAVIGATION .nav data file');
if handles.filename_n
    set(handles.NAVfile,'String',handles.filename_n);
    handles.path_filename_n = [handles.path_filename_n handles.filename_n];
else
    set(handles.NAVfile,'String','');
    msgbox('.nav file not selected!')
end
handles.handles_BackUp = handles;
guidata(hObject,handles);

% --- Executes on button press in RunButton.
function RunButton_Callback(hObject, eventdata, handles)
% hObject    handle to RunButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Run Main
DisplayPlot(hObject,handles,'1',[])
handles = handles.handles_BackUp;
PseudorangeModelSelection_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
SVToFilter_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
SmoothNumber_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

if sum(handles.filename_n ~= 0) && sum(handles.filename_o ~= 0)
    [handles] = Main_Students_PVT(handles);
else
    msgbox('Unable to run because OBS and NAV files were not properly loaded.')
end

CarrierSmoothing_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

% Default
handles.RX_Position_ENU = handles.RX_Position_ENU_NLSE_IT;
handles.RX_Position_LLH = handles.RX_Position_LLH_NLSE_IT;
handles.DOP = handles.DOP_NLSE_IT;

% set(handles.WTypeSelection,'Enable','on');
set(handles.SVSelection,'Enable','on')
set(handles.SVSelection,'String','ALL')
set(handles.SmoothNumber,'Enable','on')
set(handles.SmoothNumber,'String','100')
SVSelection_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

guidata(hObject,handles);



% --- Executes on button press in SNRButton.
function SNRButton_Callback(hObject, eventdata, handles)
% hObject    handle to SNRButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1',[])


axes(handles.Plot1)
hold off
plot(handles.Tracked_mS1(:,handles.SVList),'linewidth',1)
ylabel('dB/Hz')
xlabel('Epoch Number')
title('Signal to Noise Ratio','fontweight','bold')
legend(strcat('PRN # ', string(handles.SVTracked(handles.SVList))))
MaxEpoch_Callback(hObject, eventdata, handles);
MinEpoch_Callback(hObject, eventdata, handles);

% --- Executes on button press in OrbitsButton.
function OrbitsButton_Callback(hObject, eventdata, handles)
% hObject    handle to OrbitsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DisplayPlot(hObject,handles,'1',[])


axes(handles.Plot1)
hold off
for index = 1 : handles.INDEX
    EpochToPlotIndex = find(handles.mTracked(:,handles.SVTracked(index)) ~= 0 );
    EpochToPlot = round(( EpochToPlotIndex(1) + EpochToPlotIndex(end) ) / 2);
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
DisplayPlot(hObject,handles,'1',[])


axes(handles.Plot1)
hold off
plot(handles.Tracked_mC1(:,handles.SVList),'linestyle','none','marker','.','markersize',15)
ylabel('meters')
xlabel('Epoch Number')
title('Pseudorange','fontweight','bold')
legend(strcat('PRN # ', string(handles.SVTracked(:,handles.SVList))))
grid on
MaxEpoch_Callback(hObject, eventdata, handles);
MinEpoch_Callback(hObject, eventdata, handles);

% --- Executes on button press in CarrierPhaseButton.
function CarrierPhaseButton_Callback(hObject, eventdata, handles)
% hObject    handle to CarrierPhaseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1',[])


axes(handles.Plot1)
hold off
plot(handles.lambda * handles.Tracked_mL1(:,handles.SVList),'linestyle','none','marker','.','markersize',15)
ylabel('meters')
xlabel('Epoch Number')
title('Carrier Phase','fontweight','bold')
% legend(strcat('PRN # ', string(find(sum(handles.mL1(:,handles.SVList)) ~= 0))),'Location','BestOutside')
legend(strcat('PRN # ', string(handles.SVTracked(:,handles.SVList))))
grid on
MaxEpoch_Callback(hObject, eventdata, handles);
MinEpoch_Callback(hObject, eventdata, handles);

% --- Executes on button press in RXENUButtom.
function RXENUButtom_Callback(hObject, eventdata, handles)
% hObject    handle to RXENUButtom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'234',[])
WNLSE_Callback(hObject, eventdata, handles);
handles = guidata(hObject);

axes(handles.Plot2)
hold off;
plot(handles.RX_Position_ENU(:,1))
hold on
plot(handles.RX_Position_ENU_smoothed(:,1))
ylabel('meters','fontsize',8)
xlabel('Epoch Number','fontsize',8)
title('Error of East','fontweight','bold','fontsize',8)
legend('Unsmoothed','Carrier Smoothing')
grid on
MaxEpoch_Callback(hObject, eventdata, handles);
MinEpoch_Callback(hObject, eventdata, handles);

axes(handles.Plot3)
hold off;
plot(handles.RX_Position_ENU(:,2));
hold on
plot(handles.RX_Position_ENU_smoothed(:,2))
ylabel('meters','fontsize',8)
xlabel('Epoch Number','fontsize',8)
title('Error of North','fontweight','bold','fontsize',8)
legend('Unsmoothed','Carrier Smoothing')
grid on
MaxEpoch_Callback(hObject, eventdata, handles);
MinEpoch_Callback(hObject, eventdata, handles);

axes(handles.Plot4)
hold off;
plot(handles.RX_Position_ENU(:,3));
hold on
plot(handles.RX_Position_ENU_smoothed(:,3))
ylabel('meters','fontsize',8)
xlabel('Epoch Number','fontsize',8)
title('Error of Up','fontweight','bold','fontsize',8)
legend('Unsmoothed','Carrier Smoothing')
grid on
MaxEpoch_Callback(hObject, eventdata, handles);
MinEpoch_Callback(hObject, eventdata, handles);

% --- Executes on button press in ElevationButton.
function ElevationButton_Callback(hObject, eventdata, handles)
% hObject    handle to ElevationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1',[])


axes(handles.Plot1)
hold off
plot(handles.elevation_SV(:,handles.SVList),'linestyle','none','marker','.','markersize',15)
ylabel('Degrees')
xlabel('Epoch Number')
title('Satellite ELEVATION with respect to RX','fontweight','bold')
legend(strcat('PRN # ', string(handles.SVTracked(handles.SVList))))
grid on
MaxEpoch_Callback(hObject, eventdata, handles);
MinEpoch_Callback(hObject, eventdata, handles);

% --- Executes on button press in AzimuthButton.
function AzimuthButton_Callback(hObject, eventdata, handles)
% hObject    handle to AzimuthButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1',[])


axes(handles.Plot1)
hold off
plot(handles.azimuth_SV(:,handles.SVList),'linestyle','none','marker','.','markersize',15)
ylabel('Degrees')
xlabel('Epoch Number')
title('Satellite AZIMUTH with respect to RX','fontweight','bold')
legend(strcat('PRN # ', string(handles.SVTracked(handles.SVList))))
grid on
MaxEpoch_Callback(hObject, eventdata, handles);
MinEpoch_Callback(hObject, eventdata, handles);

% --- Executes on button press in LATLONButton.
function LATLONButton_Callback(hObject, eventdata, handles)
% hObject    handle to LATLONButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


DisplayPlot(hObject,handles,'1',[])


axes(handles.Plot1)
hold off

for index = 1 : handles.INDEX
    EpochToPlotIndex = find(handles.mTracked(:,handles.SVTracked(index)) ~= 0 );
    EpochToPlot = round(( EpochToPlotIndex(1) + EpochToPlotIndex(end) ) / 2);
    scatter((handles.SV(index).llh(EpochToPlot,2)),(handles.SV(index).llh(EpochToPlot,1)),50,'d','filled','DisplayName',strcat('PRN # ', num2str(handles.SVTracked(index))));
    legend('-DynamicLegend','Location','BestOutside')
    hold all
end
Legend = get(gca,'Legend');
Legend = Legend.String;
for index = 1 : handles.INDEX
    plot((handles.SV(index).llh(:,2)),(handles.SV(index).llh(:,1)),'k.');%,'DisplayName',sprintf("PRN %d", PRN));
    hold on
end
grid on
legend(Legend);
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
title('Footprint Latitude - Longitude','fontweight','bold')


% --- Executes on button press in DoPButton.
function DoPButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoPButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1',[]);
WNLSE_Callback(hObject, eventdata, handles);
handles = guidata(hObject);

RMS_Errors = [sqrt(handles.DOP.EDOP.^2 + handles.DOP.NDOP.^2); ...
    sqrt(handles.DOP.VDOP.^2);...
    sqrt(handles.DOP.EDOP.^2 + handles.DOP.NDOP.^2 + handles.DOP.VDOP.^2);...
    sqrt(handles.DOP.EDOP.^2 + handles.DOP.NDOP.^2 + handles.DOP.VDOP.^2 + handles.DOP.TDOP.^2); ...
    handles.DOP.TDOP];

% RMS_Errors_CS = [sqrt(handles.DOP_smoothed.EDOP.^2 + handles.DOP_smoothed.NDOP.^2); ...
%     sqrt(handles.DOP_smoothed.VDOP.^2);...
%     sqrt(handles.DOP_smoothed.EDOP.^2 + handles.DOP_smoothed.NDOP.^2 + handles.DOP_smoothed.VDOP.^2);...
%     sqrt(handles.DOP_smoothed.EDOP.^2 + handles.DOP_smoothed.NDOP.^2 + handles.DOP_smoothed.VDOP.^2 + handles.DOP_smoothed.TDOP.^2); ...
%     handles.DOP_smoothed.TDOP];

axes(handles.Plot1)
hold off
plot(RMS_Errors.' ,'linewidth',1);
legend('HDOP','VDOP','PDOP','GDOP','TDOP')
title('Dilution of Precision)')
xlabel('Epoch Number')
grid on;
MaxEpoch_Callback(hObject, eventdata, handles);
MinEpoch_Callback(hObject, eventdata, handles);

% --- Executes on button press in RMSButton.
function RMSButton_Callback(hObject, eventdata, handles)
% hObject    handle to RMSButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1',[])
WNLSE_Callback(hObject, eventdata, handles);
handles = guidata(hObject);
% set(handles.WTypeSelection,'Enable','on');


RMS_Errors = [sqrt(handles.DOP.EDOP.^2 + handles.DOP.NDOP.^2); ...
    sqrt(handles.DOP.EDOP.^2 + handles.DOP.NDOP.^2 + handles.DOP.VDOP.^2);...
    sqrt(handles.DOP.EDOP.^2 + handles.DOP.NDOP.^2 + handles.DOP.VDOP.^2 + handles.DOP.TDOP.^2); ...
    handles.DOP.TDOP];

axes(handles.Plot1)
hold off
plot(RMS_Errors.' ,'linewidth',1);
hold on
legend('HDOP','PDOP','GDOP','TDOP')
title('Dilusion of Precision (without considering RMS errors (CS, PD, MP and noise))')
xlabel('Epoch Number')
ylabel('Meters')
grid on;

% --- Executes on button press in CodeMinusCarrierButton.
function CodeMinusCarrierButton_Callback(hObject, eventdata, handles)
% hObject    handle to CodeMinusCarrierButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DisplayPlot(hObject,handles,'1',[])

CMC = handles.Tracked_mC1 - handles.Tracked_mL1*handles.lambda;
CMC_smoothed = handles.Tracked_smoothed - handles.Tracked_mL1*handles.lambda;

axes(handles.Plot1)
hold off
plot(CMC(:,handles.SVList),'linewidth',1);
ylabel('meters')
xlabel('Epoch Number')
title('Code-minus-Carrier','fontweight','bold')
Legend = strcat('PRN # ', string(handles.SVTracked(handles.SVList)));
hold on
plot(CMC_smoothed(:,handles.SVList),'linewidth',1);
legend(Legend)
grid on
MaxEpoch_Callback(hObject, eventdata, handles);
MinEpoch_Callback(hObject, eventdata, handles);


% --- Executes on button press in PositionScatterDisplayButton.
function PositionScatterDisplayButton_Callback(hObject, eventdata, handles)
% hObject    handle to PositionScatterDisplayButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1',[])


WNLSE = get(handles.WNLSE,'value');
WTypeSelection_Callback(hObject, eventdata, handles);
handles = guidata(hObject);
if WNLSE
    String = get(handles.WTypeSelection,'String');
    ValueToString = get(handles.WTypeSelection,'Value');
else
    String = {'1'};
    ValueToString = 1;
end
LEGEND = [{'Raw','with IT','Carrier Smoothing'}];
Type = get(handles.WTypeSelection,'Value');

East = handles.RX_Position_ENU(:,1);
North = handles.RX_Position_ENU(:,2);
std_East = std(East);
std_North = std(North);

East_CS = handles.RX_Position_ENU_smoothed(:,1);
North_CS = handles.RX_Position_ENU_smoothed(:,2);
std_East_CS = std(East_CS);
std_North_CS = std(North_CS);

if ~WNLSE
    East_nonIT = handles.RX_Position_ENU_NLSE(:,1);
    North_nonIT = handles.RX_Position_ENU_NLSE(:,2);
    string = 'LSE';
elseif WNLSE
    East_nonIT = handles.RX_Position_ENU_W(Type).NLSE(:,1);
    North_nonIT = handles.RX_Position_ENU_W(Type).NLSE(:,2);
    string = 'Weighted LSE';
end
std_East_nonIT = std(East_nonIT);
std_North_nonIT = std(North_nonIT);

sigmas = [std_East, std_East_nonIT, std_East_CS, std_North, std_North_nonIT, std_North_CS];
Titles = [{'East'}, {'East no IT'}, {'East Carrier Smoothing'}, {'North'}, {'North non IT'}, {'North Carrier Smoothing'}];
if WNLSE
    Names = String(Type);
else
    Names = {'NLSE'};
end
for index = 1 : size(sigmas,1)
    Rows(index,:) = num2cell(sigmas(index,:));
end

clc
fprintf('Error Standard Deviation with %s algorithm:', String{ValueToString});
disp(' ')
fprintf('------------------------------------------------------');
STD = [{'Weight \ Coordinate'},Titles;Names,Rows]

axes(handles.Plot1)
hold off
plot(East_nonIT,North_nonIT,'b.','linewidth',1);
hold on
plot(East,North,'g.','linewidth',1);
hold on
plot(East_CS,North_CS,'r.','linewidth',1);
title(sprintf('Recevier Position for %s and weight %s (Raw and I/T Corrected).',string, String{ValueToString}))
xlabel('East (m)')
ylabel('North (m)')
legend(LEGEND)
% xmin = max(min(East_nonIT),min(East));
% xmax = min(max(East_nonIT),max(East));
% ymin = max(min(North_nonIT),min(North));
% ymax = min(max(North_nonIT),max(North));
% xmin = max(xmin,min(East_CS));
% xmax = min(xmax,max(East_CS));
% ymin = max(ymin,min(East_CS));
% ymax = min(ymax,max(East_CS));
% x = max(abs(xmin), abs(xmax));
% y = max(abs(ymin), abs(ymax));
% xlim([-x x])
% ylim([-y y])
grid on



% --- Executes on button press in WNLSE.
function WNLSE_Callback(hObject, eventdata, handles)
% hObject    handle to WNLSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of WNLSE

WNLSE = get(handles.WNLSE,'value');

if ~WNLSE
    set(handles.WTypeSelection,'Enable','off');
    handles.RX_Position_ENU = handles.RX_Position_ENU_NLSE_IT;
    handles.RX_Position_LLH = handles.RX_Position_LLH_NLSE_IT;
    handles.DOP = handles.DOP_NLSE_IT;
    handles.RX_Position_ENU_smoothed = handles.RX_Position_ENU_UWsmoothed;
    handles.RX_Position_LLH_smoothed = handles.RX_Position_LLH_UWsmoothed;
    handles.DOP_smoothed = handles.DOP_UWsmoothed;
else
    Type = get(handles.WTypeSelection,'value');
    set(handles.WTypeSelection,'Enable','on');
    handles.RX_Position_ENU = handles.RX_Position_ENU_W(Type).NLSE_IT;
    handles.RX_Position_LLH = handles.RX_Position_LLH_W(Type).NLSE_IT;
    handles.DOP = handles.DOP_W(Type).NLSE_IT;
    handles.RX_Position_ENU_smoothed = handles.RX_Position_ENU_Wsmoothed(Type).sm;
    handles.RX_Position_LLH_smoothed = handles.RX_Position_LLH_Wsmoothed(Type).sm;
    handles.DOP_smoothed = handles.DOP_Wsmoothed(Type).sm;
end

guidata(hObject,handles);



function DisplayPlot(hObject,handles,PlotsToDisplay,Call)

if strcmp(PlotsToDisplay,'1')
    set(handles.Plot1,'visible','on');
    set(handles.Plot2,'visible','off');
    cla(handles.Plot2);
    set(handles.Plot3,'visible','off');
    cla(handles.Plot3);
    set(handles.Plot4,'visible','off');
    cla(handles.Plot4);
elseif strcmp(PlotsToDisplay,'234')
    set(handles.Plot1,'visible','off');
    cla(handles.Plot1);
    set(handles.Plot2,'visible','on');
    set(handles.Plot3,'visible','on');
    set(handles.Plot4,'visible','on');
end

% set(handles.WTypeSelection,'Enable','on');

if ~strcmp(Call,'PVT_OpeningFcn')
    if strcmp(get(hObject,'string'),'Elevation')...
            || strcmp(get(hObject,'string'),'Azimuth')...
            || strcmp(get(hObject,'string'),'Pseudorange')...
            || strcmp(get(hObject,'string'),'Carrier Phase')...
            || strcmp(get(hObject,'string'),'Code - Carrier')...
            || strcmp(get(hObject,'string'),'SNR')...
        set(handles.MinEpoch,'Enable','on')
        set(handles.MaxEpoch,'Enable','on')
        set(handles.SVSelection,'Enable','on')
    else
        set(handles.MinEpoch,'Enable','off')
        set(handles.MaxEpoch,'Enable','off')
        set(handles.SVSelection,'Enable','off')
    end
    if strcmp(get(hObject,'string'),'Code - Carrier')...
            || strcmp(get(hObject,'string'),'Position Estimates')...
            || strcmp(get(hObject,'string'),'ENU')
        set(handles.SmoothNumber,'Enable','on')
    else
        set(handles.SmoothNumber,'Enable','off')
    end
end

%     if strcmp(get(hObject,'string'),'DoP') || strcmp(get(hObject,'string'),'RMS')
%     else
%         set(handles.WTypeSelection,'Enable','off');
%         set(handles.WTypeSelection,'Value',1);
%     end



guidata(hObject,handles);


% --- Executes on selection change in WTypeSelection.
function WTypeSelection_Callback(hObject, eventdata, handles)
% hObject    handle to WTypeSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns WTypeSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from WTypeSelection
string = get(handles.WTypeSelection,'String');
value = get(handles.WTypeSelection,'Value');
handles.WType = string{value};
WNLSE_Callback(hObject, eventdata, handles);
handles = guidata(hObject);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function WTypeSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WTypeSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MinEpoch_Callback(hObject, eventdata, handles)
% hObject    handle to MinEpoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinEpoch as text
%        str2double(get(hObject,'String')) returns contents of MinEpoch as a double
CurrentLimits = xlim;
MinEpoch = str2num(get(handles.MinEpoch,'String'));
if ~numel(MinEpoch)
    MinEpoch = 1;
end
xlim([MinEpoch, CurrentLimits(2)])

% --- Executes during object creation, after setting all properties.
function MinEpoch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinEpoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxEpoch_Callback(hObject, eventdata, handles)
% hObject    handle to MaxEpoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxEpoch as text
%        str2double(get(hObject,'String')) returns contents of MaxEpoch as a double
CurrentLimits = xlim;
MaxEpoch = str2num(get(handles.MaxEpoch,'String'));
if ~numel(MaxEpoch)
    MaxEpoch = handles.Nb_Epoch;
end
xlim([CurrentLimits(1), MaxEpoch])

% --- Executes during object creation, after setting all properties.
function MaxEpoch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxEpoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SVSelection_Callback(hObject, eventdata, handles)
% hObject    handle to SVSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SVSelection as text
%        str2double(get(hObject,'String')) returns contents of SVSelection as a double

% handles = guidata(hObject);

STR = get(handles.SVSelection,'string');
Index = find(STR == ',');
if numel(Index)
    SVList(1) = str2num(STR(1:Index(1)-1));
    for count = 1 : length(Index) - 1
        SVList(count + 1) = str2num(STR(Index(count)+1:Index(count+1)-1));
    end
    SVList(end + 1) = str2num(STR(Index(end)+1:end));
else
    SVList = str2num(STR);
    if ~numel(SVList)
        SVList = handles.SVTracked;
        set(handles.SVSelection,'String','ALL')
    end
end

index = 0;
for counter = 1 : length(SVList)
    Lag = find(handles.SVTracked == SVList(counter));
    if Lag
        index = index + 1;
        List(index) = Lag;
    end
end

if ~exist('List','var')
    SVList = handles.SVTracked;
    set(handles.SVSelection,'String','ALL')
    index = 0;
    for counter = 1 : length(SVList)
        Lag = find(handles.SVTracked == SVList(counter));
        if Lag
            index = index + 1;
            List(index) = Lag;
        end
    end
end
handles.SVList = List;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function SVSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SVSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadButton.
function LoadButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.LoadFile,handles.path_LoadFile] = uigetfile('.mat','Select file to load');
set(handles.MessageBox,'String','Loading data...')
pause(0.001)
if handles.LoadFile
    handles.path_LoadFile = [handles.path_LoadFile handles.LoadFile];
    LoadedFile = handles.LoadFile;
    load(handles.path_LoadFile)
    set(handles.MessageBox,'String',['LOADED ' LoadedFile]);
else
    set(handles.MessageBox,'String','');
    msgbox('File to load not selected!')
end

if handles.LoadFile
    FieldNames = fieldnames(HandlesSaved);
    for index = 1 : numel(FieldNames)
        handles = setfield(handles,cell2mat(FieldNames(index)),getfield(HandlesSaved,cell2mat(FieldNames(index))));
    end
    set(handles.MessageBox,'String',[LoadedFile ' Loaded!'])
    pause(0.001)
    % Default
    handles.RX_Position_ENU = handles.RX_Position_ENU_NLSE_IT;
    handles.RX_Position_LLH = handles.RX_Position_LLH_NLSE_IT;
    handles.DOP = handles.DOP_NLSE_IT;
    
    set(handles.WTypeSelection,'Enable','on');
    set(handles.SVSelection,'Enable','on')
    set(handles.SVSelection,'String','ALL')
    SVSelection_Callback(hObject, eventdata, handles)
    handles = guidata(hObject);
    
    guidata(hObject,handles);
    OrbitsButton_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in SaveButton.
% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FieldNames = fieldnames(handles);
HandlesSaved = struct();

StartToSave = 0;
for index = 1 : numel(FieldNames)
    if strcmp(cell2mat(FieldNames(index)),'ref_pos')
        StartToSave = 1;
    end
    if StartToSave && ~strcmp(cell2mat(FieldNames(index)),'handles_BackUp')
        HandlesSaved = setfield(HandlesSaved,cell2mat(FieldNames(index)),getfield(handles,cell2mat(FieldNames(index))));
    end
end

[handles.SaveDirectory] = uigetdir('Select directory');
DataSetName = {''};

if handles.SaveDirectory
    DataSetName = inputdlg('Write a name for the data to save');
    if numel(cell2mat(DataSetName))
        set(handles.MessageBox,'String','Saving data...')
        pause(0.001)
        save([handles.SaveDirectory '\' cell2mat(DataSetName) '.mat'],'HandlesSaved');
        set(handles.MessageBox,'String','Data Saved!')
        filename = [handles.SaveDirectory '\' cell2mat(DataSetName) '.kml'];
        iconFilename = fullfile(pwd, 'icon18.png');
        kmlwritepoint(filename,handles.RX_Position_LLH_W(1).NLSE_IT(:,1),...
            handles.RX_Position_LLH_W(2).NLSE_IT(:,2),handles.RX_Position_LLH_W(2).NLSE_IT(:,3),...
            'Icon', iconFilename, 'IconScale',0.2, 'Name', blanks(handles.Nb_Epoch), 'Color', 'blue')
        pause(0.001)
    end
end
if ~numel(handles.SaveDirectory) || ~numel(cell2mat(DataSetName))
    msgbox('Data not saved')
end
guidata(hObject,handles);

% --- Executes on button press in GoogleMapsButton.
function GoogleMapsButton_Callback(hObject, eventdata, handles)
% hObject    handle to GoogleMapsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
String = ['https://www.google.fr/maps/@' num2str(handles.RX_Position_LLH(end,1),'%.17f') ',' num2str(handles.RX_Position_LLH(end,2),'%.17f')];
web(String)


function SVToFilter_Callback(hObject, eventdata, handles)
% hObject    handle to SVToFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SVToFilter as text
%        str2double(get(hObject,'String')) returns contents of SVToFilter as a double

STR = get(handles.SVToFilter,'string');
Index = find(STR == ',');
if numel(Index)
    SVListFilter(1) = str2num(STR(1:Index(1)-1));
    for count = 1 : length(Index) - 1
        SVListFilter(count + 1) = str2num(STR(Index(count)+1:Index(count+1)-1));
    end
    SVListFilter(end + 1) = str2num(STR(Index(end)+1:end));
else
    SVListFilter = str2num(STR);
    if ~numel(SVListFilter)
        SVListFilter = 0;
        set(handles.SVToFilter,'String','NONE')
    end
end

handles.SVListFilter = SVListFilter;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function SVToFilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SVToFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SkyPlot.
function SkyPlot_Callback(hObject, eventdata, handles)
% hObject    handle to SkyPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1',[]);

axes(handles.Plot1)
hold off
skyPlot(handles)



function SmoothNumber_Callback(hObject, eventdata, handles)
% hObject    handle to SmoothNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SmoothNumber as text
%        str2double(get(hObject,'String')) returns contents of SmoothNumber as a double
handles.smoothnumber = str2num(get(handles.SmoothNumber,'String'));
if ~numel(handles.smoothnumber)
    handles.smoothnumber = 100;
    set(handles.SmoothNumber,'String','100')
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function SmoothNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SmoothNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CarrierSmoothing.
function CarrierSmoothing_Callback(hObject, eventdata, handles)
% hObject    handle to CarrierSmoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Caculate Carrier Smoothing
handles.CMC = handles.mC1 - handles.mL1*handles.lambda;
[handles.smoothed,handles.cycle_slip] = Carrier_Smoothing(handles.mC1,handles.mL1,handles.smoothnumber,handles.CMC);

%Recaculate Tracked_smoothed
handles.Tracked_smoothed = zeros(handles.Nb_Epoch,length(handles.SVTracked));
for num_SV=1:length(handles.SVTracked)
    handles.Tracked_smoothed(:,num_SV) = handles.smoothed(:,handles.SVTracked(num_SV));
end
handles.Tracked_smoothed(handles.Tracked_smoothed==0) = nan;


%Smoothed Non Linear LSE
[handles.RX_Position_XYZ_UWsmoothed, handles.RX_ClockError_UWsmoothed, handles.Matrix_UWsmoothed] = RX_Position_and_Clock(handles.Result,handles.smoothed,handles.mS1,handles.Nb_Epoch,handles.vNb_Sat,'NLSE',0,handles.Tiono,handles.Ttropo,7,[],handles.MessageBox);
[handles.RX_Position_LLH_UWsmoothed, handles.RX_Position_ENU_UWsmoothed, handles.Matrix_UWsmoothed, handles.DOP_UWsmoothed] = RX_Position_LLH_ENU(handles.RX_Position_XYZ_UWsmoothed,handles.Nb_Epoch,handles.Matrix_UWsmoothed);

% Smoothed Weighted (SNR) Non Linear LSE - ATMOSPHERIC CORRECTION
[handles.RX_Position_XYZ_Wsmoothed(1).sm, handles.RX_ClockError_Wsmoothed(1).sm, handles.Matrix_Wsmoothed(1).sm] = RX_Position_and_Clock(handles.Result,handles.smoothed,handles.mS1,handles.Nb_Epoch,handles.vNb_Sat,'NWLSE',1,handles.Tiono,handles.Ttropo,8,[],handles.MessageBox);
[handles.RX_Position_LLH_Wsmoothed(1).sm, handles.RX_Position_ENU_Wsmoothed(1).sm, handles.Matrix_Wsmoothed(1).sm, handles.DOP_Wsmoothed(1).sm] = RX_Position_LLH_ENU(handles.RX_Position_XYZ_Wsmoothed(1).sm,handles.Nb_Epoch,handles.Matrix_Wsmoothed(1).sm);

% Smoothed Weighted (SNR + ELEVATION) Non Linear LSE - ATMOSPHERIC CORRECTION
[handles.RX_Position_XYZ_Wsmoothed(2).sm, handles.RX_ClockError_Wsmoothed(2).sm, handles.Matrix_Wsmoothed(2).sm] = RX_Position_and_Clock(handles.Result,handles.smoothed,handles.mS1,handles.Nb_Epoch,handles.vNb_Sat,'NWLSE',2,handles.Tiono,handles.Ttropo,9,handles.Elevation_Azimuth,handles.MessageBox);
[handles.RX_Position_LLH_Wsmoothed(2).sm, handles.RX_Position_ENU_Wsmoothed(2).sm, handles.Matrix_Wsmoothed(2).sm, handles.DOP_Wsmoothed(2).sm] = RX_Position_LLH_ENU(handles.RX_Position_XYZ_Wsmoothed(2).sm,handles.Nb_Epoch,handles.Matrix_Wsmoothed(2).sm);

guidata(hObject,handles);


