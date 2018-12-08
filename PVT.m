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

% Last Modified by GUIDE v2.5 05-Dec-2018 19:56:30

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
% Weighting Function Options
handles.StringWType = {'SNR', 'SNR + SV GEO'};

c = 299792458;
f = 1575.42e6;
handles.lambda = c/f;

DisplayPlot(hObject,handles,'1')

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

handles = handles.handles_BackUp;
PseudorangeModelSelection_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

if sum(handles.filename_n ~= 0) && sum(handles.filename_o ~= 0)
    [handles] = Main_Students_PVT(handles);
else
    msgbox('Unable to run because OBS and NAV files were not properly loaded.')
end

% Default
handles.RX_Position_ENU = handles.RX_Position_ENU_NLSE_IT;
handles.DOP = handles.DOP_NLSE_IT;
% MeasurementJump = 0;
% 
% for index = 1 : size(handles.CMC,2)
%     tmp = find(handles.CMC(:,Index)>1e3);
%     MeasurementJump(Index) = tmp(1);
% end
% handles.MeasurementJump = min(MeasurementJump);

% WTypeSelection_Callback(hObject, eventdata, handles)
% handles = guidata(hObject);

% guidata(hObject,handles);

set(handles.WTypeSelection,'Enable','on');
set(handles.SVSelection,'Enable','on')
set(handles.SVSelection,'String','ALL')
SVSelection_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

guidata(hObject,handles);



% --- Executes on button press in SNRButton.
function SNRButton_Callback(hObject, eventdata, handles)
% hObject    handle to SNRButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1')


axes(handles.Plot1)
hold off
plot(handles.Tracked_mS1(:,handles.SVList),'linewidth',1)
ylabel('dB/Hz')
xlabel('Epoch')
title('Signal to Noise Ratio','fontweight','bold')
legend(strcat('PRN # ', string(handles.SVTracked(handles.SVList))),'Location','BestOutside')

% --- Executes on button press in OrbitsButton.
function OrbitsButton_Callback(hObject, eventdata, handles)
% hObject    handle to OrbitsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DisplayPlot(hObject,handles,'1')


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
DisplayPlot(hObject,handles,'1')


axes(handles.Plot1)
hold off
plot(handles.Tracked_mC1(:,handles.SVList),'linewidth',1,'linestyle','none','marker','.')
ylabel('m')
xlabel('Epoch')
title('Pseudorange','fontweight','bold')
legend(strcat('PRN # ', string(handles.SVTracked(:,handles.SVList))))


% --- Executes on button press in CarrierPhaseButton.
function CarrierPhaseButton_Callback(hObject, eventdata, handles)
% hObject    handle to CarrierPhaseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1')


axes(handles.Plot1)
hold off
plot(handles.lambda * handles.Tracked_mL1(:,handles.SVList),'linewidth',1,'linestyle','none','marker','.')
ylabel('meters')
xlabel('Epoch')
title('Carrier Phase','fontweight','bold')
% legend(strcat('PRN # ', string(find(sum(handles.mL1(:,handles.SVList)) ~= 0))),'Location','BestOutside')
legend(strcat('PRN # ', string(handles.SVTracked(:,handles.SVList))))


% --- Executes on button press in RXENUButtom.
function RXENUButtom_Callback(hObject, eventdata, handles)
% hObject    handle to RXENUButtom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'234')
WNLSE_Callback(hObject, eventdata, handles);
handles = guidata(hObject);

axes(handles.Plot2)
hold off;
plot(handles.RX_Position_ENU(:,1));
ylabel('m')
xlabel('Epoch')
title('East','fontweight','bold')

axes(handles.Plot3)
hold off;
plot(handles.RX_Position_ENU(:,2));
ylabel('m')
xlabel('Epoch')
title('North','fontweight','bold')

axes(handles.Plot4)
hold off;
plot(handles.RX_Position_ENU(:,3));
ylabel('m')
xlabel('Epoch')
title('Up','fontweight','bold')


% --- Executes on button press in ElevationButton.
function ElevationButton_Callback(hObject, eventdata, handles)
% hObject    handle to ElevationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1')


axes(handles.Plot1)
hold off
plot(handles.elevation_SV(:,handles.SVList),'linewidth',2)
ylabel('Degrees')
xlabel('Epoch')
title('Elevation between RX and SV','fontweight','bold')
legend(strcat('PRN # ', string(handles.SVTracked(handles.SVList))),'Location','BestOutside')
grid on

% --- Executes on button press in AzimuthButton.
function AzimuthButton_Callback(hObject, eventdata, handles)
% hObject    handle to AzimuthButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1')


axes(handles.Plot1)
hold off
plot(handles.azimuth_SV(:,handles.SVList),'linewidth',2)
ylabel('Degrees')
xlabel('Epoch')
title('Azimuth between RX and SV','fontweight','bold')
legend(strcat('PRN # ', string(handles.SVTracked(handles.SVList))))
grid on

% --- Executes on button press in LATLONButton.
function LATLONButton_Callback(hObject, eventdata, handles)
% hObject    handle to LATLONButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if strcmp(get(hObject,'string'),'LAT - LON')

DisplayPlot(hObject,handles,'1')


axes(handles.Plot1)
hold off

% for index = 1 : handles.INDEX
%     for epoch = 1 : length(handles.SV(index).Result_x)
%         [handles.SV(index).llh(epoch,:)] = xyz_2_lla_PVT( [handles.SV(index).Result_x(epoch), handles.SV(index).Result_y(epoch), handles.SV(index).Result_z(epoch)] );
%     end
% end
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
title('Sky Plot Latitude - Longitude','fontweight','bold')


% --- Executes on button press in DoPButton.
function DoPButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoPButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1');
WNLSE_Callback(hObject, eventdata, handles);
handles = guidata(hObject);
set(handles.WTypeSelection,'Enable','on');

% axes(handles.Plot1)
% hold off
% plot(...
%     [handles.DOP.EDOP; handles.DOP.NDOP; handles.DOP.VDOP; handles.DOP.TDOP].'...
%     ,'linewidth',1);
% hold on
% legend('EDOP','NDOP','VDOP','TDOP')
% title('Dilusion of Precision (considering RMS errors (CS, PD, MP and noise))')
% xlabel('Epoch')
% ylabel('Meters')
% grid on

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
xlabel('Epoch')
grid on;

% --- Executes on button press in RMSButton.
function RMSButton_Callback(hObject, eventdata, handles)
% hObject    handle to RMSButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1')
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
xlabel('Epoch')
ylabel('Meters')
grid on;

% --- Executes on button press in CodeMinusCarrierButton.
function CodeMinusCarrierButton_Callback(hObject, eventdata, handles)
% hObject    handle to CodeMinusCarrierButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DisplayPlot(hObject,handles,'1')


CMC = handles.Tracked_mC1 - handles.Tracked_mL1*handles.lambda;

axes(handles.Plot1)
hold off
plot(CMC(:,handles.SVList),'linewidth',1);
ylabel('meters')
xlabel('Epoch')
title('Code-minus-Carrier','fontweight','bold')
legend(strcat('PRN # ', string(handles.SVTracked(handles.SVList))))
grid on


% --- Executes on button press in PositionScatterDisplayButton.
function PositionScatterDisplayButton_Callback(hObject, eventdata, handles)
% hObject    handle to PositionScatterDisplayButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DisplayPlot(hObject,handles,'1')


WNLSE = get(handles.WNLSE,'value');
% String = get(handles.WTypeSelection,'String');
String = handles.StringWType;
Value = 1;
if ~WNLSE
%     Value = 1;
    LEGEND = [{'Raw','with IT'}];
else
%     Value = get(handles.WTypeSelection,'Value');
    LEGEND = [{[String{1} ' Raw']}, {[String{2} ' Raw']}, {[String{1} ' with IT']}, {[String{2} ' with IT']}];
end
VectorElements = [1:numel(String)]; % 1 2 3, Value  = 2
tmp = find(VectorElements ~= Value); % 1 3
IndexVector = VectorElements(tmp); % 1 2 3 (in 1 3) = 1 3

    East(:,Value) = handles.RX_Position_ENU(:,1); % column 1
    North(:,Value) = handles.RX_Position_ENU(:,2); % column 1
    std_East(Value) = std(East);
    std_North(Value) = std(North);

if ~WNLSE
    East_nonIT(:,Value) = handles.RX_Position_ENU_NLSE(:,1);
    North_nonIT(:,Value) = handles.RX_Position_ENU_NLSE(:,2);
    std_East_nonIT(Value) = std(East_nonIT);
    std_North_nonIT(Value) = std(North_nonIT);
    string = 'NLSE';
elseif WNLSE
    Index = 0;
    East_nonIT(:,1) = handles.RX_Position_ENU_W(Value).NLSE(:,1);
    North_nonIT(:,1) = handles.RX_Position_ENU_W(Value).NLSE(:,2);
    std_East_nonIT(1) = std(East_nonIT);
    std_North_nonIT(1) = std(North_nonIT);
    for Type = IndexVector % [1 3] + 1 = [2 4]
        East(:,Type) = handles.RX_Position_ENU_W(Type).NLSE_IT(:,1);
        North(:,Type) = handles.RX_Position_ENU_W(Type).NLSE_IT(:,2);
        East_nonIT(:,Type) = handles.RX_Position_ENU_W(Type).NLSE(:,1);
        North_nonIT(:,Type) = handles.RX_Position_ENU_W(Type).NLSE(:,2);
        
        std_East(Type) = std(East(:,Type));
        std_North(Type) = std(North(:,Type));
        std_East_nonIT(Type) = std(East_nonIT(:,Type));
        std_North_nonIT(Type) = std(North_nonIT(:,Type));
    end
    string = 'Weighted NLSE';
end

sigmas = [std_East.', std_East_nonIT.', std_North.', std_North_nonIT.'];
Titles = [{'East'}, {'East no IT'}, {'North'}, {'North non IT'}];
if WNLSE
    Names = [String(Value), String(IndexVector)];
else
    Names = [{'NLSE'}];
end
for index = 1 : size(sigmas,1)
    Rows(index,:) = num2cell(sigmas(index,:));
end

clc
fprintf('Error Standard Deviation with %s algorithm:', string);
disp(' ')
fprintf('------------------------------------------------------');
STD = [{'Weight \ Coordinate'},Titles;Names.',Rows]

axes(handles.Plot1)
% for k = 1 : 3
hold off
plot(East_nonIT,North_nonIT,'.','linewidth',1);
hold on
plot(East,North,'.','linewidth',1);
legend(LEGEND)
title(sprintf('Recevier Position for %s (Raw and I/T Corrected).', string))
xlabel('East (m)')
ylabel('North (m)')
grid on
% end



% --- Executes on button press in WNLSE.
function WNLSE_Callback(hObject, eventdata, handles)
% hObject    handle to WNLSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of WNLSE

WNLSE = get(handles.WNLSE,'value');
% handles.WType = 'SNR';
% switch handles.WType
%     case 'SNR'
%         Type = 1;
%     case 'SNR + SV GEO'
%         Type = 2;
% end
if ~WNLSE
    set(handles.WTypeSelection,'Enable','off');
    handles.RX_Position_ENU = handles.RX_Position_ENU_NLSE_IT;
    handles.DOP = handles.DOP_NLSE_IT;
else
    Type = get(handles.WTypeSelection,'value');
    set(handles.WTypeSelection,'Enable','on');
    handles.RX_Position_ENU = handles.RX_Position_ENU_W(Type).NLSE_IT;
    handles.DOP = handles.DOP_W(Type).NLSE_IT;
end

guidata(hObject,handles);

% --- Executes on button press in AtmosphericCorrections.
function AtmosphericCorrections_Callback(hObject, eventdata, handles)
% hObject    handle to AtmosphericCorrections (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AtmosphericCorrections


% % --- Executes on button press in PositionEstimatesWarmMapButton.
% function PositionEstimatesWarmMapButton_Callback(hObject, eventdata, handles)
% % hObject    handle to PositionEstimatesWarmMapButton (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% DisplayPlot(hObject,handles,'5678')
% 
% WNLSE = get(handles.WNLSE,'value');
% 
%     
% %     Up = handles.RX_Position_ENU(:,3);
% 
% if ~WNLSE
%     East_nonIT = handles.RX_Position_ENU_NLSE(:,1);
%     North_nonIT = handles.RX_Position_ENU_NLSE(:,2);
% %     Up_nonIT = handles.RX_Position_ENU_NLSE(:,3);
%     string = 'NLSE';
% elseif WNLSE
%     
%     East_nonIT = handles.RX_Position_ENU_W(Type).NLSE(:,1);
%     North_nonIT = handles.RX_Position_ENU_W(Type).NLSE(:,2);
% %     Up_nonIT = handles.RX_Position_ENU_NWLSE(:,3);
%     string = 'Weighted NLSE';
% end
% 
% title(sprintf('Recevier Position for %s (Raw [UP] and I/T Corrected [BOTTOM])', string))
% axes(handles.Plot5)
% hold off
% WarmPlot = hist3([East_nonIT,North_nonIT],[20,20]);
% Plot5 = pcolor(WarmPlot);
% set(Plot5,'edgecolor','none')
% xlabel('East (m)')
% ylabel('North (m)')
% 
% % axes(handles.Plot6)
% % hold off
% % hist3([East_nonIT,North_nonIT],[20,20])
% % xlabel('East(m)')
% % ylabel('North (m)')
% % zlabel('Density')
% % grid on
% 
% axes(handles.Plot7)
% hold off
% [WarmPlot,xy] = hist3([East,North],[20,20]);
% Plot7 = pcolor(WarmPlot);
% xlim
% set(Plot7,'edgecolor','none')
% xlabel('East (m)')
% ylabel('North (m)')
% grid on
% colormap('jet')
% 
% % axes(handles.Plot8)
% % hold off
% % hist3([East_nonIT,North_nonIT],[20,20])
% % xlabel('East(m)')
% % ylabel('North (m)')
% % zlabel('Density')
% % grid on

function DisplayPlot(hObject,handles,PlotsToDisplay)

if strcmp(PlotsToDisplay,'1')
    set(handles.Plot1,'visible','on');
    set(handles.Plot2,'visible','off');
    cla(handles.Plot2);
    set(handles.Plot3,'visible','off');
    cla(handles.Plot3);
    set(handles.Plot4,'visible','off');
    cla(handles.Plot4);
    set(handles.Plot5,'visible','off');
    cla(handles.Plot5);
    set(handles.Plot6,'visible','off');
    cla(handles.Plot6);
    set(handles.Plot7,'visible','off');
    cla(handles.Plot7);
    set(handles.Plot8,'visible','off');
    cla(handles.Plot8);
elseif strcmp(PlotsToDisplay,'234')
    set(handles.Plot1,'visible','off');
    cla(handles.Plot1);
    set(handles.Plot2,'visible','on');
    set(handles.Plot3,'visible','on');
    set(handles.Plot4,'visible','on');
    set(handles.Plot5,'visible','off');
    cla(handles.Plot5);
    set(handles.Plot6,'visible','off');
    cla(handles.Plot6);
    set(handles.Plot7,'visible','off');
    cla(handles.Plot7);
    set(handles.Plot8,'visible','off');
    cla(handles.Plot8);
elseif strcmp(PlotsToDisplay,'5678') %????????????????????????????????????
    set(handles.Plot1,'visible','off');
    cla(handles.Plot1);
    set(handles.Plot2,'visible','off');
    cla(handles.Plot2);
    set(handles.Plot3,'visible','off');
    cla(handles.Plot3);
    set(handles.Plot4,'visible','off');
    cla(handles.Plot4);
    set(handles.Plot5,'visible','on');
    cla(handles.Plot5);
%     set(handles.Plot6,'visible','on');
    cla(handles.Plot6);
    set(handles.Plot7,'visible','on');
    cla(handles.Plot7);
%     set(handles.Plot8,'visible','on');
    cla(handles.Plot8);
end

set(handles.WTypeSelection,'Enable','on');

    if isfield(hObject,'string')
    if strcmp(get(hObject,'string'),'Orbits') || strcmp(get(hObject,'string'),'LAT - LON') || strcmp(get(hObject,'string'),'Position Estimates')
        set(handles.MinEpoch,'Enable','off')
        set(handles.MaxEpoch,'Enable','off')
        set(handles.SVSelection,'Enable','off')
    else
        set(handles.MinEpoch,'Enable','on')
        set(handles.MaxEpoch,'Enable','on')
        set(handles.SVSelection,'Enable','on')
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
WNLSE_Callback(hObject, eventdata, handles)
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
