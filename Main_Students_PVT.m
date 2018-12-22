function [handles] = Main_Students_PVT(handles)

clc;
% clear;
% close all;

%%-------------------------------------------------------------------------
%% Add path to the folders/subfolders that contain the useful functions

addpath(genpath('./Library'));

%%-------------------------------------------------------------------------
%% Define collected data files to process: .obs & .nav (Rinex 2.11)

filename_o = handles.path_filename_o; %'./Data/test/COM44_150205_094639.obs';
filename_n = handles.path_filename_n; %'./Data/test/COM44_150205_094639.nav';
ref_pos = handles.ref_pos; %[4627537.2739   119698.4035  4373317.5742];
ref_LLA = f_xyz_2_llh(ref_pos);

%%-------------------------------------------------------------------------
%% Initialize Variables

% Maximum number of epochs to process
nEpoch_max = handles.nEpoch_max; %10000;
% Figure index
iFig = 0;
% Color grid
cmap=colormap(jet(32));
% ENAC reference position (ECEF) - [x,y,z]
ENAC_xyz = handles.ENAC_xyz; %1e6*[4.627536601003540,0.119700014080275,4.373318373560944];
% ENAC reference position (ECEF) - [latitude, longitude, heigth]
ENAC_llh = handles.ENAC_llh; %[43.564758116,1.48173363,203.8171];

%%-------------------------------------------------------------------------
%% Read RINEX files

% Read observation file (.obs)
[OUTPUT io_flag_o] = read_RINEX_OBS_v2(filename_o);
if (io_flag_o ~= 0)
    fprintf(1,'Error reading observation file\n');
    return;
end
HEADER_O = OUTPUT.HEADER; DATA_O = OUTPUT.DATA;
% Read navigation file (.nav)
[OUTPUT io_flag_n] = read_RINEX_NAV_v2(filename_n);
if (io_flag_n ~= 0)
    fprintf(1,'Error reading navigation file\n');
    return;
end
HEADER_N = OUTPUT.HEADER; DATA_N = OUTPUT.DATA;
fprintf('\nEnd of reading the RINEX files.\n');

%%-------------------------------------------------------------------------
%% Extract data (see functions header for variables details)

% Extract GPS navigation message data
[Iono_a, Iono_b, Ephem] = ExtractData_N(HEADER_N, DATA_N);
% Extract GPS observation data
[mEpoch, Nb_Epoch, vNb_Sat, Total_Nb_Sat, mTracked, handles.mC1, handles.mL1, handles.mD1, handles.mS1]=ExtractData_O(DATA_O, nEpoch_max);
fprintf('\nEnd of extracting the data.\n');

%There exist some epochs that can receive the satellite data yet cannot
%receive the ephemeris data
for epoch=1:Nb_Epoch
    redundant_SV = setdiff(find(mTracked(epoch,:)~=0),Ephem(:,1));
    %redundant_SV is the satellite that is in mTracked yet not in Ephem
    if redundant_SV
        vNb_Sat(epoch) = vNb_Sat(epoch)-length(setdiff(find(mTracked(epoch,:)~=0),Ephem(:,1)));
        mTracked(epoch,redundant_SV) = 0;
        handles.mC1(epoch,redundant_SV) = 0;
        handles.mL1(epoch,redundant_SV) = 0;
        handles.mD1(epoch,redundant_SV) = 0;
        handles.mS1(epoch,redundant_SV) = 0;
    end
end
Total_Nb_Sat = length(find(sum(mTracked)~=0));

handles.Nb_Epoch = Nb_Epoch;
handles.vNb_Sat = vNb_Sat;

%%-------------------------------------------------------------------------
%% Pseudorange Model

switch handles.PseudorangeModel
    case 'Code'
        handles.mC1 = handles.mC1;
    case 'CodeAndCarrier'
        %         mC1 = CODE + CARREIR;
end
%%-------------------------------------------------------------------------

%% Avoid SVs tracked
if handles.SVListFilter ~= 0
    mTracked(:,handles.SVListFilter) = zeros(Nb_Epoch,length(handles.SVListFilter));
    handles.mC1(:,handles.SVListFilter) = zeros(Nb_Epoch,length(handles.SVListFilter));
    handles.mD1(:,handles.SVListFilter) = zeros(Nb_Epoch,length(handles.SVListFilter));
    handles.mS1(:,handles.SVListFilter) = zeros(Nb_Epoch,length(handles.SVListFilter));
    handles.mL1(:,handles.SVListFilter) = zeros(Nb_Epoch,length(handles.SVListFilter));
    for epoch = 1 : Nb_Epoch
        vNb_Sat(epoch) = sum(mTracked(epoch,:));
    end
end


%%-------------------------------------------------------------------------
%% Compute Transmission Time  &  SV Position and ClockCorrection

global c;
global f;
global lambda;
c = 299792458.0; %(m/s) speed of light
f = 1575.42e6; %(Hz) L1 frequency
lambda = c/f; %(m) L1 wavelength
Result(Nb_Epoch) = struct(); % This way automatically allocate memory for all Nb_Epoch structures.
Result_Info(Nb_Epoch) = struct();

Titles = {'SV # PRN','SV_X','SV_Y','SV_Z','SV_ClockError_second','SV_ClockError_meter'};


for epoch=1:Nb_Epoch
    Result_Info(epoch).INFO = Titles;
    iUser_NoS = mEpoch(:,3); %user time(NoS)
    iUser_SoW = mEpoch(:,2); %user time(SoW)
    count = 1;
    for PRN=1:length(mTracked(1,:))
        if mTracked(epoch,PRN)==1 %Use mTracked(i,j) to decide if PRN j is tracked at epoch i
            f_C1 = handles.mC1(epoch,PRN); %code range measurement of j satellite at epoch i
            %Select Ephermis (set the best fits SV iPRN at tiem iUser_Nos)
            [vEphemeris] = SelectEphemeris(Ephem, PRN, iUser_NoS(epoch));
            [iDelta_t_SV, iT_SV] = ComputeTransmissionTime(vEphemeris, iUser_SoW(epoch), f_C1);
            [vSatellite_xyz, fSV_ClockCorr] = SV_Position_and_ClockCorrection(vEphemeris, iT_SV, iUser_SoW(epoch), iDelta_t_SV);
            fSV_ClockCorr_meter = fSV_ClockCorr*c;
            Result(epoch).SV(count,:) =[PRN, vSatellite_xyz, fSV_ClockCorr, fSV_ClockCorr_meter];
            Result_Info(epoch).INFO = [Result_Info(epoch).INFO; num2cell([PRN, vSatellite_xyz,fSV_ClockCorr, fSV_ClockCorr_meter])];
            count = count+1;
            
        end
    end
    
    clc;
    ProcessingCompleted = round(epoch/(Nb_Epoch)*100);
    fprintf("\n SV Position Processig... ");
    fprintf(" \n Total Completed - %.2f %% \n",ProcessingCompleted);
    if mod(ProcessingCompleted,10)==0
        set(handles.MessageBox,'string',sprintf("\n SV Position Processig...\n Total Completed - %.2f %%",ProcessingCompleted));
        pause(0.01)
    end
end

%%-------------------------------------------------------------------------
%% Plot Satellite Orbit

handles.mTracked = mTracked;
handles.SVTracked = find(sum(handles.mTracked)~=0);

SV(Total_Nb_Sat) = struct();
index = 0;
for PRN=handles.SVTracked
    index = index + 1;
    SV(index).PRN = PRN;
    for epoch=1:Nb_Epoch
        INDEX = find(Result(epoch).SV(:,1) == PRN);
        if mTracked(epoch,PRN)
            SV(index).Result_x(epoch) = Result(epoch).SV(INDEX,2); % We have to plot a vector once directly, instead of ploting point by point. Thus, we need to pass from the struct to a vector.
            SV(index).Result_y(epoch) = Result(epoch).SV(INDEX,3); % In every epoch we have 1 position for each satellite, this forces us to plot point by point and may take several time (Ni plots).
            SV(index).Result_z(epoch) = Result(epoch).SV(INDEX,4); % So we change rewrite as Ni positions for each satellite, going around all its epochs.
        else
            SV(index).Result_x(epoch) = NaN; % We have to plot a vector once directly, instead of ploting point by point. Thus, we need to pass from the struct to a vector.
            SV(index).Result_y(epoch) = NaN; % In every epoch we have 1 position for each satellite, this forces us to plot point by point and may take several time (Ni plots).
            SV(index).Result_z(epoch) = NaN;
        end
    end
    
end

handles.Result = Result;
handles.SV = SV;
handles.INDEX = length(handles.SVTracked);

axes(handles.Plot1)
hold off
for index = 1 : length(handles.SVTracked)
    EpochToPlotIndex = find(mTracked(:,handles.SVTracked(index)) ~= 0 );
    EpochToPlot = round(( EpochToPlotIndex(1) + EpochToPlotIndex(end) ) / 2);
    scatter3(SV(index).Result_x(EpochToPlot)/1000,SV(index).Result_y(EpochToPlot)/1000,SV(index).Result_z(EpochToPlot)/1000,50,'d','filled','DisplayName',strcat('SV # ', num2str(handles.SVTracked(index))));
    legend('-DynamicLegend')
    hold all
end
Legend = get(gca,'Legend');
Legend = Legend.String;
earth_sphere("km")
hold on
for index = 1 : length(handles.SVTracked)
    plot3(SV(index).Result_x/1000,SV(index).Result_y/1000,SV(index).Result_z/1000,'k.');%,'DisplayName',sprintf("PRN %d", PRN));
    hold on
end
grid on
title('Satellites orbits during data collection','fontweight','bold')
legend(Legend);


%%-------------------------------------------------------------------------
%% Recaculate Tracked_mS1,Tracked_mC1,Tracked_mL1

handles.Tracked_mS1 = zeros(Nb_Epoch,length(handles.SVTracked));
for num_SV=1:length(handles.SVTracked)
    handles.Tracked_mS1(:,num_SV) = handles.mS1(:,handles.SVTracked(num_SV));
end
handles.Tracked_mS1(handles.Tracked_mS1==0) = nan;

handles.Tracked_mC1 = zeros(Nb_Epoch,length(handles.SVTracked));
for num_SV=1:length(handles.SVTracked)
    handles.Tracked_mC1(:,num_SV) = handles.mC1(:,handles.SVTracked(num_SV));
end
handles.Tracked_mC1(handles.Tracked_mC1==0) = nan;

handles.Tracked_mL1 = zeros(Nb_Epoch,length(handles.SVTracked));
for num_SV=1:length(handles.SVTracked)
    handles.Tracked_mL1(:,num_SV) = handles.mL1(:,handles.SVTracked(num_SV));
end
handles.Tracked_mL1(handles.Tracked_mL1==0) = nan;

%%-------------------------------------------------------------------------

%% Caculate Receiver Position and Receiver Clock Error

Tiono = zeros(Nb_Epoch,Total_Nb_Sat);
Ttropo = zeros(Nb_Epoch,Total_Nb_Sat);

% Non Linear LSE - NO ATMOSPHERIC CORRECTION
[handles.RX_Position_XYZ_NLSE, handles.RX_ClockError_NLSE, handles.Matrix_NLSE] = RX_Position_and_Clock(Result,handles.mC1,handles.mS1,Nb_Epoch,vNb_Sat,'NLSE',0,Tiono,Ttropo,1,[],handles.MessageBox);
[handles.RX_Position_LLH_NLSE, handles.RX_Position_ENU_NLSE, handles.Matrix_NLSE, handles.DOP_NLSE] = RX_Position_LLH_ENU(handles.RX_Position_XYZ_NLSE,Nb_Epoch,handles.Matrix_NLSE);

% SV Latitude, Longitude and Height.

for epoch = 1 : Nb_Epoch
    for index = 1 : length(handles.SVTracked)
        [handles.SV(index).llh(epoch,:)] = xyz_2_lla_PVT( [handles.SV(index).Result_x(epoch), handles.SV(index).Result_y(epoch), handles.SV(index).Result_z(epoch)] );
        handles.SV(index).llh(epoch,1) = rad2deg(handles.SV(index).llh(epoch,1));
        handles.SV(index).llh(epoch,2) = rad2deg(handles.SV(index).llh(epoch,2));
    end
end

%% Calculate and Plot Elevation & Azimuth

Elevation_Azimuth(Nb_Epoch) = struct();
for epoch=1:Nb_Epoch
    Elevation_Azimuth(epoch).SV(:,1) = Result(epoch).SV(:,1);
    for epoch_sv=1:vNb_Sat(epoch)
        [fElevation, fAzimuth] = elevation_azimuth(handles.RX_Position_XYZ_NLSE(epoch,:), Result(epoch).SV(epoch_sv,2:4));
        Elevation_Azimuth(epoch).SV(epoch_sv,2) = rad2deg(fElevation);
        Elevation_Azimuth(epoch).SV(epoch_sv,3) = rad2deg(fAzimuth);
    end
end

azimuth_SV = zeros(Nb_Epoch,length(handles.SVTracked));
elevation_SV = zeros(Nb_Epoch,length(handles.SVTracked));
for epoch=1:Nb_Epoch
    for num_SV=1:length(handles.SVTracked)
        find_SV = find(Elevation_Azimuth(epoch).SV(:,1)==handles.SVTracked(num_SV));
        if find_SV
            elevation_SV(epoch,num_SV) = Elevation_Azimuth(epoch).SV(find_SV,2);
            azimuth_SV(epoch,num_SV) = Elevation_Azimuth(epoch).SV(find_SV,3);
        end
    end
end
elevation_SV(elevation_SV==0) = nan;
azimuth_SV(azimuth_SV==0) = nan;

handles.Elevation_Azimuth = Elevation_Azimuth;
handles.elevation_SV = elevation_SV;
handles.azimuth_SV = azimuth_SV;


% Weighted (SNR) Non Linear LSE - NO ATMOSPHERIC CORRECTION
[handles.RX_Position_XYZ_W(1).NLSE, handles.RX_ClockError_W(1).NLSE, handles.Matrix_W(1).NLSE] = RX_Position_and_Clock(Result,handles.mC1,handles.mS1,Nb_Epoch,vNb_Sat,'NWLSE',1,Tiono,Ttropo,2,[],handles.MessageBox);
[handles.RX_Position_LLH_W(1).NLSE, handles.RX_Position_ENU_W(1).NLSE, handles.Matrix_W(1).NLSE, handles.DOP_W(1).NLSE] = RX_Position_LLH_ENU(handles.RX_Position_XYZ_W(1).NLSE,Nb_Epoch,handles.Matrix_W(1).NLSE);
%%-------------------------------------------------------------------------

% Weighted (SNR + ELEVATION) Non Linear LSE - NO ATMOSPHERIC CORRECTION
[handles.RX_Position_XYZ_W(2).NLSE, handles.RX_ClockError_W(2).NLSE, handles.Matrix_W(2).NLSE] = RX_Position_and_Clock(Result,handles.mC1,handles.mS1,Nb_Epoch,vNb_Sat,'NWLSE',2,Tiono,Ttropo,3,Elevation_Azimuth,handles.MessageBox);
[handles.RX_Position_LLH_W(2).NLSE, handles.RX_Position_ENU_W(2).NLSE, handles.Matrix_W(2).NLSE, handles.DOP_W(2).NLSE] = RX_Position_LLH_ENU(handles.RX_Position_XYZ_W(2).NLSE,Nb_Epoch,handles.Matrix_W(2).NLSE);
%%-------------------------------------------------------------------------



%%-------------------------------------------------------------------------
%% Ionosphere Correction & Troposhere Correction

[Tiono] = Ionospheric_Correction(Iono_a, Iono_b, Elevation_Azimuth, handles.RX_Position_LLH_NLSE, mEpoch, Total_Nb_Sat);
[Ttropo] = Tropospheric_Correction(handles.RX_Position_LLH_NLSE, mEpoch, Elevation_Azimuth, Total_Nb_Sat);

handles.Tiono = Tiono;
handles.Ttropo = Ttropo;

%Non Linear LSE
[handles.RX_Position_XYZ_NLSE_IT, handles.RX_ClockError_NLSE_IT, handles.Matrix_NLSE_IT] = RX_Position_and_Clock(Result,handles.mC1,handles.mS1,Nb_Epoch,vNb_Sat,'NLSE',0,Tiono,Ttropo,4,[],handles.MessageBox);
[handles.RX_Position_LLH_NLSE_IT, handles.RX_Position_ENU_NLSE_IT, handles.Matrix_NLSE_IT, handles.DOP_NLSE_IT] = RX_Position_LLH_ENU(handles.RX_Position_XYZ_NLSE_IT,Nb_Epoch,handles.Matrix_NLSE_IT);

% Weighted (SNR) Non Linear LSE - ATMOSPHERIC CORRECTION
[handles.RX_Position_XYZ_W(1).NLSE_IT, handles.RX_ClockError_W(1).NLSE, handles.Matrix_W(1).NLSE_IT] = RX_Position_and_Clock(Result,handles.mC1,handles.mS1,Nb_Epoch,vNb_Sat,'NWLSE',1,Tiono,Ttropo,5,[],handles.MessageBox);
[handles.RX_Position_LLH_W(1).NLSE_IT, handles.RX_Position_ENU_W(1).NLSE_IT, handles.Matrix_W(1).NLSE_IT, handles.DOP_W(1).NLSE_IT] = RX_Position_LLH_ENU(handles.RX_Position_XYZ_W(1).NLSE_IT,Nb_Epoch,handles.Matrix_W(1).NLSE_IT);

% Weighted (SNR + ELEVATION) Non Linear LSE - ATMOSPHERIC CORRECTION
[handles.RX_Position_XYZ_W(2).NLSE_IT, handles.RX_ClockError_W(2).NLSE_IT, handles.Matrix_W(2).NLSE_IT] = RX_Position_and_Clock(Result,handles.mC1,handles.mS1,Nb_Epoch,vNb_Sat,'NWLSE',2,Tiono,Ttropo,6,Elevation_Azimuth,handles.MessageBox);
[handles.RX_Position_LLH_W(2).NLSE_IT, handles.RX_Position_ENU_W(2).NLSE_IT, handles.Matrix_W(2).NLSE_IT, handles.DOP_W(2).NLSE_IT] = RX_Position_LLH_ENU(handles.RX_Position_XYZ_W(2).NLSE_IT,Nb_Epoch,handles.Matrix_W(2).NLSE_IT);
%%-------------------------------------------------------------------------

%%-------------------------------------------------------------------------

end
