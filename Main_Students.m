clc;
clear;
close all;

%%-------------------------------------------------------------------------
%% Add path to the folders/subfolders that contain the useful functions

addpath(genpath('./Library'));

%%-------------------------------------------------------------------------
%% Define collected data files to process: .obs & .nav (Rinex 2.11)

filename_o = './Data/GYM/COM5_181004_125632.obs';
filename_n = './Data/GYM/COM5_181004_125632.nav';

% filename_o = './Data/Main Hall/COM6_181001_071919.obs';
% filename_n = './Data/Main Hall/COM6_181001_071919.nav';

% filename_o = './Data/test/COM44_150205_094639.obs';
% filename_n = './Data/test/COM44_150205_094639.nav';
ref_pos = [4627537.2739   119698.4035  4373317.5742];
ref_LLA = f_xyz_2_llh(ref_pos);

%%-------------------------------------------------------------------------
%% Initialize Variables

% Maximum number of epochs to process
nEpoch_max = 10000;
% Figure index
iFig = 0;
% Color grid
cmap=colormap(jet(32));
% ENAC reference position (ECEF) - [x,y,z]
ENAC_xyz = 1e6*[4.627536601003540,0.119700014080275,4.373318373560944];
% ENAC reference position (ECEF) - [latitude, longitude, heigth]
ENAC_llh = [43.564758116,1.48173363,203.8171];

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
[mEpoch, Nb_Epoch, vNb_Sat, Total_Nb_Sat, mTracked, mC1, mL1, mD1, mS1]=ExtractData_O(DATA_O, nEpoch_max);
fprintf('\nEnd of extracting the data.\n');

%There exist some epochs that can receive the satellite data yet cannot
%receive the ephemeris data
for epoch=1:Nb_Epoch
    redundant_SV = setdiff(find(mTracked(epoch,:)~=0),Ephem(:,1));
    %redundant_SV is the satellite that is in mTracked yet not in Ephem
    if redundant_SV
    vNb_Sat(epoch) = vNb_Sat(epoch)-length(setdiff(find(mTracked(epoch,:)~=0),Ephem(:,1)));
    mTracked(epoch,redundant_SV) = 0;
    mC1(epoch,redundant_SV) = 0;
    mL1(epoch,redundant_SV) = 0;
    mD1(epoch,redundant_SV) = 0;
    mS1(epoch,redundant_SV) = 0;
    end
end
Total_Nb_Sat = length(find(sum(mTracked)~=0));



%%-------------------------------------------------------------------------
%% Compute Transmission Time  &  SV Position and ClockCorrection

global c;
global f;
c = 299792458.0; %(m/s) speed of light
f = 1575.42e6; %(Hz) L1 frequency

Result(Nb_Epoch) = struct(); % This way automatically allocate memory for all Nb_Epoch structures.
Result_Info(Nb_Epoch) = struct();

Titles = {'SV # PRN','SV_X','SV_Y','SV_Z', 'SV_X_der','SV_Y_der','SV_Z_der','SV_ClockError_second','SV_ClockError_meter'};

for epoch=1:Nb_Epoch
    Result_Info(epoch).INFO = Titles;
    iUser_NoS = mEpoch(:,3); %user time(NoS)
    iUser_SoW = mEpoch(:,2); %user time(SoW)
    count = 1;
    for PRN=1:length(mTracked(1,:))
        if mTracked(epoch,PRN)==1 %Use mTracked(i,j) to decide if PRN j is tracked at epoch i
            f_C1 = mC1(epoch,PRN); %code range measurement of j satellite at epoch i
            %Select Ephermis (set the best fits SV iPRN at tiem iUser_Nos)
            [vEphemeris] = SelectEphemeris(Ephem, PRN, iUser_NoS(epoch));
            [iDelta_t_SV, iT_SV] = ComputeTransmissionTime(vEphemeris, iUser_SoW(epoch), f_C1);
            [rSatellite_xyz, vSatellite_xyz, fSV_ClockCorr] = SV_Position_and_ClockCorrection(vEphemeris, iT_SV, iUser_SoW(epoch), iDelta_t_SV);
            fSV_ClockCorr_meter = fSV_ClockCorr*c;
            Result(epoch).SV(count,:) =[PRN, rSatellite_xyz, vSatellite_xyz, fSV_ClockCorr, fSV_ClockCorr_meter];
            
            Result_Info(epoch).INFO = [Result_Info(epoch).INFO; num2cell([PRN, rSatellite_xyz, vSatellite_xyz,fSV_ClockCorr, fSV_ClockCorr_meter])];
            count = count+1;
        end
    end
    clc;
    ProcessingCompleted = epoch/(Nb_Epoch)*100;
    fprintf("\n SV Position Processig... ");
    fprintf(" \n Total Completed - %.2f %% \n",ProcessingCompleted);
end

%%-------------------------------------------------------------------------
%% Plot Satellite Orbit

SVTracked = find(sum(mTracked)~=0);

SV(Total_Nb_Sat) = struct();
index = 0;
for PRN=SVTracked
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

hold off
for index = 1 : length(SVTracked)
    EpochToPlotIndex = find(mTracked(:,SVTracked(index)) ~= 0 );
    EpochToPlot = round(( EpochToPlotIndex(1) + EpochToPlotIndex(end) ) / 2);
    scatter3(SV(index).Result_x(EpochToPlot)/1000,SV(index).Result_y(EpochToPlot)/1000,SV(index).Result_z(EpochToPlot)/1000,50,'d','filled','DisplayName',strcat('SV # ', num2str(SVTracked(index))));
    legend('-DynamicLegend')
    hold all
end
Legend = get(gca,'Legend');
Legend = Legend.String;
earth_sphere("km")
hold on
for index = 1 : length(SVTracked)
    plot3(SV(index).Result_x/1000,SV(index).Result_y/1000,SV(index).Result_z/1000,'k.');%,'DisplayName',sprintf("PRN %d", PRN));
    hold on
end
grid on
title('Satellites orbits during data collection','fontweight','bold')
legend(Legend);



%%-------------------------------------------------------------------------
%% Receiver Position and Receiver Clock Error
Tiono = zeros(Nb_Epoch,Total_Nb_Sat);
Ttropo = zeros(Nb_Epoch,Total_Nb_Sat);
[RX_Position_XYZ, RX_ClockError, Matrix] = RX_Position_and_Clock(Result,mC1,mD1,mS1,Nb_Epoch,vNb_Sat,'NWLSE',0,Tiono,Ttropo,2,[],[],SVTracked);
[RX_Position_LLH, RX_Position_ENU, Matrix, DOP] = RX_Position_LLH_ENU(RX_Position_XYZ,Nb_Epoch,Matrix);

mean_LLH = mean(RX_Position_LLH);
stdev_LLH = std(RX_Position_LLH);
mean_ENU = mean(RX_Position_ENU);
stdev_ENU = std(RX_Position_ENU);

for epoch = 1 : Nb_Epoch
    for index = 1 : length(SVTracked)
        [SV(index).llh(epoch,:)] = xyz_2_lla_PVT( [SV(index).Result_x(epoch), SV(index).Result_y(epoch), SV(index).Result_z(epoch)] );
        SV(index).llh(epoch,1) = rad2deg(SV(index).llh(epoch,1));
        SV(index).llh(epoch,2) = rad2deg(SV(index).llh(epoch,2));
    end
end

%%-------------------------------------------------------------------------
%% plot error ellispe ???????????????????????????????????????????????????

Sigma_URE = 6;
RMS_Errors = Sigma_URE * [sqrt(DOP.EDOP.^2 + DOP.NDOP.^2); ...
    sqrt(DOP.EDOP.^2 + DOP.NDOP.^2 + DOP.VDOP.^2);...
    sqrt(DOP.EDOP.^2 + DOP.NDOP.^2 + DOP.VDOP.^2 + DOP.TDOP.^2); ...
    DOP.TDOP];

% for epoch=1:Nb_Epoch
%     PCOV = Matrix(epoch).COV(1:2,1:2);
%     P_deltaz = Sigma_URE*eye(2);
%     P_deltap = PCOV*P_deltaz;
%     [eigenvector,eigenvalue] = eig(inv(P_deltap));
%     ellipse_angle = atan2(eigenvector(1,2),eigenvector(1,1));
%     H_ellipse = ellipse(eigenvalue(1,1),eigenvalue(2,2),ellipse_angle...
%         ,RX_Position_ENU(epoch,1),RX_Position_ENU(epoch,2));
%     hold on
%     plot(RX_Position_ENU(epoch,1),RX_Position_ENU(epoch,2),'.','Markersize',5)
%     hold on
% end

%%-------------------------------------------------------------------------
%% Calculate and Plot Elevation & Azimuth

clear Elevation_Azimuth;
Elevation_Azimuth(Nb_Epoch) = struct();
for epoch=1:Nb_Epoch
    Elevation_Azimuth(epoch).SV(:,1) = Result(epoch).SV(:,1);
    for epoch_sv=1:vNb_Sat(epoch)
        [fElevation, fAzimuth] = elevation_azimuth(RX_Position_XYZ(epoch,:), Result(epoch).SV(epoch_sv,2:4));
        Elevation_Azimuth(epoch).SV(epoch_sv,2) = rad2deg(fElevation);
        Elevation_Azimuth(epoch).SV(epoch_sv,3) = rad2deg(fAzimuth);
    end
end

azimuth_SV = zeros(Nb_Epoch,length(SVTracked));
elevation_SV = zeros(Nb_Epoch,length(SVTracked));
for epoch=1:Nb_Epoch
    for num_SV=1:length(SVTracked)
        find_SV = find(Elevation_Azimuth(epoch).SV(:,1)==SVTracked(num_SV));
        if find_SV
            elevation_SV(epoch,num_SV) = Elevation_Azimuth(epoch).SV(find_SV,2);
            azimuth_SV(epoch,num_SV) = Elevation_Azimuth(epoch).SV(find_SV,3);
        end
    end
end
elevation_SV(elevation_SV==0) = nan;
azimuth_SV(azimuth_SV==0) = nan;

% figure;
% plot(elevation_SV,'linewidth',2)
% ylabel('m')
% xlabel('Epoch')
% title('Elevation between RX and SV','fontweight','bold')
% legend(strcat('PRN # ', string(SVTracked)),'Location','BestOutside')
%
% figure;
% plot(azimuth_SV,'linewidth',2)
% ylabel('m')
% xlabel('Epoch')
% title('Azimuth between RX and SV','fontweight','bold')
% legend(strcat('PRN # ', string(SVTracked)),'Location','BestOutside')

%%-------------------------------------------------------------------------
%% Ionosphere Correction & Troposhere Correction

[Tiono] = Ionospheric_Correction(Iono_a, Iono_b, Elevation_Azimuth, RX_Position_LLH, mEpoch, Total_Nb_Sat);
[Ttropo] = Tropospheric_Correction(RX_Position_LLH, mEpoch, Elevation_Azimuth, Total_Nb_Sat);
[RX_Position_XYZ_IT, RX_ClockError_IT, Matrix_IT] = RX_Position_and_Clock(Result,mC1,mS1,Nb_Epoch,vNb_Sat,'NWLSE',[],Tiono,Ttropo,2,[],[],SVTracked);
[RX_Position_LLH_IT, RX_Position_ENU_IT, Matrix_IT, DOP_IT] = RX_Position_LLH_ENU(RX_Position_XYZ_IT,Nb_Epoch,Matrix_IT);

%%-------------------------------------------------------------------------
%% Receiver Velocity and Receiver Clock Error Rate ???????????????????????

% new Rseult & new Result_Info (SV Velocity and Clock Correction Rate)
[Result, Result_Info] = SV_Velocity_and_ClockCorrectionRate(Result,Result_Info,mTracked,SVTracked);
[RX_Velocity_XYZ, RX_ClockErrorRate] = RX_Velocity_and_Clock_Rate(Result,mD1,Nb_Epoch,vNb_Sat,Matrix_IT);


