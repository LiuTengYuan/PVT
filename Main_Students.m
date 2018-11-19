

clc;
clear;
close all;

%%-------------------------------------------------------------------------
%% Add path to the folders/subfolders that contain the useful functions

addpath(genpath('./Library'));

%%-------------------------------------------------------------------------
%% Define collected data files to process: .obs & .nav (Rinex 2.11)

filename_o = './Data/test/COM44_150205_094639.obs';
filename_n = './Data/test/COM44_150205_094639.nav';
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

%%-------------------------------------------------------------------------
%% Compute Transmission Time  &  SV Position and ClockCorrection

global c;
Result(Nb_Epoch) = struct(); % This way automatically allocate memory for all Nb_Epoch structures.
Result_Info(Nb_Epoch) = struct();
Epoch_SV_Number = zeros(Nb_Epoch,1); %the number of satellite for every epoch

Titles = {'SV # PRN','SV_X','SV_Y','SV_Z','SV_ClockError_second','SV_ClockError_meter'};

for epoch=1:Nb_Epoch
    Epoch_SV_Number(epoch) = sum(mTracked(epoch,:));
    Result_Info(epoch).INFO = Titles;
    iUser_NoS = mEpoch(:,3); %user time(NoS)
    iUser_SoW = mEpoch(:,2); %user time(SoW)
    count = 1;
    for PRN=1:32
        if mTracked(epoch,PRN)==1 %Use mTracked(i,j) to decide if PRN j is tracked at epoch i
            f_C1 = mC1(epoch,PRN); %code range measurement of j satellite at epoch i
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
    ProcessingCompleted = epoch/(Nb_Epoch)*100;
    fprintf("\n SV Position Processig... ");
    fprintf(" \n Total Completed - %.2f %% \n",ProcessingCompleted);
end

%%-------------------------------------------------------------------------
%% Plot Satellite Orbit

SV = struct();
SVTracked = zeros();
INDEX = 0;
LastEpochCheck = 1;
for PRN=1:32
    ChangeSV = 1;
    for epoch=1:Nb_Epoch
        if mTracked(epoch,PRN)
            if ChangeSV, INDEX = INDEX + 1; FirstEpoch = epoch; ChangeSV = 0; end
            
            findprn = find(Result(epoch).SV(:,1)==PRN);
            SV(INDEX).Result_x(epoch - FirstEpoch + 1) = Result(epoch).SV(findprn,2);
            SV(INDEX).Result_y(epoch - FirstEpoch + 1) = Result(epoch).SV(findprn,3);
            SV(INDEX).Result_z(epoch - FirstEpoch + 1) = Result(epoch).SV(findprn,4);
            % We have to plot a vector once directly, instead of ploting point by point. Thus, we need to pass from the struct to a vector.
            % In every epoch we have 1 position for each satellite, this forces us to plot point by point and may take several time (Ni plots).
            % So we change rewrite as Ni positions for each satellite, going around all its epochs.
            SVTracked(INDEX) = PRN;
        end
    end
    
end
for index = 1 : INDEX
    EpochToPlot = round(length(SV(index).Result_x)/2);
    scatter3(SV(index).Result_x(EpochToPlot)/1000,SV(index).Result_y(EpochToPlot)/1000,SV(index).Result_z(EpochToPlot)/1000,50,'d','filled','DisplayName',strcat('SV # ', num2str(SVTracked(index))));
    legend('-DynamicLegend')
    hold all
end
Legend = get(gca,'Legend');
Legend = Legend.String;
earth_sphere("km")
hold on
for index = 1 : INDEX
    plot3(SV(index).Result_x/1000,SV(index).Result_y/1000,SV(index).Result_z/1000,'k.');
    %,'DisplayName',sprintf("PRN %d", PRN));
    hold on
end
grid on
legend(Legend);

figure;
plot(mS1);

%%-------------------------------------------------------------------------
%% Caculate Receiver Position and Receiver Clock Error

[RX_Position_XYZ, RX_ClockError, Matrix] = RX_Position_and_Clock(Result,mC1,Nb_Epoch,Epoch_SV_Number,'NLSE');
[RX_Position_LLH, RX_Position_ENU] = RX_Position_LLH_ENU(RX_Position_XYZ,Nb_Epoch);

% figure;
% plot3(RX_Position_XYZ(:,1),RX_Position_XYZ(:,2),RX_Position_XYZ(:,3),'r.');
% grid on;
% figure;
% plot(iUser_NoS,RX_ClockError);
% figure;
% plot(iUser_NoS,RX_Position_XYZ(:,1));
% figure;
% plot(iUser_NoS,RX_Position_XYZ(:,2));
% figure;
% plot(iUser_NoS,RX_Position_XYZ(:,3));

% figure;
% plot(iUser_NoS,RX_Position_LLH(:,1));
% figure;
% plot(iUser_NoS,RX_Position_LLH(:,2));
% figure;
% plot(iUser_NoS,RX_Position_LLH(:,3));

mean_LLH = mean(RX_Position_LLH);
stdev_LLH = std(RX_Position_LLH);
mean_ENU = mean(RX_Position_ENU);
stdev_ENU = std(RX_Position_ENU);
% figure;
% plot(RX_Position_ENU(:,1));
% figure;
% plot(RX_Position_ENU(:,2));
% figure;
% plot(RX_Position_ENU(:,3));


%%-------------------------------------------------------------------------
%% Calculate and Plot Elevation & Azimuth

% clear Elevation_Azimuth;
Elevation_Azimuth(Nb_Epoch) = struct();
for epoch=1:Nb_Epoch
    Elevation_Azimuth(epoch).SV(:,1) = Result(epoch).SV(:,1);
    for epoch_sv=1:Epoch_SV_Number(epoch)
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

figure;
plot(elevation_SV,'linewidth',2)
ylabel('m')
xlabel('Epoch')
title('Elevation between RX and SV','fontweight','bold')
legend(strcat('PRN # ', string(SVTracked)),'Location','BestOutside')

figure;
plot(azimuth_SV,'linewidth',2)
ylabel('m')
xlabel('Epoch')
title('Azimuth between RX and SV','fontweight','bold')
legend(strcat('PRN # ', string(SVTracked)),'Location','BestOutside')




