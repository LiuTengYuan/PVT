function skyPlot(handles)

rad2deg=180/pi;
deg2rad=pi/180;

prn=handles.SVTracked;               %Satellite PRN
for index = 1 : handles.INDEX
    EpochToPlotIndex = find(handles.mTracked(:,handles.SVTracked(index)) ~= 0 );
    EpochToPlot(index) = EpochToPlotIndex(round(length(EpochToPlotIndex)/2));
    
    %% Satellite Information
    azi(index)=handles.azimuth_SV(EpochToPlot(index),index);  %Azimuth in degrees
    el(index)=handles.elevation_SV(EpochToPlot(index),index); %Elevation angle in degrees
end

%% Plot Figure
a=azi*deg2rad;                                  %Convert degrees to radians
r=90-el;                                        %Convert elevation angle to zenith

for i=1:size(azi,2), 
    svx(i)=r(i)*cos(a(i))  ; svy(i)=r(i)*sin(a(i)); %Calculate polar co-ordinates
end
polarhg([30 60])                                %Prerequisite script used to format axis
hold on
plot( svx,svy,'.','markers',50);               %Plot satellite location
hold off

%% Format output 
for i=1:length(prn),
    text(svx(i)+7,svy(i),strcat('SV',num2str(prn(i))), 'FontSize' ,10,'color','b','fontweight','bold') ; %Add PRN labels to each point
end

title('Sky Plot','fontweight','bold')
axis('square')
grid on;  
set(gcf, 'Color', 'w');                        %Change background of figure from grey to white
ti = get(gca,'TightInset')   ;                 %Remove extra spacing around figure
set(gca, 'LooseInset', [0,0,0,0.01]);          %Depending on the figure, you may need to add extra
                                               %spacing [left bottom width height])
% print( '-dtiff',  ['skyPlot'], '-r600');       %Change "-r600" to the required DPI

 