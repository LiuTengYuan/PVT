function [Ttropo] = Tropspheric_Correction(RX_Position_LLH, mEpoch, Elevation_Azimuth, Total_Nb_Sat)

Nb_Epoch = length(mEpoch);
Ttropo = zeros(Nb_Epoch,Total_Nb_Sat);
LATRAD = deg2rad(RX_Position_LLH(:,1));
HEIGHTM = RX_Position_LLH(:,3);
DAYOYEAR = mEpoch(:,5)*30+mEpoch(:,6);
ELEVRAD = zeros(Nb_Epoch,Total_Nb_Sat);
SVnum = zeros(Nb_Epoch,1);
for epoch=1:Nb_Epoch
    SVnum(epoch) = length(Elevation_Azimuth(epoch).SV(:,1));
    ELEVRAD(epoch,1:SVnum(epoch)) = deg2rad(Elevation_Azimuth(epoch).SV(:,2)');
end

RTROP = zeros(Nb_Epoch,Total_Nb_Sat);
DRATE = zeros(Nb_Epoch,Total_Nb_Sat);
T = zeros(Nb_Epoch,Total_Nb_Sat);
P = zeros(Nb_Epoch,Total_Nb_Sat);
E = zeros(Nb_Epoch,Total_Nb_Sat);
TM = zeros(Nb_Epoch,Total_Nb_Sat);
for epoch=1:Nb_Epoch
    for num_SV = 1:SVnum(epoch)
        [RTROP(epoch,num_SV) HZD HMF WZD WMF]=UNB3M(LATRAD(epoch),HEIGHTM(epoch),DAYOYEAR(epoch),ELEVRAD(epoch,num_SV));
        [DRATE(epoch,num_SV) HZD DHMFDEL WZD DWMFDEL]=UNB3MR(LATRAD(epoch),HEIGHTM(epoch),DAYOYEAR(epoch),ELEVRAD(epoch,num_SV));
        [T(epoch,num_SV) P(epoch,num_SV) E(epoch,num_SV) TM(epoch,num_SV)]=UNB3MM(LATRAD(epoch),HEIGHTM(epoch),DAYOYEAR(epoch));
    end
end
Ttropo = RTROP;

end
