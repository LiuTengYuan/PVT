function [vSatellite_xyz, fSV_ClockCorr] = SV_Position_and_ClockCorrection(vEphemeris, iT_SV, iUser_SoW, iDelta_t_SV)
%--------------------------------------------------------------------------
% Copyright © ENAC, 2015.
% ENAC : http://www.enac.fr/.
% signav@recherche.enac.fr
%
% This function computes SV PRN Earth-fixed position and SV PRN clock
% correction (in meter), from ephemeris data.
% Input Variables:
% 1) vEphemeris satellite ephemeris data - 1x29 vector
% 2) iT_SV signal transmission time in GPS system time [s]
% 3) iUser_SoW user time (Second of Week - SoW)
% 4) iDelta_t_SV SV PRN code phase offset [s] (do not include the
% relativistic correction term)
% Output Variables:
% 1) vSatellite_xyz satellite position in rectangular coordinates in
% ECEF - 1x4 vector
% vSatellite_xyz(1) x coordinate
% vSatellite_xyz(2) y coordinate
% vSatellite_xyz(3) z coordinate
% 2) fSV_ClockCorr SV PRN clock correction [m]
%
% Reference:
% ICD-GPS-200C, 10 OCT 1993, page 98
%--------------------------------------------------------------------------

%Initialize Variables
global c;
OmegaEDot = 7.2921151467e-5; %(r/s)
mu = 3.986005e+14; %(m^3/s^2)
F = -4.442807633e-10; %(s/m^0.5)

toe = vEphemeris(13); %Referenece Time Ephemeris
e = vEphemeris(15); %Eccentricity
sqrtA = vEphemeris(16); %Square Root of the Semi-Major Axis
Omega0 = vEphemeris(17); %Longitude of ascending node of orbital plane at weekly epoch
i0 = vEphemeris(18); %Inclination angle at reference time
IDOT = vEphemeris(19); % Rate of Inclination Angle
omega = vEphemeris(20); %Argument of perigee
OmegaDot = vEphemeris(21); %Rate of right ascension
M0 = vEphemeris(22); %Mean Anomaly at Reference Time
deltaN = vEphemeris(23); %Mean Motion Difference From Computed Value
Crs = vEphemeris(24); %Amplitude of the Sine Harmonic Correction Term to the Arguemnt of Latitude
Crc = vEphemeris(25); %Amplitude of the Cosine Harmonic Correction Term to the Arguemnt of Latitude
Cus = vEphemeris(26); %Amplitude of the Sine Harmonic Correction Term to the Orbit Radius
Cuc = vEphemeris(27); %Amplitude of the Cosine Harmonic Correction Term to the Orbit Radius
Cis = vEphemeris(28); %Amplitude of the Sine Harmonic Correction Term to the Angle of Inclination
Cic = vEphemeris(29); %Amplitude of the Cosine Harmonic Correction Term to the Angle of Inclination

%Caculate Satellite Position
A = sqrtA^2; %Semi-Major Axis
n0 = sqrt(mu/A^3); %Computed Mean Motion(rad/sec)
delta_tr = 0; %Relativistic Correction Term

while(1)
    iT_SV_corrected = iT_SV-delta_tr; %GPS Time at Time of Transmission
    tk = iT_SV_corrected-toe; %Time from Ephemeris Reference Epoch
    if tk>=302400
        tk = tk-604800;
    elseif tk<-302400
            tk = tk+604800;
    end
    n = n0+deltaN; %Corrected Mean Motion
    Mk = M0+n*tk; %Mean Anomaly
    E0 = Mk;
    Ek = Mk+e*sin(E0); %Eccentric Anomaly
    
    %Iteration for Kepler's Equation
    while abs(Ek-E0)>=(10^-12)
        E0 = Ek;
        Ek = Mk+e*sin(E0);
    end
    
    delta_tr_new = F*(e*(sqrtA))*sin(Ek);
    fSV_ClockCorr = iDelta_t_SV+delta_tr;
    %Iteration for delta_tr
    if abs(delta_tr_new-delta_tr)<=(10^-12)
        break;
    end
    delta_tr = delta_tr_new;
end
vk = atan2(sqrt(1-e^2)*sin(Ek),(cos(Ek)-e)); %True Anomaly
Ek = acos((e+cos(vk))/(1+e*cos(vk))); %Eccentric Anomaly
thetak = vk+omega; %Argument of Latitude

%Second Harmonic Perturbations
delta_uk = Cus*sin(2*thetak)+Cuc*cos(2*thetak); %Argument of Latitude Correction
delta_rk = Crs*sin(2*thetak)+Crc*cos(2*thetak); %Radius Correction
delta_ik = Cis*sin(2*thetak)+Cic*cos(2*thetak); %Inclination Correction

uk = thetak+delta_uk; %Corrected Argument of Latitude
rk = A*(1-e*cos(Ek))+delta_rk; %Corrected Radius
ik = i0+delta_ik+IDOT*tk; %Corrected Inclination

%Position in Orbital Plane
xko = rk*cos(uk);
yko = rk*sin(uk);

OmegaK = Omega0+(OmegaDot-OmegaEDot)*tk-OmegaEDot*toe; %Corrected Longitude of Ascending Node

%Earth-Fixed Coordinates of Satellite(Transmission)
xk = xko*cos(OmegaK)-yko*cos(ik)*sin(OmegaK);
yk = xko*sin(OmegaK)+yko*cos(ik)*cos(OmegaK);
zk = yko*sin(ik);

%Earth-Fixed Coordinates of Satellite(Reception)
%ECI is equal to ECEF at transmission time
%Now we have to calculate the SV positon at receiver time in ECEF(changed)
theta = OmegaEDot*(iUser_SoW-iT_SV+delta_tr);
R_ECEF2ECI = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
R_ECI2ECEF = inv(R_ECEF2ECI);
SV_Pos_ECI = [xk;yk;zk];
SV_Pos_ECEF = R_ECI2ECEF*SV_Pos_ECI; % = inv(R_ECEF2ECI)*SV_Pos_ECI
%x = A\b is computed differently than x = inv(A)*b and is recommended for solving systems of linear equations.

vSatellite_xyz = SV_Pos_ECEF.'; %B = A.' -> B = transpose(A)



