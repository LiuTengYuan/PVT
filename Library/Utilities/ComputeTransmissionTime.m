function [iDelta_t_SV, iT_SV] = ComputeTransmissionTime(vEphemeris, iUser_SoW, f_C1)
%--------------------------------------------------------------------------
% Copyright © ENAC, 2015.
% ENAC : http://www.enac.fr/.
% signav@recherche.enac.fr
%
% This function computes SV PRN signal transmission time in GPS system
% time.
% data (L1 signal).
% Input Variables:
% 1) vEphemeris SV PRN ephemeris data - 1x29 vector
% 2) iUser_SoW user time (Second of Week - SoW)
% 3) f_C1 code range measurement
% Output Variables:
% 1) iDelta_t_SV SV PRN code phase offset [s] (do not include the
% relativistic correction term)
% 2) iT_SV signal transmission time in GPS system time [s]
%
% Reference:
% ICD-GPS-200C, 10 OCT 1993, page 88, page 90
%--------------------------------------------------------------------------

% Initialize Variables
iDelta_t_SV = 0;
iT_SV = 0;
c = 299792458.0; %speed of light(m/s)

%Caculate iDelta_t_SV(without relativistic correction term)
t = iUser_SoW - f_C1/c;
delta_t = t-vEphemeris(4);
af0 = vEphemeris(6); af1 = vEphemeris(7); af2 = vEphemeris(8);
Tgd = vEphemeris(9);
iDelta_t_SV = af0+af1*delta_t+af2*delta_t^2-Tgd;

%Caculate iT_SV
iT_SV = t-iDelta_t_SV;







