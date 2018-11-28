function [Tiono] = Ionospheric_Correction(Iono_a, Iono_b, Elevation_Azimuth, RX_Position_LLH, mEpoch, Total_Nb_Sat)

%--------------------------------------------------------------------------
% Satllite Transmitted Terms
%     (1)alpha_n: the coefficeient of a cubic equation representing the amplitude
%                of the verticcal delay (4 coefficients-8 bits each)
%     (2)beta_n: the coefficeient of a cubic equation representing the period of the
%                model (4 coefficients-8 bits each)
% Receiver generated Terms
%     (1)E: elevation angle between the user and satellite (semi-circles)
%     (2)A: azimuth angle between the user and satellite, measured clockwise positive
%           from the true north (semi-circles)
%     (3)phi_u: user geodetic latitude (semi-circles) WGS 84
%     (4)lambda_u: user geodetic longtitude (semi-circles) WGS 84
%     (5)GPS_time: receiver computed system time
% Computed Terms
%     (1)X: phase (radians)
%     (2)F: obliquity factor (dimensionless)
%     (3)t: local time (sec)
%     (4)phi_m: geomagnetic latitude of the earth projection of the ionospheric
%               intersection point (mean ionospheric height assumed 350km) (semi circles)
%     (5)lambda_i: geodetic longitude of the earth projection of the ionospheric
%                  intersection point (semi-circles)
%     (6)phi_i: geodetic latitude of the earth projection of the ionospheric intersection
%               point (semi-circles)
%     (7)psi: earth's central angle between the user position and the earth projection of
%             ionospheric intersection point (semi-circles)
%--------------------------------------------------------------------------

global c;
GPS_time = mEpoch(:,2);
alpha_n = Iono_a;   beta_n = Iono_b;
phi_u = deg2rad(RX_Position_LLH(:,1));   lambda_u = deg2rad(RX_Position_LLH(:,2));
Nb_Epoch = length(mEpoch);
E = zeros(Nb_Epoch,Total_Nb_Sat);   A = zeros(Nb_Epoch,Total_Nb_Sat);
SVnum = zeros(Nb_Epoch,1);
for epoch=1:Nb_Epoch
    SVnum(epoch) = length(Elevation_Azimuth(epoch).SV(:,1));
    E(epoch,1:SVnum(epoch)) = deg2rad(Elevation_Azimuth(epoch).SV(:,2)');
    A(epoch,1:SVnum(epoch)) = deg2rad(Elevation_Azimuth(epoch).SV(:,3)');
end

psi = zeros(Nb_Epoch,Total_Nb_Sat);
phi_i = zeros(Nb_Epoch,Total_Nb_Sat);
lambda_i = zeros(Nb_Epoch,Total_Nb_Sat);
phi_m = zeros(Nb_Epoch,Total_Nb_Sat);
t = zeros(Nb_Epoch,Total_Nb_Sat);
AMP = zeros(Nb_Epoch,Total_Nb_Sat);
PER = zeros(Nb_Epoch,Total_Nb_Sat);
X = zeros(Nb_Epoch,Total_Nb_Sat);
F = zeros(Nb_Epoch,Total_Nb_Sat);
Tiono = zeros(Nb_Epoch,Total_Nb_Sat);

for epoch=1:Nb_Epoch
    for num_SV = 1:SVnum(epoch)
        %Calculate the earth-centred angle (elevation  in semicircles)
        psi(epoch,num_SV) = 0.0137/(E(epoch,num_SV)+0.11)-0.022;
        %Compute the latitude of the Ionospheric Pierce Point (IPP)
        phi_i(epoch,num_SV) = phi_u(epoch)+psi(epoch,num_SV)*cos(A(epoch,num_SV));
        if phi_i(epoch,num_SV)>0.416
            phi_i(epoch,num_SV) = 0.416;
        end
        if phi_i(epoch,num_SV)<-0.416
            phi_i(epoch,num_SV) = -0.416;
        end
        %Compute the longitude of the IPP
        lambda_i(epoch,num_SV) = lambda_u(epoch)+psi(epoch,num_SV)*sin(A(epoch,num_SV))/cos(phi_i(epoch,num_SV));
        %Find the geomagnetic latitude of the IPP
        phi_m(epoch,num_SV) = phi_i(epoch,num_SV)+0.064*cos(lambda_i(epoch,num_SV)-1.617);
        %Find the local time at the IPP
        t(epoch,num_SV) = mod(4.32e4*lambda_i(epoch,num_SV)+GPS_time(epoch),86400);
        %Compute the amplitude of ionospheric delay
        for n = 1:4
            AMP(epoch,num_SV) = AMP(epoch,num_SV)+alpha_n(n)*phi_m(epoch,num_SV)^(n-1);
        end
        if AMP(epoch,num_SV)<0
            AMP(epoch,num_SV) = 0;
        end
        %Compute the period of ionospheric delay
        for n = 1:4
            PER(epoch,num_SV) = PER(epoch,num_SV)+beta_n(n)*phi_m(epoch,num_SV)^(n-1);
        end
        if PER(epoch,num_SV)<72000
            PER(epoch,num_SV) = 72000;
        end
        %Compute the phase of ionospheric delay
        X(epoch,num_SV) = 2*pi*(t(epoch,num_SV)-50400)/PER(epoch,num_SV);
        %Compute the slant factor (elevation  in semicircles)
        F(epoch,num_SV) = 1.0+16.0*(0.53-E(epoch,num_SV))^3;
        %Compute the slant factor (elevation  in semicircles)
        if abs(X(epoch,num_SV))<1.57
            Tiono(epoch,num_SV) = F(epoch,num_SV)*(5e-9+AMP(epoch,num_SV)*(1-X(epoch,num_SV)^2/2+X(epoch,num_SV)^4/24));
        else
            Tiono(epoch,num_SV) = F(epoch,num_SV)*5e-9;
        end
        Tiono(epoch,num_SV) = Tiono(epoch,num_SV)*c;
    end
end
end


