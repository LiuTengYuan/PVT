function [RX_Position_LLH, RX_Position_ENU, Matrix, DOP] = RX_Position_LLH_ENU(RX_Position_XYZ, Nb_Epoch, Matrix)

%transforms the rectangular coordinates (x,y,z) in ECEF frame to the geodetic coordinates(Latitude,Longitude,Height) above the reference ellipsoid)
RX_Position_LLH = zeros(Nb_Epoch,3);
for epoch=1:Nb_Epoch
    vXYZ = RX_Position_XYZ(epoch,:);
    vLLH = xyz_2_lla_PVT(vXYZ);
    vLLH(1) = vLLH(1)*180/pi; vLLH(2) = vLLH(2)*180/pi;
    RX_Position_LLH(epoch,:) = vLLH;
end

% % ---
% RX_Position_LLH_test = zeros(Nb_Epoch,3);
% for epoch=1:Nb_Epoch
%     xyz = RX_Position_XYZ(epoch,:);
%     [lat, lon, alt] = wgsxyz2lla(xyz);
%     RX_Position_LLH_test(epoch,:) = [lat, lon, alt];
% end
% % ---

mean_XYZ = mean(RX_Position_XYZ);
stdev_XYZ = std(RX_Position_XYZ);

RX_Position_ENU = zeros(Nb_Epoch,3);
vXYZ0 = mean_XYZ;
for epoch=1:Nb_Epoch
    vXYZw = RX_Position_XYZ(epoch,:);
    vXYZw_dir = vXYZw-vXYZ0;
    [vXYZl, mTRANSF] = delta_wgs84_2_local(vXYZw_dir', vXYZ0');
    RX_Position_ENU(epoch,:) = [vXYZl(2) vXYZl(1) vXYZl(3)];
    mTRANSF_extended = [mTRANSF zeros(size(mTRANSF,1),1); zeros(1,size(mTRANSF,2)) 1];
    G = Matrix(epoch).H*mTRANSF_extended.'; % Mx4x(4x4).'=Mx4->By making.'every H entry multiplied by ENU coefficients.
    Matrix(epoch).COV = inv(G.' * Matrix(epoch).W * G);
    tmp = sqrt(diag(Matrix(epoch).COV));
    DOP.EDOP(epoch) = tmp(1);
    DOP.NDOP(epoch) = tmp(2);
    DOP.VDOP(epoch) = tmp(3);
    DOP.TDOP(epoch) = tmp(4);
end

% % ---
% RX_Position_ENU_test = zeros(Nb_Epoch,3);
% [reflat, reflon, refalt] = wgsxyz2lla(mean_XYZ);
% for epoch=1:Nb_Epoch
%     xyz = RX_Position_XYZ(epoch,:);
%     enu = wgsxyz2enu(xyz', reflat, reflon, refalt);
%     RX_Position_ENU_test(epoch,:) = enu';
% end
% % ---

end