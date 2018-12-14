function [RX_Position_XYZ, RX_ClockError,Matrix] = RX_Position_and_Clock(Result,mC1,mS1,Nb_Epoch,vNb_Sat,HMode,WType,Tiono,Ttropo,CallNumber,Elevation_Azimuth,handles)

% Result = handles.Result;
% mC1 = handles.mC1;

RX_Position_XYZ = zeros(Nb_Epoch,3);
RX_ClockError = zeros(Nb_Epoch,1);
Matrix = struct();
for epoch=1:Nb_Epoch
    RX_Xbar = 0; RX_Ybar = 0; RX_Zbar = 0; trx_meter = 0;
    delta_X = 1; delta_Y = 1; delta_Z = 1; delta_trx = 1;
    H = zeros(vNb_Sat(epoch),4);
    Pseudorange_difference = zeros(vNb_Sat(epoch),1);
    while (abs(delta_X) > 0.0001 || abs(delta_Y) > 0.0001 || abs(delta_Z) > 0.0001) % norm(Corrected_delta(1:3))
        weights = zeros(1,vNb_Sat(epoch));
        for SV_num = 1:vNb_Sat(epoch)
            SV_X = Result(epoch).SV(SV_num,2); SV_Y = Result(epoch).SV(SV_num,3); SV_Z = Result(epoch).SV(SV_num,4);
            tsv_meter = Result(epoch).SV(SV_num,6);
            Pseudorangebar = sqrt((SV_X-RX_Xbar)^2+(SV_Y-RX_Ybar)^2+(SV_Z-RX_Zbar)^2);
            Pseudorange = mC1(epoch,Result(epoch).SV(SV_num,1));
            %Pseudorange = TrueRange + trx_meter - tsv_meter + Ionosphere_delay + Troposphere_delay + Multipath + Noise
            Pseudorange_difference(SV_num) = Pseudorange-Pseudorangebar+tsv_meter-trx_meter-Tiono(epoch,SV_num)-Ttropo(epoch,SV_num);
            H(SV_num,:) = [(RX_Xbar-SV_X)/Pseudorangebar (RX_Ybar-SV_Y)/Pseudorangebar (RX_Zbar-SV_Z)/Pseudorangebar 1];
            if CallNumber > 1
                if (WType == 2) && (strcmp(HMode,'NWLSE'))
                    prn = find(handles.SVTracked == Result(epoch).SV(SV_num,1));
                    % 1 - 15    index 1 - 15              1 - 32
%                     weights(SV_num) = 1/norm(Elevation_Azimuth(prn).llh(epoch,1:2) - RXLLH(epoch,1:2)); % In position SV_num, the 1/norm between SV prn (corresponding to index SV_num) and RX
                    weights(SV_num) = (Elevation_Azimuth(epoch).SV(SV_num,2));
                else
                    weights(SV_num) = 1;
                end
            end
            
        end
        switch HMode
            case 'NLSE'
                W = eye(size(H,1));
            case 'NWLSE'
                %??????????????????????????????????????????????????????
                %Measurement Error Covariance Matrix(Cp)
                %W = inv(Cp)
                SNR_Weighted = mS1(epoch,mS1(epoch,:)~=0) .* weights;
                W = diag(SNR_Weighted);
        end
        Corrected_delta = (H'*W*H)\H'*W*Pseudorange_difference;
        Matrix(epoch).H = H;
        Matrix(epoch).W = W;
        
        delta_X = Corrected_delta(1); delta_Y = Corrected_delta(2); delta_Z = Corrected_delta(3);
        delta_trx = Corrected_delta(4);
        RX_Xbar = RX_Xbar + delta_X; RX_Ybar = RX_Ybar + delta_Y; RX_Zbar = RX_Zbar + delta_Z;
        trx_meter = trx_meter + delta_trx;
    end
    RX_Position_XYZ(epoch,:) = [RX_Xbar RX_Ybar RX_Zbar];
    RX_ClockError(epoch) =  trx_meter;
    clc;
    ProcessingCompleted = round(epoch/(Nb_Epoch)*100);
    if mod(ProcessingCompleted,10)==0
        set(handles.MessageBox,'string',sprintf(" SV Position Completed\n RX Position Processing...\n Total Completed - %.2f %%\n   Call %d out of %d.",ProcessingCompleted,CallNumber,6));
        pause(0.001)
    end
    fprintf("\n RX Position Processig... ");
    fprintf(" \nTotal Completed - %.2f , Call %d out of %d. \n",ProcessingCompleted, CallNumber, 6);
end

clc
fprintf("\n SV Position Processing completed. ");
fprintf("\n RX Position Processing completed. \n");
set(handles.MessageBox,'string',sprintf("\n SV Position Processing completed.\n RX Position Processing completed. \n"));


end