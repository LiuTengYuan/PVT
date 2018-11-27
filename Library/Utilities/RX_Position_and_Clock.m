function [RX_Position_XYZ, RX_ClockError,Matrix] = RX_Position_and_Clock(Result,mC1,Nb_Epoch,Epoch_SV_Number,HMode)

% Result = handles.Result;
% mC1 = handles.mC1;

RX_Position_XYZ = zeros(Nb_Epoch,3);
RX_ClockError = zeros(Nb_Epoch,1);
Matrix = struct();
for epoch=1:Nb_Epoch
    RX_Xbar = 0; RX_Ybar = 0; RX_Zbar = 0; trx_meter = 0;
    delta_X = 1; delta_Y = 1; delta_Z = 1; delta_trx = 1;
    H = zeros(Epoch_SV_Number(epoch),4);
    cc = 0;
    Pseudorange_difference = zeros(Epoch_SV_Number(epoch),1);
    while (abs(delta_X) > 0.0001 || abs(delta_Y) > 0.0001 || abs(delta_Z) > 0.0001) % norm(Corrected_delta(1:3))
        for SV_num = 1:Epoch_SV_Number(epoch)
            SV_X = Result(epoch).SV(SV_num,2); SV_Y = Result(epoch).SV(SV_num,3); SV_Z = Result(epoch).SV(SV_num,4);
            tsv_meter = Result(epoch).SV(SV_num,6);
            Pseudorangebar = sqrt((RX_Xbar-SV_X)^2+(RX_Ybar-SV_Y)^2+(RX_Zbar-SV_Z)^2);
            Pseudorange = mC1(epoch,Result(epoch).SV(SV_num,1));
            %Pseudorange = TrueRange + trx_meter - tsv_meter + Ionosphere_delay + Troposphere_delay + Multipath + Noise
            Pseudorange_difference(SV_num) = Pseudorange-Pseudorangebar+tsv_meter-trx_meter;
            H(SV_num,:) = [(RX_Xbar-SV_X)/Pseudorangebar (RX_Ybar-SV_Y)/Pseudorangebar (RX_Zbar-SV_Z)/Pseudorangebar 1];
        end
        switch HMode
            case 'NLSE'
%                 W = eye(size(H,2));
            case 'NWLSE'
%                 W = SNR WEIGHTED;
%                 Matrix(epoch).W = W;
        end
        Corrected_delta = (H'*H)\H'*Pseudorange_difference;
        Matrix(epoch).H = H;
        
        delta_X = Corrected_delta(1); delta_Y = Corrected_delta(2); delta_Z = Corrected_delta(3);
        delta_trx = Corrected_delta(4);
        RX_Xbar = RX_Xbar + delta_X; RX_Ybar = RX_Ybar + delta_Y; RX_Zbar = RX_Zbar + delta_Z;
        trx_meter = trx_meter + delta_trx;
        cc = cc+1;
    end
    RX_Position_XYZ(epoch,:) = [RX_Xbar RX_Ybar RX_Zbar];
    RX_ClockError(epoch) =  trx_meter;
    clc;
    ProcessingCompleted = epoch/(Nb_Epoch)*100;
    fprintf("\n SV Position Processing completed. ");
    fprintf("\n RX Position Processig... ");
    fprintf(" \n Total Completed - %.2f %% \n",ProcessingCompleted);
end

end