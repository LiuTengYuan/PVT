function [RX_Velocity_XYZ, RX_ClockErrorRate] = RX_Velocity_and_Clock_Rate(Result,mD1,Nb_Epoch,vNb_Sat,Matrix)

global c;
frequency_L1 = 1.57542e9; %1575.42 MHz
wavelength_L1 = c/frequency_L1;

RX_Velocity_XYZ = zeros(Nb_Epoch,3);
RX_ClockErrorRate = zeros(Nb_Epoch,1);

for epoch=1:Nb_Epoch
    RX_VXbar = 0; RX_VYbar = 0; RX_VZbar = 0; trx_rate = 0;
    delta_VX = 1; delta_VY = 1; delta_VZ = 1; delta_trx_rate = 1;
    H = zeros(vNb_Sat(epoch),4);
    Pseudorange_rate_difference = zeros(vNb_Sat(epoch),1);
    while (abs(delta_VX) > 0.0001 || abs(delta_VY) > 0.0001 || abs(delta_VZ) > 0.0001) % norm(Corrected_delta(1:3))
        for SV_num = 1:vNb_Sat(epoch)
            SV_VX = Result(epoch).SV(SV_num,7); SV_VY = Result(epoch).SV(SV_num,8); SV_VZ = Result(epoch).SV(SV_num,9);
            tsv_rate = Result(epoch).SV(SV_num,10);
            Pseudorange_rate_bar = -(Matrix(epoch).H(SV_num,1)*(SV_VX-RX_VXbar)+Matrix(epoch).H(SV_num,2)*(SV_VY-RX_VYbar)+Matrix(epoch).H(SV_num,3)*(SV_VZ-RX_VZbar));
            Pseudorange_rate = mD1(epoch,Result(epoch).SV(SV_num,1))*wavelength_L1;
            Pseudorange_rate_difference(SV_num) = Pseudorange_rate-Pseudorange_rate_bar+tsv_rate-trx_rate;
        end
        
        H = Matrix(epoch).H;
        W = Matrix(epoch).W;
        Corrected_delta = (H'*W*H)\H'*W*Pseudorange_rate_difference;
        
        delta_VX = Corrected_delta(1); delta_VY = Corrected_delta(2); delta_VZ = Corrected_delta(3);
        delta_trx_rate = Corrected_delta(4);
        RX_VXbar = RX_VXbar + delta_VX; RX_VYbar = RX_VYbar + delta_VY; RX_VZbar = RX_VZbar + delta_VZ;
        trx_rate = trx_rate + delta_trx_rate;
    end
    RX_Velocity_XYZ(epoch,:) = [RX_VXbar RX_VYbar RX_VZbar];
    RX_ClockErrorRate(epoch) =  trx_rate;
    clc;
    ProcessingCompleted = epoch/(Nb_Epoch)*100;
    fprintf("\n SV Position Processing completed. ");
    fprintf("\n RX Position Processing completed. ");
    fprintf("\n RX Velocity Processig... ");
    fprintf(" \n Total Completed - %.2f %% \n",ProcessingCompleted);
end

clc
fprintf("\n SV Position Processing completed. ");
fprintf("\n RX Position Processing completed. ");
fprintf("\n RX Velocity Processing completed. \n");


end