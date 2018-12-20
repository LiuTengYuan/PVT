function [smoothed,cycle_slip] = Carrier_Smoothing(mC1,mL1,N,CMC)

global lambda;
cycle_slip = Cycle_Slip_Detection(mC1,mL1,N,CMC);
smoothed = zeros(size(mC1));
for sv=1:length(mC1(1,:))
    if sum(mC1(:,sv))
        tracked_epoch = find(mC1(:,sv)~=0);
        first_tracked_epoch = tracked_epoch(1);
        last_tracked_epoch = tracked_epoch(end);
        count = 1;
        for epoch=first_tracked_epoch:last_tracked_epoch
            if cycle_slip(epoch,sv) %cycle slip detection
                count = 1;
            end   
%             % simple way to implenment cycle clip detection
%             if mC1(epoch,sv)==0 || mL1(epoch,sv)==0
%                 count = 1;
%             end
%             if epoch>1 && (mC1(epoch-1,sv)==0 || mL1(epoch-1,sv)==0)
%                 count = 1;
%             end
            if count==1
                smoothed(epoch,sv) = mC1(epoch,sv);
            elseif count>1 && count<N
                smoothed(epoch,sv) = 1/count*mC1(epoch,sv)+(count-1)/count*...
                    (smoothed(epoch-1,sv)+lambda*(mL1(epoch,sv)-mL1(epoch-1,sv)));
            else %count>=N
                smoothed(epoch,sv) = 1/N*mC1(epoch,sv)+(N-1)/N*...
                    (smoothed(epoch-1,sv)+lambda*(mL1(epoch,sv)-mL1(epoch-1,sv)));
            end
            count = count+1;
        end
    end
end

end