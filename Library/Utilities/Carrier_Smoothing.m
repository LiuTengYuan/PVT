function [smoothed] = Carrier_Smoothing(mC1,mL1,N)

smoothed = zeros(size(mC1));
for sv=1:length(mC1(1,:))
    if sum(mC1(:,sv))
        tracked_epoch = find(mC1(:,sv)~=0);
        first_tracked_epoch = tracked_epoch(1);
        last_tracked_epoch = tracked_epoch(end);
        smoothed(first_tracked_epoch,sv) = mC1(first_tracked_epoch,sv);
        for epoch=first_tracked_epoch:last_tracked_epoch
            
        end
    end
end

end