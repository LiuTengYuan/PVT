function [cycle_slip] = Cycle_Slip_Detection(mC1,mL1,N,CMC)

cycle_slip = zeros(size(CMC));
md = zeros(size(CMC));
md_2 = zeros(size(CMC));
Sd = zeros(size(CMC));
nt = 6; %confidence range of threshold
So = 1; %predefined initial value for sigma
for sv=1:length(CMC(1,:))
    if sum(CMC(:,sv))
        tracked_epoch = find(CMC(:,sv)~=0);
        first_tracked_epoch = tracked_epoch(1);
        last_tracked_epoch = tracked_epoch(end);
        count = 1;
        for epoch=first_tracked_epoch:last_tracked_epoch
            if mC1(epoch,sv)==0 || mL1(epoch,sv)==0
                cycle_slip(epoch,sv)=1;
                count = 1;
            else
                if count==1
                    cycle_slip(epoch,sv) = 1;
                end
                if count>1 && count<=N
                    for n = 1:(count-1)
                        md(epoch,sv) = md(epoch,sv)+CMC(epoch-n,sv);
                        md_2(epoch,sv) = md_2(epoch,sv)+CMC(epoch-n,sv)^2;
                    end
                    md(epoch,sv) = md(epoch,sv)/(count-1);
                    md_2(epoch,sv) = md_2(epoch,sv)/(count-1);
                    Sd(epoch,sv) = sqrt(md_2(epoch,sv)-md(epoch,sv)^2);
                    Sd(epoch,sv) = sqrt((count-1)/count*Sd(epoch,sv)^2+1/count*So^2);
                end
                if count>N
                    for n = 1:N
                        md(epoch,sv) = md(epoch,sv)+CMC(epoch-n,sv);
                        md_2(epoch,sv) = md_2(epoch,sv)+CMC(epoch-n,sv)^2;
                    end
                    md(epoch,sv) = md(epoch,sv)/N;
                    md_2(epoch,sv) = md_2(epoch,sv)/N;
                    Sd(epoch,sv) = sqrt(md_2(epoch,sv)-md(epoch,sv)^2);
                end
                threshold = nt*Sd(epoch,sv);
                if count>1 && norm(CMC(epoch,sv)-md(epoch,sv))>threshold
                    cycle_slip(epoch,sv) = 1;
                    count = 0;
                end
                count = count+1;
            end
        end
    end
end


end
