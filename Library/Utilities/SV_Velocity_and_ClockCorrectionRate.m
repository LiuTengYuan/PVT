function [Result, Result_Info] = SV_Velocity_and_ClockCorrectionRate(Result,Result_Info,mTracked,SVTracked)

for sv=SVTracked
    for epoch=1:length(mTracked)
        epoch_tracked = find(mTracked(:,sv)~=0)';
        if mTracked(epoch,sv)
            now = find(epoch_tracked==epoch);
            target_sv = find(Result(epoch).SV(:,1)==sv);
            if now~=1 && now~=length(epoch_tracked)
                before_epoch = epoch_tracked(now-1);
                after_epoch = epoch_tracked(now+1);
                
                target_sv_before = find(Result(before_epoch).SV(:,1)==sv);
                target_sv_after = find(Result(after_epoch).SV(:,1)==sv);
                
                distance = Result(after_epoch).SV(target_sv_after,2:4)-Result(before_epoch).SV(target_sv_before,2:4);
                diff_sv_clock = Result(after_epoch).SV(target_sv_after,6)-Result(before_epoch).SV(target_sv_before,6);
                time = after_epoch-before_epoch;
            elseif now==1
                after_epoch = epoch_tracked(now+1);
                
                target_sv_after = find(Result(after_epoch).SV(:,1)==sv);
                
                distance = Result(after_epoch).SV(target_sv_after,2:4)-Result(epoch).SV(target_sv,2:4);
                diff_sv_clock = Result(after_epoch).SV(target_sv_after,6)-Result(epoch).SV(target_sv,6);
                time = after_epoch-epoch;
            else %now==length(epoch_tracked)
                before_epoch = epoch_tracked(now-1);
                
                target_sv_before = find(Result(before_epoch).SV(:,1)==sv);
                
                distance = Result(epoch).SV(target_sv,2:4)-Result(before_epoch).SV(target_sv_before,2:4);
                diff_sv_clock = Result(epoch).SV(target_sv,6)-Result(before_epoch).SV(target_sv_before,6);
                time = epoch-before_epoch;
            end
            Result(epoch).SV(target_sv,5:7) = distance/time;
            Result(epoch).SV(target_sv,10) = diff_sv_clock/time;
        end
    end
end

new_Titles = {'SV_XVelocity','SV_YVelocity','SV_ZVelocity','SV_ClockError_Rate'};
new_Result_Info(length(mTracked)) = struct();
for epoch=1:length(mTracked)
    new_Result_Info(epoch).INFO = new_Titles;
    new_Result_Info(epoch).INFO = [new_Result_Info(epoch).INFO; num2cell([Result(epoch).SV(:,7:10)])];
end
for epoch=1:length(mTracked)
    Result_Info(epoch).INFO = horzcat(Result_Info(epoch).INFO,new_Result_Info(epoch).INFO);
end

end