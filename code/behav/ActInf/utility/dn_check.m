load([root_path, 'data\behav\action_fix_time.mat'], 'action_time', 'fix_time');
action_time(isnan(action_time)) = 1.5;
load([root_path, 'data\behav\behav_measure.mat'], 'lo_feed', 'lo_cue', 're_feed', 're_cue');

feed_time = [lo_feed, re_feed] - 7.5;
cue_time = [lo_cue, re_cue] - 7.5;

dn = []; wn = [];
for ti = 1:40
    dn(:,ti) = MDP_best(ti).dn;   % [iter(16) * T(2) x 1]
    wn(:,ti) = MDP_best(ti).wn;
end

%%
si = 1; Precision = [];
for ti = 1:40
    if ti == 1
        continue
    else
        for wti = 1:size(dn,1)
            if (wti > 16) && (dn(wti, ti - 1) >= 0)
                Precision(round((cue_time(ti)+(action_time(ti,si)/16)*(wti-1))*10):round((cue_time(ti)+(action_time(ti,si)/16)*(wti-1))*10)+9) = repelem(dn(wti,ti - 1),10);
            end
        end
    end
end

figure
area(1:length(Precision), Precision)