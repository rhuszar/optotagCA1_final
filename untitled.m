
% Here, we picked a session with decent behavior, extracted the theta
% cycles during maze running (>2cm/s running speed), and extract a bunch of
% theta cycles


powerThresh = 1;
samplingRate = 1250;
passband = [6 14];

basepath = 'C:\Users\Roman\Google Drive\buzz_workspace\optotag\DATA\e14\e14_2m3\e14_2m3_201118';
% Load behavior data + LFP

basename = bz_BasenameFromBasepath(basepath); 

load( fullfile(basepath, [basename '.ripples.events.mat']) )
load( fullfile(basepath, [basename '.Behavior.mat']) )

data = ripples.detectorinfo.detectionparms.lfp;
timestamps  = [ 0:(1/1250):length( data ) * (1/1250) ]'; timestamps(end) = [];

    lfp = struct();
    lfp.data = data;
    lfp.timestamps = timestamps;
    lfp.channels = 1;

lfp_th = bz_Filter( lfp, 'passband', passband );
%%
tvec_beh = behavior.trials.position_trcat.timestamps;
v = behavior.trials.position_trcat.v;
trial_ints = behavior.trials.trial_ints;
cycle_duration =  0.1;
intervals = [];
for kp = 1:length(trial_ints)
    ind = InIntervals(tvec_beh, trial_ints(kp, :));
    % plot(tvec_beh(ind), v(ind) )
    idx.states = double( v(ind) > 2 );
    idx.statenames = {'FAST'};
    idx.timestamps = tvec_beh(ind);
    int = bz_IDXtoINT( idx );
    fast_dur = diff( int.FASTstate' );
    % Only consider intervals with at least 5 cycles
    for jp = find( fast_dur / cycle_duration > 5 )
        
        
        % intervals = [intervals ; int.FASTstate(1,:)];
        ip = find( InIntervals(lfp_th.timestamps, int.FASTstate(jp,:)) );
        % Identify the theta peaks
        [~, jj] = findpeaks(-abs( lfp_th.phase(ip) ));
        for c = 3:5:length(jj)
            intervals = [intervals ; [ip( jj(c) )-250 ip( jj(c) )+250 ]];
        end
%         if length( jj ) > 20
%             intervals = [ intervals ; [ip( jj(10) )  - 625 ip( jj(10) )  + 625] ];
%         end
%         subplot(2,1,1)
%         plot( lfp_th.timestamps(ip), lfp_th.data(ip) )
%        subplot(2,1,2)
%         plot( lfp_th.timestamps(ip), lfp_th.phase(ip) )
    end
%     shg
%     uiwait
end
%%
cycles_data = []; phase_data = [];
for kp = 1:length(intervals)
    cycles_data = [cycles_data ;  double( data(intervals(kp,1):intervals(kp,2))') ];
    phase_data = [phase_data ; lfp_th.phase(intervals(kp,1):intervals(kp,2))'];
end

%%
subplot(2,1,1)
plot(nanmean(cycles_data))
subplot(2,1,2)
imagesc(phase_data)

%%

inds = 110:386;
lfp_mu = nanmean(cycles_data(:, inds));

[~, ip] = min( abs(  phase_data(:,inds(1)) ) + abs(  phase_data(:,inds(end)) ) );
phase = phase_data(ip,inds); phase(1) = abs(phase(1));
ph = [ rad2deg( mod( phase(1:142), 2*pi ) ) 360+rad2deg( mod( phase(143:end), 2*pi ) ) ];

plot(ph, lfp_mu)
%%


%%
inds = 326:627;
lfp_mu = nanmean(cycles_data([ 1:4 6:10], inds));
phase_mu = nanmean( phase_data([ 1:4 6:10],inds) );
% ph = [ rad2deg( mod( phase_mu(1:143), 2*pi ) ) 360+rad2deg( mod( phase_mu(144:end), 2*pi ) ) ];
phase_mu = phase_data(1,inds);
ph = [ rad2deg( mod( phase_mu(1:141), 2*pi ) ) 360+rad2deg( mod( phase_mu(142:end-1), 2*pi ) ) ];

plot( [ ph ph(end)+2.2484] ,  lfp_mu)
xlim([0 720])

save

%%



% See if there is behavior
try
    load( fullfile(basepath, [basename '.MergePoints.events.mat']) )
    beh_ts = MergePoints.timestamps(2,:);
catch
    beh_ts = [];
end

% Use the ripple LFP
data = ripples.detectorinfo.detectionparms.lfp;
lfp = struct();
lfp.data = data;

% filter the LFP in the desired band
filtered = bz_Filter(lfp, 'passband', passband);

% Define recording interval, get rid of pulses
intervals = [0 filtered.timestamps(end)];
intervals = SubtractIntervals(intervals, pulses.intsPeriods);

power = filtered.amp;
hilb = filtered.hilb;
lfpphase = mod( filtered.phase,2*pi);

thresh = mean(power) + std(power)*powerThresh;
minWidth = (samplingRate./passband(2)) * 2; % set the minimum width to two cycles

below = find(power < thresh);
if max(diff(diff(below))) == 0
    below_thresh = [below(1) below(end)];
elseif length(below)>0;
    ends=find(diff(below)~=1);
    ends(end+1)=length(below);
    ends=sort(ends);
    lengths=diff(ends);
    stops=below(ends)./samplingRate;
    starts=lengths./samplingRate;
    starts = [1; starts];
    below_thresh(:,2)=stops;
    below_thresh(:,1)=stops-starts;
else
    below_thresh=[];
end
% now merge interval sets from input and power threshold
intervals = SubtractIntervals(intervals,below_thresh);
intervals = intervals(diff(intervals')>minWidth./samplingRate,:); 

%%
%Treat each interval like an event - apply the ripple routine
theta_cycle_ints = [];
[~, indices,~] = InIntervals( filtered.timestamps, intervals );
for kp = 1481:length(intervals)
    
    ts = filtered.timestamps(indices == kp);
    phase = filtered.phase(indices == kp);
    
    [~, peak_indices] = findpeaks( -abs( 0 - phase') );
    
    
    % Add missing peak at the start
    if phase(1)  > 0 && phase(1) < 0.4
        peak_ts = [ ts(1) ts(peak_indices) ];
        peak_indices = [1 peak_indices];
    elseif (peak_indices(1) ~= 1 && peak_indices(1) ~= 2) && ...
            phase(1) < 0 && phase(2) > 0
        tmp = [ phase(1) phase(2) ];
        [~, ii] = min( abs(tmp) );
        peak_ts = [ ts(ii) ts(peak_indices) ];
        peak_indices = [ii peak_indices];
    else
        peak_ts = ts(peak_indices);
    end

     % Add missing peak at the end
    if phase(end) < 0 && phase(end) > -0.4
        peak_ts = [ ts(peak_indices) ts(end) ];
        peak_indices = [peak_indices length(ts) ];
    elseif (peak_indices(end) ~= length(ts) && peak_indices(end) ~= length(ts)-1) && ...
            phase(end-1) < 0 && phase(end) > 0
        tmp = [ phase(end) phase(end-1) ];
        [~, ii] = min( abs(tmp) );
        peak_ts = [ ts(peak_indices) ts(end-(ii-1)) ];
        peak_indices = [peak_indices length(ts)-(ii-1) ];
    else
        peak_ts = ts(peak_indices);
    end
    
    if length(peak_ts) >= 4
        theta_cycle_ints = [theta_cycle_ints ; peak_ts(1)];
    end
    
    
%     plot(ts, phase)
%     hold on
%     plot(peak_ts, phase(peak_indices), 'kx')
%     plot(ts, filtered.data(indices)*0.01)
%     uiwait
%     
end