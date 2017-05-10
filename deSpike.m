function [dout, spiket] = deSpike(din, Fs, varargin)
%------------------------------------------------------------------------
% [dout, spiket] = deSpike(in, varargin)
%------------------------------------------------------------------------
% IntraHEKAToolbox
%------------------------------------------------------------------------
% 
% removes spikes from intracellular trace 
% 
% *requires AudioToolbox and UtilitiesToolbox from TytoLogy project*
%
%------------------------------------------------------------------------
% Input Arguments:
%	din	raw spike trace
%	Fs		sample rate for din trace
% 
% 	Optional:
%		'SpikeThreshold', <spike threshold, volts>, default = -0.015
% 		'DV1', <dV/dt onset threshold, volts/sec>, default = 2 
% 		'DV2', <dV/dt offset threshold, volts/sec>, default = -1.5
%		'SnipTime', <snip size for spline fit, msec>, default = 20
%		'RefractoryTime', <spike refractory time, ms>, default = 2 ms
% 
%		Filtering options (see lpfilter command)
%	 		'MeanTime', <time to compute mean, ms>, default = 1 ms
%			'WinTime',	<time to add to start and end of snippet window, ms>, 
%								default = 5 ms
%
%
% Output Arguments:
%	dout		despiked trace (at sampling rate of Fs)
%	spiket	times of detected spikes (seconds)
%------------------------------------------------------------------------
% See also: filtfilt, butter, lpfilter, spikeschmitt3
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 1 May, 2012 (SJS)
%
% Revisions:
%	17 Oct 2016 (SJS): some cleanup
%	31 Oct 2016 (SJS): adding options
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%--------------------------------------------------------------
% settings
%--------------------------------------------------------------
% # of points to compute mean (based on time, ms)
Meantime = 1;
Meanpts = ms2samples(Meantime, Fs);
% window size (ms) and points to add to start and end
Wintime = 5;
Winpts = ms2samples(Wintime, Fs);
% threshold for spikes, V
spikethresh = -0.015;
% dV/dt onset, offset threshold, V/s
dvthresh1 = 2.0;
dvthresh2 = -1.50;
% size of snip to use
sniptime = 20;
% length of snip in samples
snipwin = ms2samples(sniptime, Fs);
% spike refractory period, ms
refractorytime = 2;

%--------------------------------------------------------------
% process args
%--------------------------------------------------------------
if ~isempty(varargin)
	optargs = length(varargin);
	% parse varargin args
	n = 1;
	while n < optargs
		switch upper(varargin{n})
			case 'SPIKETHRESHOLD'
				spikethresh = varargin{n+1};
				n = n + 2;
			case 'DV1'
				dvthresh1 = varargin{n+1};
				n = n + 2;
			case 'DV2'
				dvthresh2 = varargin{n+1};
				n = n + 2;
			case 'MEANTIME'
				Meantime = varargin{n+1};
				Meanpts = ms2samples(Meanpts, Fs);
				n = n + 2;
			case 'WINTIME'
				Wintime = varargin{n+1};
				Winpts = ms2samples(Winpts, Fs);
				n = n + 2;
			case 'SNIPTIME'
				sniptime = varargin{n+1};
				snipwin = ms2samples(sniptime, Fs);
				n = n + 2;
			case 'REFRACTORYTIME'
				refractorytime = varargin{n+1};
				n = n + 2;
			otherwise
				error('%s: bad arg %s', mfilename, varargin{n});
		end
	end
end


%--------------------------------------------------------------
% filter data
%--------------------------------------------------------------
% first, define a low-pass filter  
% sample interval
dt = 1 ./ Fs;
% low-pass cutoff frequency
fc = 250;
% order of filter
forder = 1;
% then filter data
fdata = lpfilter(din, Fs, fc, 'MeanPad', ...
							[Winpts Meanpts], 'FilterOrder', forder);

%--------------------------------------------------------------
% find spike times (bins)
%--------------------------------------------------------------
spiket = spikeschmitt2(force_row(din) - spikethresh, ...
					0, refractorytime, Fs);
nspikes = length(spiket);

%--------------------------------------------------------------
% remove detected spikes
%--------------------------------------------------------------
% default: assign input vector to output
dout = din;

% if no spikes detected, return
if isempty(spiket) || (nspikes == 0)
	return
end

% loop through spikes
snipstart = zeros(nspikes, 1);
snipend = zeros(nspikes, 1);
snippts = cell(nspikes, 1);
snip = cell(nspikes, 1);
maxval = zeros(nspikes, 1);
maxpt = zeros(nspikes, 1);
dsnip = cell(nspikes, 1);
dmaxval = zeros(nspikes, 1);
dmaxpt = zeros(nspikes, 1);
dminval = zeros(nspikes, 1);
dminpt = zeros(nspikes, 1);
pstart = zeros(nspikes, 1);
pend = zeros(nspikes, 1);
spikestart = zeros(nspikes, 1);
spikeend = zeros(nspikes, 1);
m =  zeros(nspikes, 1);
midpt_x = zeros(nspikes, 1);
midpt_y =  zeros(nspikes, 1);
pp = cell(nspikes, 1);

for s = 1:length(spiket)

	% start and end of snippet re: full trace
	snipstart(s) = spiket(s) - snipwin;
	snipend(s) = spiket(s) + snipwin;

	% detect if snipstart or snipend are out of range
	if snipstart(s) < 1
		snipstart(s) = 1;
	end
	if snipend(s) > length(fdata)
		snipend(s) = length(fdata);
	end

	% vector of snip points (indices and values) re: full trace
	snippts{s} = snipstart(s):snipend(s);
	snip{s} = fdata(snippts{s});

	% find max val of snip
	[maxval(s), maxpt(s)] = max(snip{s});

	% get gradient (differentiate) the snip
	dsnip{s} = gradient(snip{s}, dt);
	% find max, min values and points of dsnip, re: snip
	[dmaxval(s), dmaxpt(s)] = max(dsnip{s});
	[dminval(s), dminpt(s)] = min(dsnip{s});

	% start of spike is first point of spike before peak of dspike that
	% crosses dvthresh1, re: snippet
	tmp = find(dsnip{s}(1:dmaxpt(s)) >= dvthresh1, 1, 'first');
	if isempty(tmp)
		pstart(s) = 1;
	else
		pstart(s) = tmp;
	end
	clear tmp
	% end of spike is first point after minimum dspike that crosses dvthresh2
	% dminpt(s) is added to the found point in order to align pend with the 
	% sample indices of the snippet snip{s}, re: snippet
	tmp = dminpt(s) + find(dsnip{s}(dminpt(s):end) >= dvthresh2, 1, 'first');
	if isempty(tmp)
		pend(s) = dminpt(s);
	else
		pend(s) = tmp;
	end
	if pend(s) > length(dsnip{s})
		pend(s) = length(dsnip{s});
	end
	clear tmp

	% spikestart and spikeend are the spike start/end points relative to the 
	% full fdata trace, d
	spikestart(s) = snipstart(s) + pstart(s);
	spikeend(s) = snipstart(s) + pend(s);
	if spikeend(s) > length(fdata)
		spikeend(s) = length(fdata);
	end
	if spikestart(s) < 1
		spikestart(s) = 1;
	end

	try
		% plot full trace with points/snippet identified
		figure(99)
		subplot(211);
		plot(fdata);
		hold on
			plot(snippts{s}, snip{s}, 'y.');
			plot(spiket, fdata(spiket), 'k.');
			plot(spikestart(s):spikeend(s), fdata(spikestart(s):spikeend(s)), 'g.');
			plot(spikeend(s), fdata(spikeend(s)), 'r.');
		hold off
	catch errObj
		fprintf('%s: trapped error\n', mfilename)
		fprintf('%s\n', errObj.message);
		keyboard
	end
	%----------------------------------------------------
	% fit spline to snipped out fdata
	%----------------------------------------------------
	% get the relevant points in front of and after the spike, including a
	% midpoint that has a y-value on a line between the start and end of the 
	% spike

	try
		y2 = snip{s}(pend(s)); %#ok<*NASGU>
		y1 = snip{s}(pstart(s));
		x2 = pend(s);
		x1 = pstart(s);
		m(s) = (snip{s}(pend(s)) - snip{s}(pstart(s))) / (pend(s) - pstart(s));
		midpt_x(s) = dminpt(s);
		midpt_y(s) = (m(s) * (midpt_x(s) - pstart(s))) + snip{s}(pstart(s));
% 			midpt_y(s) = snip{s}(pend(s))
% 			ppts = snipstart(s) + [(1:pstart(s)) midpt_x(s) (pend(s):length(snip{s}))];
% 			ppts = [(1:pstart(s)) midpt_x(s) (pend(s):length(snip{s}))];
		splinewin = 0;
		s1 = (pstart(s) - splinewin):pstart(s);
		s2 = pend(s):(pend(s) + splinewin);
		ppts = [s1 midpt_x(s) s2];
% 			xsnip = snippts{s}(ppts);
		xsnip = ppts;
		ysnip = [snip{s}(s1) midpt_y(s) snip{s}(s2)];

		% use this to get spline fit
		pp{s} = spline(xsnip, ysnip);

		% get fit points
		xx = pstart(s):pend(s);
		yy = ppval(pp{s}, xx);

		dout(spikestart(s):spikeend(s)) = yy;
	catch errObj
		fprintf('%s\n', errObj.message);
		fprintf('error in spike removal, aborting this spike...\n')
	end

end		% END SPIKE LOOP (s)

subplot(212)
plot(din, 'g');
hold on
	plot(dout);
hold off
drawnow




