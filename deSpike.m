function [dout, spiket] = deSpike(din, Fs, varargin)
%-----------------------------------------------------------------------------
% [dout, spiket] = deSpike(in, varargin)
%-----------------------------------------------------------------------------
% IntraHEKAToolbox
%-----------------------------------------------------------------------------
% 
% removes spikes from intracellular trace 
% 
%-----------------------------------------------------------------------------
% Input Arguments:
%		
% Output Arguments:
%
%-----------------------------------------------------------------------------
% See also: filtfilt, butter
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 1 May, 2012 (SJS)
%
% Revisions:
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------


%--------------------------------------------------------------
% define a low-pass filter  
%--------------------------------------------------------------
% sample interval
dt = 1 ./ Fs;
% low-pass cutoff frequency
fc = 250;
% order of filter
forder = 1;
% # of points to compute mean
Meanpts = 10;
% window size (ms) and points to add to start and end
Winsize = 5;
Winpts = ms2samples(Winsize, Fs);

% threshold for spikes, V
spikethresh = -0.015;
% dV/dt onset, offset threshold, V/s
dvthresh1 = 2.0;
dvthresh2 = -1.50;
% size of snip to use
sniptime = 20;
% length of snip in samples
snipwin = ms2samples(sniptime, Fs);

% filter data
fdata = lpfilter(din, Fs, fc, 'MeanPad', [Winpts Meanpts], 'FilterOrder', forder);

% find spike times (bins)
spiket = spikeschmitt2(force_row(din) - spikethresh, 0, 2, Fs);

dout = din;

if isempty(spiket)
	return
else
	% loop through spikes
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
			pend(s) = dminpt(s)
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
		subplot(211);
		plot(fdata);
		hold on
			plot(snippts{s}, snip{s}, 'y.');
			plot(spiket, fdata(spiket), 'k.');
			plot(spikestart(s):spikeend(s), fdata(spikestart(s):spikeend(s)), 'g.');
			plot(spikeend(s), fdata(spikeend(s)), 'r.');
		hold off
		catch
		keyboard
		end
		%----------------------------------------------------
		% fit spline to snipped out fdata
		%----------------------------------------------------
		% get the relevant points in front of and after the spike, including a
		% midpoint that has a y-value on a line between the start and end of the 
		% spike
		
		try
			y2 = snip{s}(pend(s));
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
		catch
			fprintf('error in spike removal, aborting this spike...\n')
		end
		
	end		% END SPIKE LOOP (s)

	subplot(212)
	plot(din, 'g');
	hold on
	plot(dout);
	hold off
end




