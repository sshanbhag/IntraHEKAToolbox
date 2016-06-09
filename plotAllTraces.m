function aH = plotAllTraces(trialdata, varargin)
%-----------------------------------------------------------------------------
% aH = plotAllTraces(trialdata, varargin)
%-----------------------------------------------------------------------------
% 
% Given input data of all trial sweeps (trialdata), plots all traces
% overlaid
%-----------------------------------------------------------------------------
% Input Arguments:
%
% Output Arguments:
%
%-----------------------------------------------------------------------------
% See also: readPGF, readPUL, readHEKA, readDAT
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 4 April, 2012 (SJS)
%
% Revisions:
%	18 Apr 2012 (SJS): 
% 	 -	added PlotsPerColumn as input variable, but still need to 
% 		implement functionality to plot several columns of subplots
% 		to make datasets with large numbers of sweeps more visible
% 	19 Apr 2012 (SJS):
% 	 -	added Scale input argument to scale traces in plots by constant
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------


%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% Some settings/constants
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

PLOTSTIM = 0;
TVEC = 0;
COLORS = 1;
COLORLIST = { ...
	[0.00	0.00	1.00]	, ...
	[0.00	1.00	0.00]	, ...
	[1.00	0.00	0.00]	, ...
	[0.75	0.75	0.25]	, ...
	[1.00	0.00	1.00]	, ...
	[0.00	1.00	1.00]	, ...
	[0.00	0.75	0.75]	, ...
	[0.25	0.00	0.50]	, ...
	[0.00	0.75	0.00]	, ...
	[0.75	0.25	1.00]	, ...
	[0.00	0.00	0.75]	, ...
	[0.75	0.00	0.25] , ...
	[0.25	0.10	0.50]	, ...
	[0.10	0.50	0.50]	, ...
	[0.50	0.10	0.25]	, ...
};
SCALE = 1;
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% Parse inputs
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

optargin = size(varargin,2);
stdargin = nargin - optargin;

if stdargin ~= 1
	help plotAllTraces;
	aH = [];
	return
end

if optargin
	% parse varargin args
	errFlag = 0;
	n = 1;
	while n < optargin
		if ischar(varargin{n})
			if strcmpi(varargin{n}, 'Stimulus')
				PLOTSTIM = 1;
				stim = varargin{n+1};
				n = n + 2;
			elseif strcmpi(varargin{n}, 'SampleRate')
				Fs = varargin{n+1};
				TVEC = 1;
				n = n + 2;
			elseif strcmpi(varargin{n}, 'SampleInterval')
				Fs = 1/varargin{n+1};
				TVEC = 1;
				n = n + 2;
			elseif strcmpi(varargin{n}, 'Colors')
				if (n == optargin)
					COLORS = 1;
					n = n + 1;
				else
					COLORS = 1;
					if iscell(varargin{n + 1})
						COLORLIST = varargin{n + 1};
					else
						COLORLIST = {varargin{n + 1}};
					end
					n = n + 2;
				end
			elseif strcmp(varargin{n}, 'PlotsPerColumn')
				PlotsPerColumn = varargin{n+1};
				n = n + 1;
			elseif strcmp(varargin{n}, 'Scale')
				SCALE = varargin{n+1};
				n = n + 1;
			else
				error('%s: unknown option %s', mfilename, varargin{n});
			end
		else
			n = n + 1;
		end
	end
end

NCOLORS = length(COLORLIST);

ntrials = length(trialdata);

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% plot!
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

if ~exist('PlotsPerColumn', 'var')

	cIndex = 1;
	for p = 1:ntrials
		if cIndex > NCOLORS
			cIndex = 1;
		end
		if TVEC
			tvec = 1000 .* ((0:(length(trialdata{p})-1)) ./ Fs);
		else
			tvec = 1:length(trialdata{p});
		end

		if p ~= 1
			hold on
		end
		plot(tvec, (ntrials - (p-1)) + SCALE * trialdata{p}, 'Color', COLORLIST{cIndex});
		if p ~= 1
			hold off
		end
		cIndex = cIndex + 1;
	end

	if PLOTSTIM
		if TVEC
			tvec = 1000 .* ((0:(length(stim)-1)) ./ Fs);
		else
			tvec = 1:length(stim);
		end
		hold on
		plot(tvec, ntrials + 1.25 + (ntrials/2)*stim, 'k');
		hold off
	end

	if TVEC
		xlabel('Time (ms)');
	else
		xlabel('Sample');
	end

	ylim([0 ntrials + 2])
	set(gca, 'YTickLabel', [], 'YTick', [], 'Box', 'off');
	set(gca, 'YColor', get(gca, 'Color'));

	aH = gca;
	
	return
end

%%%% deal with columns... still needs to be done
%{

Ncols = fix(ntrials/PlotsPerColumn)


colIndex = 1;

cIndex = 1;
for p = 1:ntrials
	
	if cIndex > NCOLORS
		cIndex = 1;
	end
	if TVEC
		tvec = 1000 .* ((0:(length(trialdata{p})-1)) ./ Fs);
	else
		tvec = 1:length(trialdata{p});
	end

	if p ~= 1
		hold on
	end
	plot(tvec, (ntrials - (p-1)) + trialdata{p}, 'Color', COLORLIST{cIndex});
	if p ~= 1
		hold off
	end
	cIndex = cIndex + 1;
end

if PLOTSTIM
	if TVEC
		tvec = 1000 .* ((0:(length(stim)-1)) ./ Fs);
	else
		tvec = 1:length(stim);
	end
	hold on
	plot(tvec, ntrials + 1.25 + (ntrials/2)*stim, 'k');
	hold off
end

if TVEC
	xlabel('Time (ms)');
else
	xlabel('Sample');
end

ylim([0 ntrials + 2])
set(gca, 'YTickLabel', [], 'YTick', [], 'Box', 'off');
set(gca, 'YColor', get(gca, 'Color'));

aH = gca;


% 
%}
