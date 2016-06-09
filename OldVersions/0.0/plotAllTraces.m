function aH = plotAllTraces(trialdata, varargin)
%-----------------------------------------------------------------------------
% aH = plotAllTraces(stim, trialdata, varargin)
%-----------------------------------------------------------------------------
% 
% 
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
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% Read information from file
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------

optargin = size(varargin,2);
stdargin = nargin - optargin;

if stdargin ~= 1
	help plotAllTraces;
	aH = [];
	return
end

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


