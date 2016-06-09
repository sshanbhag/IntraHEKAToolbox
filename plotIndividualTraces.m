function plotH = plotIndividualTraces(trialdata, varargin)
%-----------------------------------------------------------------------------
% aH = plotIndividualTraces(trialdata, varargin)
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
	help plotIndividualTraces;
	plotH = [];
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

HGAP = 0.05;
VGAP = 0.00;
PGAP = 0.0125;

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
			elseif strcmpi(varargin{n}, 'Hgap')
				HGAP = varargin{n + 1};
				n = n + 2;
			elseif strcmpi(varargin{n}, 'Vgap')
				VGAP = varargin{n + 1};
				n = n + 2;
			elseif strcmpi(varargin{n}, 'Pgap')
				PGAP = varargin{n + 1};
				n = n + 2;
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


Ncols = 1;
Nrows = ntrials;
plotwidth = (1 - ((Ncols+1) * HGAP)) / Ncols;
plotheight = (1 - ((Nrows+1) * VGAP)) / (Nrows);

pos1 = cell(Nrows, Ncols);
for rr = 1:Nrows
	ypos1(rr) = 1 - rr*plotheight - (rr-1)*VGAP;
	for cc = 1:Ncols
		xpos(cc) = HGAP;
		pos1{rr, cc} = [xpos(cc) ypos1(rr) plotwidth plotheight];
	end
end

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
	
	plotH(p) = subplot('Position', pos1{p, 1});
	plot(tvec, trialdata{p}, 'Color', COLORLIST{cIndex});
	if (p == 1) && PLOTSTIM
		if TVEC
			tvec = 1000 .* ((0:(length(stim)-1)) ./ Fs);
		else
			tvec = 1:length(stim);
		end		
		hold on
		plot(tvec, stim, 'k');
		hold off
	end
	
	set(gca, 'XTickLabel', [], 'YTickLabel', []);
	set(gca, 'Box', 'off');
	set(gca, 'Visible', 'off');

	cIndex = cIndex + 1;
end



