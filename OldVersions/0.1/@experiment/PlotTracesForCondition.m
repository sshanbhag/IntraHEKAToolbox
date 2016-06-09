function h = PlotTracesForCondition(obj, Condition, varargin)
%-----------------------------------------------------------------------------
% experiment.PlotTracesForCondition(Condition, Decifactor, fHandle)
%-----------------------------------------------------------------------------
% experiment method
%-----------------------------------------------------------------------------
% 
% plot traces for condition
%
%-----------------------------------------------------------------------------
% Input Arguments:
% 	Condition			Condition id #
% 
% 	Optional Inputs:
% 		Decifactor			Decimation factor for traces (default is 1)
% 		fHandle				Figure handle for plot (default is current figure)
% 
% Output Arguments:
% 	h						figure handle
%
%-----------------------------------------------------------------------------
% See also: @experiment class definitions
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 14 April, 2012 (SJS)
%
% Revisions:
%	20 Apr 2012 (SJS):
% 	 -	added comments, Decifactor, varargin format
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------


% make sure object is initialized
if ~obj.isInitialized
	warning('%s: object not initialized', mfilename);
	return
end


% check inputs
switch length(varargin)
	case 0
		Decifactor = 1;
		figure
	case 1
		Decifactor = varargin{1};
		figure
	case 2
		Decifactor = varargin{1};
		fHandle = varargin{2};
		figure(fHandle)
	otherwise
		error('%s: strange inputs... ', mfilename);
end

% check that Condition is inbounds
if ~between(Condition, 1, obj.Info.Nconditions)
	warning('%s: Condition %d not within bounds (Nconditions: %d)', ...
					mfilename, Condition, obj.Info.Nconditions);
	return
end

% get the sweeps for the indicated Condition
sinfo = obj.Info.Stimulus(Condition);		
sweeplist = obj.GetSweepListForCondition(Condition);
[sampleinterval, t1, t2] = obj.GetResampledTrace(sweeplist, Decifactor);
htmp = plotAllTraces(	t2, ...
								'Stimulus', 0.1*normalize(t1{1}), ...
								'SampleInterval', sampleinterval, ...
								'Colors', obj.PlotOpts.Color, ...
								'Scale', obj.PlotOpts.Scale);
titlestr =	{	obj.BaseName, ...
					obj.Info.Animal, ...
					[sinfo.AuditoryStimulus ' ' sinfo.OtherStimulus] ...
				};
title(titlestr)

if nargout
	h = htmp;
end

