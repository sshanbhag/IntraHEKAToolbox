rawpath = '/Users/sshanbhag/Work/Data/IntracellularAmygdala/RawData/';
logfile = [rawpath 'IntracellDataLog.csv'];
dpath = [rawpath '2012-03-12'];
basename ='2012-03-12-5'; 
datname = [basename '.dat'];
dscpath = pwd;
pgfname = [basename '.pgf'];
pulname = [basename '.pul'];
dscname = [basename '.dsc'];


% read info from pgf file
pgf = readPGF(fullfile(dpath, pgfname));
% read info from pul file
pul = struct('tr', [], 'rr', []);
[pul.tr, pul.rr] = readPUL(fullfile(dpath, pulname));
% build dsc struct
dsc = buildDSC(fullfile(dpath, datname), pgf, pul.rr);

% test creation of sweep object
a = sweep;
a.Nsamples = dsc.size(1);
a.Rate = dsc.rate(1);
a.fp1 = dsc.trace1_fp(1);
a.fp2 = dsc.trace2_fp(1);
[trace1, trace2] = a.ReadData(fullfile(dpath, datname));

% test experiment object
e = experiment(dpath, basename, 'initialize')

% test get trace method for experiment class
[t1, t2] = e.GetTrace(1:e.Nsweeps);


% load log file
logData = readLog(logfile);

% set information in experiment object from logData
e.SetInfoFromLog(logData);

e.PlotOpts.Scale = 15;

%{
for C = 1:1
	figure

	% get the sweeps for the indicated Condition
	sinfo = e.GetStimulusForCondition(C);
	e.PlotTracesForCondition(C);
	titlestr =	{	e.BaseName, ...
						e.Info.Animal, ...
						[sinfo.AuditoryStimulus ' ' sinfo.OtherStimulus ' ' sinfo.Comments] ...
					};
	title(titlestr)
end
%}

[rate, trace1, trace2] = e.GetResampledTrace(1, 10);
stim = e.GetStimParamForCondition(1);

t = rate * ((1:length(trace1)) - 1);
plot(t, trace1)
line(stim.start .* [1 1], ylim, 'Color', 'g')
line(stim.end .* [1 1], ylim, 'Color', 'r')




