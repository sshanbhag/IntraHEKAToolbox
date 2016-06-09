function dsc = buildDSC(datfile, pgf, pul_rr)
%-----------------------------------------------------------------------------
% dsc = buildDSC(datfile, pgf, pul_rr)
%-----------------------------------------------------------------------------
% 
% Reads data from HEKA .dat file, generates dsc data struct
% 
%-----------------------------------------------------------------------------
% Input Arguments:
% 	datfile		path/filename
%	pgf			(from readPGF)
% 	pul_rr		(from readPUL)
%
% Output Arguments:
% 	dsc			data description struct 
%
%-----------------------------------------------------------------------------
% See also: readPGF, readPUL, readHEKA, readDAT
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 4 April, 2012 (SJS)
% 				based on Visual Basic code clsOData.vb by Olga Galazyuk
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

dsc.fname = datfile;

% open file
fp = fopen(datfile, 'r', 'ieee-le');

ns = 0;
for i = 1:pul_rr.number_children

	% for each series
	for j = 1:pul_rr.gr(i).number_children

		% for each sweep
		 for k = 1:pul_rr.gr(i).srr(j).number_children
			ns = ns + 1;

			dsc.size(ns) = pul_rr.gr(i).srr(j).swr(k).total_points;
			dsc.sweep_index(ns) = i + j + k;
			dsc.rate(ns) = pgf.smr(j).sample_interval;
			dfactor1(ns) = 1 ./ 1000000;
			dfactor2(ns) = pul_rr.gr(i).srr(j).swr(k).data_factor2;

			% read the sound data (trace1) first
			dsc.trace1_fp(ns) = ftell(fp);
			trace1 = fread(fp, dsc.size(ns), 'int16');

			% read the data (trace2)
			dsc.trace2_fp(ns) = ftell(fp);
			trace2 = fread(fp, dsc.size(ns), 'int16');
			
		 end
	end
end

fclose(fp);

dsc.dfactor1 = unique(dfactor1);
dsc.dfactor2 = unique(dfactor2);

if length(dsc.dfactor1) > 1
	warning('%s: multiple dfactor1 values in file %s', mfilename, datfile);
	fprintf('dfactor1 = \n')
	for n = 1:length(dsc.dfactor1)
		fprintf('\t%d\t\t%f\n', n, dsc.dfactor1(n));
	end
end

if length(dsc.dfactor2) > 1
	warning('%s: multiple dfactor2 values in file %s', mfilename, datfile);
	fprintf('dfactor2 = \n')
	for n = 1:length(dsc.dfactor2)
		fprintf('\t%d\t\t%f\n', n, dsc.dfactor2(n));
	end
end

dsc.nsweeps = ns;




