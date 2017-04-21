function in_data = readDAT(datfile, pgf, pul_rr)
%-----------------------------------------------------------------------------
% in_data = readDAT(datfile, pgf, pul)
%-----------------------------------------------------------------------------
% 
% Reads data from HEKA .dat file
% 
%-----------------------------------------------------------------------------
% Input Arguments:
% 	datfile		path/filename
%	pgf			(from readPGF)
% 	pul_rr		(from readPUL)
%
% Output Arguments:
% 	in_data		trace data structure
%
%-----------------------------------------------------------------------------
% See also: readPGF, readPUL, readHEKA
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 22 March, 2012 (SJS)
% 				based on Visual Basic code clsOData.vb by Olga Galazyuk
%
% Revisions:
%	27 Mar 2012 (SJS)
% 		-	added file location information for each sweep
% 		-	removed id substruct to simplify trace storage
%	9 Jun 2016 (SJS): updated comments, added check for failed fopen()
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% Read information from file
%-----------------------------------------------------------------------------
%-----------------------------------------------------------------------------
% open file
fp = fopen(datfile, 'r', 'ieee-le');

% check if error opening file
if fp == -1
	error('%s: could not open file %s for reading', mfilename, pgffile);
end

ns = 1;
for i = 1:pul_rr.number_children

	% for each series
	for j = 1:pul_rr.gr(i).number_children

		% for each sweep
		 for k = 1:pul_rr.gr(i).srr(j).number_children

			in_data(ns).size = pul_rr.gr(i).srr(j).swr(k).total_points; %#ok<*AGROW>
			in_data(ns).sweep = i + j + k;
			in_data(ns).rate = pgf.smr(j).sample_interval;
			in_data(ns).fname = datfile;

			% read the sound data first
			in_data(ns).trace1_fp = ftell(fp);
			in_data(ns).trace1 = fread(fp, in_data(ns).size, 'int16') ./ 1000000;

			% read the data and scale using data_factor2 (??)
			in_data(ns).trace2_fp = ftell(fp);
			tmp = fread(fp, in_data(ns).size, 'int16');
			in_data(ns).trace2 = tmp .* pul_rr.gr(i).srr(j).swr(k).data_factor2;
			
			ns = ns + 1;
		 end
	end
end

fclose(fp);


