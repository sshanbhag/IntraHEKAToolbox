dpath = '/Users/sshanbhag/Work/Data/IntracellularAmygdala/RawData/2012-03-12';
dname = '2012-03-12-7.dat';
dname = '2012-03-12-9.dat';
pgfname = '2012-03-12-9.pgf';
pulname = '2012-03-12-9.pul';
dscname = '2012-03-12-9.dsc';

pgf = readPGF(fullfile(dpath, pgfname));

[pul_tr, pul_rr] = readPUL(fullfile(dpath, pulname));

d = readDAT(fullfile(dpath, dname), pgf, pul_rr);

[D, PGF, PUL, DSC] = readHEKA(fullfile(dpath, dname));

save(dscname, 'DSC', '-MAT');
