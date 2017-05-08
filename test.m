y = -1 * [40 40 40 40 40];
x = (1:length(y));
errlohi = [5 10 15 20 25];

[x y errlohi]

figure(10)


%------------------------------------------------------------------------
% force x and y to column vector
%------------------------------------------------------------------------
x = force_col(x);
y = force_col(y);
nx = length(x);
ny = length(y);
%------------------------------------------------------------------------
% convert errlohi to column vector or matrix if needed
% assume if ncols in errlohi == 2, that it is ok
%------------------------------------------------------------------------
[nrows, ncols] = size(errlohi); 
% if ncols > 2, force to column vector using transpose
if ncols > 2
	errlohi = errlohi';
elseif ncols == 1
	errlohi = force_col(errlohi);
end
%------------------------------------------------------------------------
% make sure lengths match
%------------------------------------------------------------------------
nerr = length(errlohi);
if nerr ~= nx || nerr ~= ny || nx ~= ny
	error('%s: x, y, errlohi must be same length', mfilename);
end

%------------------------------------------------------------------------
% matrix for error area
%------------------------------------------------------------------------
if isvector(errlohi)
	errcoords = [ (y-errlohi) (y+errlohi)];
elseif ismatrix(errlohi)
	errcoords = [ (y-errlohi(:, 1)) (y+errlohi(:, 2))];
end



% plot area
ha = area(x, errcoords);
set(ha(1), 'LineStyle', 'none');
set(ha(1), 'FaceColor', 'none');
set(ha(2), 'LineStyle', 'none');
set(ha(2), 'FaceColor', 0.75 * [1 1 1]);

% plot line
hold on
hb = plot(x, y);
hold off

figure(11)
shadedErrorBar(x, y, errlohi,  {'MarkerFaceColor', [0 0.4470 0.7410]})