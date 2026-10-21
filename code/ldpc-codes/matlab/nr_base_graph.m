function V = nr_base_graph(bg, ils)
%NR_BASE_GRAPH Base matrix of 5G NR BG1 (46 x 68) or BG2 (42 x 52) for the
% lifting-set index ils (0..7), read from the TS 38.212 tables (csv files).
if bg == 1, V = -ones(46, 68); else, V = -ones(42, 52); end
fid = fopen(fullfile(fileparts(mfilename('fullpath')), sprintf('5G_bg%d.csv', bg)));
fgetl(fid); fgetl(fid);                               % two header lines
row = 0;
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    f = strsplit(line, ';');
    if ~isempty(strtrim(f{1})), row = str2double(f{1}); end
    V(row + 1, str2double(f{2}) + 1) = str2double(f{3 + ils});
end
fclose(fid);
end
