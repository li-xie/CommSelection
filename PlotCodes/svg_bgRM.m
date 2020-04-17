function svg_bgRM(filePath,cut_file_path)
% Full path of the sample text file
% filePath = fullfile('Copy.svg');
% fid = fopen('svg_cut.svg');
fid = fopen(cut_file_path);
cut_str = fscanf(fid,'%c');
fclose(fid);
% Read the file
fid = fopen(filePath,'r');
str = fscanf(fid,'%c');
fclose(fid);
% get rid of the background objects
str2 = erase(str, cut_str);
% Save as a text file
fid2 = fopen(filePath,'w');
fprintf(fid2,'%s\n', str2);
fclose(fid2);