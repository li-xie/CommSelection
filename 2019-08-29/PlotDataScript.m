N=2000;
folder_name1='PlotData';
if ~exist(folder_name1,'dir')
    mkdir(folder_name1)
end
for i=7:9
    DataCollect(i,N);
%     filename1=[num2str(i) '/Data' num2str(i) '.mat'];
%     foldername='PlotData';
%     copyfile(filename1,foldername)
end