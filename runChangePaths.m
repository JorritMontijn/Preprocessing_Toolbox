%Changing paths to D:\
drive='D';
list=dir([drive ':\Data\Processed\imagingdata\']);
rec_date='20140314';
rec_list=dir([drive ':\Data\Processed\imagingdata\' rec_date '\']);
[num_ses s]=size(rec_list);

for i=3:num_ses
	ses_name=rec_list(i,1).name;
	fprintf('Changing file [%s%s]\n',rec_date,ses_name);
	ses_list=dir([drive ':\Data\Processed\imagingdata\' rec_date '\' ses_name '\']);
	load([drive ':\Data\Processed\imagingdata\' rec_date '\' ses_name '\' rec_date ses_name '_CD.mat'])
	sRec.sMD.strMasterDir(1)=drive;
	if isfield(sRec,'sDC')
		sRec.sDC.strRecPath(1)=drive;
	end
	save([drive ':\Data\Processed\imagingdata\' rec_date '\' ses_name '\' rec_date ses_name '_CD.mat'])
end