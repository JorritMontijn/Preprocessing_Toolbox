strPath = 'D:\Data\Processed\imagingdata\';
strSes = '20140207';
vecRecordings = 1:8;
intRecalcType = 3;
%{ 
recalc types:
1) rectify non-positive values to 0
2) set outliers to nan
3) dF/F without neuropil subtraction (only somatic)
4) dF/F from neuropil annulus
5) dF/F with neuropil subtraction
6) remove neurons (with extreme dF/F values or supplied list)
%}

for intRec=vecRecordings
	%load file
	strDir = sprintf('%s%s%sxyt%02d%s',strPath,strSes,filesep,intRec,filesep);
	cd(strDir);
	strFile = sprintf('%sxyt%02d_ses',strSes,intRec);
	fprintf('Processing [%s]...\n',strFile);
	sLoad = load(strFile);
	ses = sLoad.ses;
	[ses,dummy] = doRecalcdFoF(ses,intRecalcType);
	
	%save
	save([strFile '1'],'ses') %save temp
	movefile([strFile '.mat'],[strFile '.backup'],'f');
	movefile([strFile '1.mat'],[strFile '.mat'],'f');
	fprintf('Saved recording structure to %s\n',strFile)
end