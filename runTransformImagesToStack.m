clear all
%% input
if ~exist('strSession','var')
	strSession = '20130627';
	vecRecordings = 1;
	strDir{1} = 'D:\Data\Raw\imagingdata\20130627\xyt01';
	%strDir{2} = 'D:\Data\Raw\imagingdata\20120718\OT2';
end

%% define general metadata
sPS = loadDefaultSettingsPrePro();%structProcessingSettings
strMasterDir = 'D:\Data';
strTargetDir = '\Raw\imagingdata\';

%% create filenames
strMasterPath = [strMasterDir strTargetDir strSession filesep];
strOldPath = cd(strMasterPath);
for intRec=vecRecordings
	strRec{intRec} = [sPS.strRecording{intRec}]; %masks
end

%% run
for intRec=vecRecordings

	%get path
	strPath = strDir{intRec};
	cd(strPath);
	
	%msg
	fprintf('Starting image to stack transformation for [%s]\n',strPath);
	
	%get stack info
	sFiles=dir('*t*_ch*.tif*');
	strLast=sFiles(end).name;
	
	%get header
	strHeader = getFlankedBy(strLast,'','_t');
	if ~isempty(strHeader)
		strHeader=[strHeader '_'];
		strT = '_t';
	else
		strT = 't';
	end
	
	%get t-size
	strMaxT=getFlankedBy(strLast,strT,'_ch');
	intLengthT=length(strMaxT);
	intMaxT = str2double(strMaxT);
	
	%get ch-size
	strMaxCh=getFlankedBy(strLast,'_ch','.');
	intLengthCh=length(strMaxCh);
	intMaxCh = str2double(strMaxCh);
	
	%build im format
	strIm=[strHeader 't%0' num2str(intLengthT) '.0f_ch%0' num2str(intLengthCh) '.0f.tif'];

	%load image
	imLast=imread(strLast);
	
	%get image properties
	intX=size(imLast,2);
	intY=size(imLast,1);
	strClass=class(imLast);
	
	%pre-allocate
	matStack = zeros([intY intX intMaxCh+1 intMaxT+1],strClass);
	ptrDataColl = tic;
	
	%transform
	for intT=0:intMaxT
		for intCh=0:intMaxCh
			strImTemp=sprintf(strIm,intT,intCh);
			matStack(:,:,intCh+1,intT+1) = imread(strImTemp);
		end
		if mod(intT,1000) == 0
			fprintf('Now at image %d of %d\n',intT,intMaxT);
		end
	end
	
	%save data
	dblImLoadTime=toc(ptrDataColl);
	fprintf('Data collection completed [%.1fs]; saving...\n',dblImLoadTime); 
	ptrSaveTime = tic;
	strSaveFileName = ['imStack' strSession strRec{intRec}];
	save(strSaveFileName,'matStack','-v7.3')
	
	%msg
	sInfo = dir([strSaveFileName '*']);
	fprintf('Transformation completed; data saved [%.1fs] to %s [size: %.1f MB; dir: %s]\n',toc(ptrSaveTime),strSaveFileName,sInfo(1).bytes / 10^6,strPath);
	
	clear matStack;
	ptrLoad = tic;
	fprintf('Test loading image stack...\n');
	load([strSaveFileName '.mat'])
	dblStackLoadTime=toc(ptrLoad);
	fprintf('\b		Done! Loading took %.1f seconds [separate image load took %.1f seconds]',dblStackLoadTime,dblImLoadTime);
end
cd(strOldPath);