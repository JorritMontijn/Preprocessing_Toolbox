function sRegZ = doAcrossRecordingDriftCheck(sSourceData)
	%doDriftCheck Performs z-stack registration of a time series recording
	%   Outputs a similarity measure per z plane per time point
	%
	%Syntax:
	%sRegZ = doDriftCheck(sSourceData)
	%
	%Input must be a structure with the following fields:
	%- strMasterPath; path of master directory (session)
	%- strDirRawStackZ; path containing z-stack images
	%- strSession; name of session (date)
	%- vecRecordings; name of recording (xyt01 ... xyt99)
	
	%get source data
	strMasterPath = sSourceData.strMasterPath;
	strDirRawStackZ = sSourceData.strDirRawStackZ;
	strSession = sSourceData.strSession;
	vecRecordings = sSourceData.vecRecordings;
	
	%calc paths
	strOldPath = cd();
	strSesPath = [strMasterPath strSession];
	vecFileSeps = strfind(strDirRawStackZ(1:(end-1)),filesep);
	strDirZ = strtok(strDirRawStackZ((vecFileSeps(end)+1):end),filesep);
	
	%parameters
	intUseChannelSwitch = 2; %1 only red; 2 both
	intSampleLength = 100;
	intSubPixelDepth = 1;
	
	%msg
	fprintf('Loading data and preparing across recording Z-stack registration...\n')
	
	%% load prepro file of first recording
	strRecording = sprintf('xyt%02d',vecRecordings(1));
	load([strSesPath filesep strRecording filesep strSession strRecording '_prepro.mat']);
	
	%% prepare z-stack
	cd(strDirRawStackZ);
	sDirZ = dir('*z*_ch*.tif');
	if intUseChannelSwitch == 2
		matStackZ = nan(sRec.sProcLib.y,sRec.sProcLib.x,numel(sDirZ)/2,2);
	else
		matStackZ = nan(sRec.sProcLib.y,sRec.sProcLib.x,numel(sDirZ));
	end
	
	%% load z-stack
	%get raw data
	for intFile=1:numel(sDirZ)
		%get info
		strFile = sDirZ(intFile).name;
		intZ = str2double(getFlankedBy(strFile,'_z','_ch'))+1;
		if intUseChannelSwitch == 2
			intCh = str2double(getFlankedBy(strFile,'_ch','.tif'))+1;
		else
			intCh = 1;
		end
		
		%load image & put in stack
		imTemp = imread(strFile);
		matStackZ(:,:,intZ,intCh) = im2double(imTemp);
	end
	matStackZ = mean(matStackZ,4);
	intSizeStackZ = size(matStackZ,3);
	
	%make fourier-transformed stack
	matStackZ_fft = nan(size(matStackZ));
	for intZ=1:intSizeStackZ
		matStackZ_fft(:,:,intZ) = fft2(matStackZ(:,:,intZ));
	end
	
	%get z-stack slice interval from xml file
	strStackZXML = [strDirRawStackZ filesep 'MetaData' filesep strDirZ '_Properties.xml'];
	sData = loadXMLPrePro(strStackZXML);
	sData.dblMicronPerPlane = abs(sData.dblActualImageSizeZ)/intSizeStackZ;
	
	%pre-allocate
	matRegistrationZ = nan(length(vecRecordings),3,intSizeStackZ,4);
		
	%% loop through images & compare to z-stack
	%pre-load image buffer
	for intRecording=1:length(vecRecordings)
		%define recording
		strRecording = sprintf('xyt%02d',vecRecordings(intRecording));
		
		%get recording data
		load([strSesPath filesep strRecording filesep strSession strRecording '_prepro.mat']);
		intMaxT = sRec.sProcLib.t;
		intSizeX = sRec.sProcLib.x;
		intSizeY = sRec.sProcLib.y;
		intMaxCh = sRec.sProcLib.ch;
		intMinCh = 1;
		intLengthT = length(num2str(intMaxT));
		if intUseChannelSwitch == 1
			intMaxCh = 1;
		end
		
		%define image locations
		strTargetIm = ['t%0' num2str(intLengthT) 'd_ch%02d.tif'];
		
		for intPointT=1:3
			matSample = nan(intSizeY,intSizeX,intMaxCh,intSampleLength);
			if intPointT == 1
				intStartT = 100;
			elseif intPointT == 2
				intStartT = floor(intMaxT/2 - intSampleLength/2);
			elseif intPointT == 3
				intStartT = intMaxT - intSampleLength - 100;
			end
			for intT = 1:intSampleLength
				for intCh=intMinCh:intMaxCh
					strIm = sprintf(strTargetIm,intT+intStartT-1,intCh-1);
					
					%transform to 2D to save space
					imThis = imread([strSesPath filesep strRecording filesep 'images' filesep strIm]);
					if ndims(imThis) == 3
						imThis = mean(imThis,3);
					end
					matIm2D = im2double(imThis);
					
					%put into buffer
					matSample(:,:,intCh,intT) = matIm2D;
				end
			end
			
			%make mean image from sample for registration
			matRegIm = mean(mean(matSample,4),3);
			matReg_fft = fft2(matRegIm);
			
			%loop through z-stack to register images
			for intZ=1:intSizeStackZ
				%calculate required registration
				matRegistrationZ(intRecording,intPointT,intZ,:) = dftregistration(matStackZ_fft(:,:,intZ), matReg_fft, intSubPixelDepth);
			end
		end
		fprintf('Finished recording %s...\n',strRecording)
	end
	
	%make structure
	sRegZ = struct;
	sRegZ.matRegistrationZ = matRegistrationZ;
	sRegZ.sData = sData;
	sRegZ.strSession = strSession;
	sRegZ.vecRecordings = vecRecordings;
	sRegZ.sRec = sRec;
	
	%save data
	strFile = [strSesPath filesep strSession '_across_rec_zreg.mat'];
	save(strFile,'matRegistrationZ','sData','strSession','vecRecordings','sRec');
	fprintf('Saved recording structure to %s; time is [%s]\n',strFile,getTime)
	cd(strOldPath);
end

