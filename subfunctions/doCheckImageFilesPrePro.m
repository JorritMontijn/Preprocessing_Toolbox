function sImCheck = doCheckImageFilesPrePro(sRec)
% - load and check image consistency
% - - make universal reference library to images

	intChannels = 0;
	boolMoreChannels = true;
	intLastNumel = -1;
	strPath = [sRec.sMD.strSourceDir sRec.sMD.strImgSource sRec.strSession filesep sRec.sRawLib.strRecording filesep];
	while boolMoreChannels
		sImages = dir([strPath sRec.sRawLib.strName '_t*ch' sprintf('%02d',intChannels) '.tif']);
		intThisNumel = numel(sImages);
		if intThisNumel == 0
			boolMoreChannels = false;
		else
			if intLastNumel ~= -1
				if intLastNumel ~= intThisNumel
					warning([mfilename ':ChNumelError'],'Number of elements for all channels are not consistent: files must be missing! Trying to find which one(s) are missing..');
					boolFound=false;
					for intCh=0:intChannels
						sImages = dir([strPath sRec.sRawLib.strName '_t*ch' sprintf('%02d',intCh) '.tif']);
						fprintf('Starting search through channel %d; number of elements is %d\n',intCh,numel(sImages));
						for intT=0:(max(intThisNumel,intLastNumel)-1)
							intFileT = str2double(getFlankedBy(sImages(intT+1).name,'_t','_ch'));
							if intFileT ~= intT
								fprintf('Element %d of channel %d is inconsistent with filename %s; image %d must be missing!\n',intT,intCh,sImages(intT+1).name,intT);
								boolFound=true;
								break;
							end
						end
						if ~boolFound,fprintf('No missing files in channel %d!\n',intCh);end
					end
					if ~boolFound
						fprintf('Sorry... Could not find anything wrong:s\n');
					end
					error([mfilename ':ChNumelError'],'Please find the missing files...');
				end
			end
			intLastNumel = intThisNumel;
			intChannels = intChannels + 1;
		end
	end

	
	%find t-range
	sImages = dir([strPath sRec.sRawLib.strName '_t*ch' sprintf('%02d',0) '.tif']);
	intThisNumel = numel(sImages);
	
	strFileFirst = sImages(1).name;
	intStartT = strfind(strFileFirst,'_t') + 2;
	vecSplits = strfind(strFileFirst,'_');
	intStopT = vecSplits(find(vecSplits > intStartT,1,'first'))- 1;
	intLengthT = intStopT - intStartT + 1;

	imFirst = imread([strPath sRec.sRawLib.strName '_t' sprintf(['%0' num2str(intLengthT) 'd'],0) '_ch' sprintf('%02d',intChannels - 1) '.tif']);
	imLast = imread([strPath sRec.sRawLib.strName '_t' sprintf(['%0' num2str(intLengthT) 'd'],intThisNumel - 1) '_ch' sprintf('%02d',intChannels - 1) '.tif']);

	%check im sizes
	sizeX = size(imFirst,2);
	sizeY = size(imFirst,1);
	if min(size(imFirst) == size(imLast)) == 0
		error([mfilename ':ChNumelError'],'Number of elements in all channels is inconsistent: files must be missing!');
	end

	sImCheck.t = intThisNumel;
	sImCheck.z = 1;
	sImCheck.x = sizeX;
	sImCheck.y = sizeY;
	sImCheck.ch = intChannels;
end