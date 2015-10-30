function output = convertFolder2avi(strPath)

	%get/set folder info
	newFolder = 'G:\Programs\WinFF\';
	oldFolder = cd(newFolder);

	%retrieve folder information
	if ~strcmp(strPath(end),filesep)
		strPath(end+1) = filesep;
	end
	structDir = dir(strPath);
	intFiles = length(structDir);
	
	%check subfolder availability
	strSubfolder = ['avi' filesep];
	cd([strPath strSubfolder]);
	cd(newFolder);
	
	%supported formats
	cellSupported{1} = 'avi';
	cellSupported{2} = 'mpg';

	%encoding options
	codec{1} = 'h264';
	codec{2} = 'mpeg2video';
	codec{3} = 'libx264';
	useCodec = 3;
	strBitrate = '1200k';
	fps = 25;

	%loop through files
	for intFile = 1:intFiles
		strFile = structDir(intFile).name;
		extLoc = strfind(strFile,'.');
		if isempty(extLoc)
			continue;
		end
		strExt = strFile((extLoc(end)+1):end);
		if max(strcmp(strExt,cellSupported)) == 1
			strInput = [strPath strFile];
			strOut = [strFile(1:(extLoc(end)-1)) codec{useCodec} '.avi'];

			strOutput = sprintf('%s%s%s',strPath,strSubfolder,strOut);

			%mpeg to avi
			strExpression = ['ffmpeg -i ' strInput ' -an -sn -r ' num2str(fps) ' -b:v ' strBitrate ' -vcodec ' codec{useCodec} ' ' strOutput];
			
			fprintf('Executing command (%d of %d): %s\n',intFile,intFiles,strExpression)
			dos(strExpression)
			fprintf('File %s succesfully converted to %s!\n',strFile,strOut)
		end
	end
	cd(oldFolder);
	output = status;
end
%{


rootDir='G:\Jorrit\data\J20111104\xyt01';

TPD=tpd_newDataStruct;

tpd_initializeSession(TPD, rootDir)

%BATCH_Timeseries2Session
%}