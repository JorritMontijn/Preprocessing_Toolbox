strDir = 'D:\Data\Processed\imagingdata\';
sDir= dir(strDir);
%loop through days
for intSes=1:numel(sDir)
	strSes = sDir(intSes).name;
	if length(strSes) == 8 && str2double(strSes) > 20120101
		fprintf('%s; checking for recordings\n',strSes)
		
		%loop through recordings
		strSubDir = [strDir strSes];
		sSubDir= dir(strSubDir);
		for intRec=1:numel(sSubDir)
			strRec = sSubDir(intRec).name;
			if length(strRec) == 5 && strcmp(strRec(1:3),'xyt')
				fprintf('	%s; checking for session files\n',strRec)
				
				
				%search for file
				strRecDir = [strDir strSes filesep strRec];
				strSesFile = sprintf('%s%s_ses.mat',strSes,strRec);
				strSesPtr = [strRecDir filesep strSesFile];
				strPreProFile = sprintf('%s%s_prepro.mat',strSes,strRec);
				strPreProPtr = [strRecDir filesep strPreProFile];
				if exist(strSesPtr,'file') && exist(strPreProPtr,'file')
					fprintf('		%s and %s present; updating ses file...\n',strSesFile,strPreProFile);
					
					%load files and check if it needs to be updated
					clear ses;
					load(strSesPtr);
					if isfield(ses,'xml')
						fprintf('\b		Already updated!\n')
					else
						clear sRec;
						load(strPreProPtr);
						if sRec.sProcLog.boolXMLFound
							ses.xml = sRec.xml.sData;
							ses.date = sRec.xml.sData.strStartTime;
							
							%save file
							save(strSesPtr,'ses')
							fprintf('\b		Added XML data to ses file!\n')
						else
							fprintf('\b		No XML data present in prepro file!\n')
						end
					end
				end
			end
		end
	end
end
fprintf('Done!\n')