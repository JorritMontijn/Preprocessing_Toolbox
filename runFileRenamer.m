%function runFileRenamer(strDir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% switch dirs
if exist('strDir','var'),cd(strDir);end
strHomeDir = cd();

%get subdirs
sDir = dir();

%loop through directories
for intSubDir=5%3:numel(sDir)
	%skip if not dir
	if ~sDir(intSubDir).isdir,continue;end
	cd(strHomeDir);
	
	%get old name
	strSubDir = sDir(intSubDir).name;
	sSubName=dir([strSubDir filesep 'MetaData\*Properties.xml']);
	strOldName = sSubName(1).name(1:(end-15));
	intOldLength = length(strOldName);
	
	%get all files in parent folder & rename
	cd([strHomeDir filesep strSubDir]);
	sRename=dir([strOldName '*']);
	fprintf('\nRenaming "%s" to "%s"\n',strOldName,strSubDir);
	for intFile=1:numel(sRename)
		strNewFileName = strrep(sRename(intFile).name,strOldName,strSubDir);
		strOldFileName = sRename(intFile).name;
		if strcmpi(strNewFileName,sRename(intFile).name),continue;end
		if mod(intFile,1000) == 0,fprintf('Renaming %s [%d of %d]\n',sRename(intFile).name,intFile,numel(sRename));end
		
		%rename file
		objJava = java.io.File(strOldFileName);
		objJava.renameTo(java.io.File(strNewFileName));
	end
	
	%get all files in MetaData folder & rename
	cd([strHomeDir filesep strSubDir filesep 'MetaData'])
	sRename=dir([strOldName '*']);
	for intFile=1:numel(sRename)
		strOldFileName = sRename(intFile).name;
		strNewFileName = strrep(sRename(intFile).name,strOldName,strSubDir);
		if strcmpi(strNewFileName,sRename(intFile).name),continue;end
		
		%rename file
		objJava = java.io.File(strOldFileName);
		objJava.renameTo(java.io.File(strNewFileName));
	end
	
	fprintf('Renaming "%s" to "%s completed!"\n',strOldName,strSubDir);
end
cd(strHomeDir);

%end

