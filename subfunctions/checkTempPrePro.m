function boolTempPresent = checkTempPrePro(sRec)
	strTempFile = [sRec.sMD.strTemp sRec.strSession sRec.sProcLib.strRecording '.mat'];
	if exist(strTempFile,'file')
		boolTempPresent = true;
	else
		boolTempPresent = false;
	end
end