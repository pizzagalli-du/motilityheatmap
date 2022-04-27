xPath = 'D:\Temp\pppp.ims';
imsObj = ImarisReader(xPath);

xSize = imsObj.DataSet.SizeX;
ySize = imsObj.DataSet.SizeY;
zSize = imsObj.DataSet.SizeZ;

cIdx = 1; % Second channel index
tIdx = 1; % Second time point
dataVolume = imsObj.DataSet.GetDataVolume(cIdx, tIdx);

