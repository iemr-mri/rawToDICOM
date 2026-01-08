function [SAXnames, LAXnames] = scanNameSplitter(folderStructs)
% Detects SAX or LAX and makes two lists
SAXnames = [];
LAXnames = [];
for scan = 1:length(folderStructs)
    scanName = folderStructs(scan).name;
    if contains(scanName,'SAX')
        SAXnames = [SAXnames, string(scanName)];
    elseif contains(scanName, 'LAX')
        LAXnames = [LAXnames, string(scanName)];
    else
        warning('SG Pipeline: Scan name did not contain "SAX" or "LAX"');
    end
end