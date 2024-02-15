function Collected_Scans = Collect_Scans(savef)
%Collect_Scans
%   Joins all scans results contained in a folder into a single structure and
%   saves it to a file
%   The folder is soecified in a dialogue
%   The concataned file is creatde in the active folder with an
%   automatically generated nam
%
%   Allows breaking up a long scan in multiple files and joinging them
%   later.
%

folder = uigetdir;
list=cellstr(ls(strcat(folder,'\*.mat')));
nfiles=size(list,1);

filename=strcat(folder,'\',list{1});
load(filename,'ScanResults');
Collected_Scans=ScanResults;
if not(isfield(Collected_Scans,'nerror'))
    Collected_Scans.nerror=0;
end
for i=2:nfiles
    filename=strcat(folder,'\',list{i});
    load(filename,'ScanResults');
    try
      Collected_Scans.nscan = Collected_Scans.nscan + ScanResults.nscan;
    catch
      fprintf('Error reading nscan in %s \n', filename);
    end
    try
      Collected_Scans.nstab = Collected_Scans.nstab + ScanResults.nstab;
    catch
      fprintf('Error reading nstab in %s \n', filename);
    end
    try
      Collected_Scans.nrecord = Collected_Scans.nrecord + ScanResults.nrecord;
    catch
       fprintf('Error reading nrecord in %s \n', filename);
    end
    try
         Collected_Scans.nerror = Collected_Scans.nerror + ScanResults.nerror;
    catch
         fprintf('Error reading nerror in %s \n', filename);
    end
    try
        Collected_Scans.Values = cat(1,Collected_Scans.Values,ScanResults.Values);
    catch
        fprintf('Error reading Values in %s \n', filename);
    end 
end

filename=strcat(folder,'_Collected');
if (strcmp(savef,'Y'))
       save(filename,'Collected_Scans');
end

end
