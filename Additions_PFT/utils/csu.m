function c=csu
% Startand configure cluster object
%configCluster; This is only needed once.
c=parcluster;
% In Aurora (to which accdev-0 is directly connected), choose Aurora
% R2002a as the default cluster before running this script
% The following is for use in cosmos (LUNARC)
%c.AdditionalProperties.AccountName = 'lu2024-2-6';

%c.AdditionalProperties.QueueName ='lu48';
%c.AdditionalProperties.RequireExclusiveNode=true;
%c.AdditionalProperties.ProcsPerNode=48;

% The following applies to both LUNARC cosmos and MAX IV Aurora clusters
c.AdditionalProperties.WallTime  = '24:00:00';
c.AdditionalProperties.QueueName ='all';

% parpool(c,16*n) % This is to be run on the command line everytime the
% parpool needs to be started
% delete(gcp('nocreate')) - used to stopthe parpool.
end