function c=csucosmos
% Configures cluster in cosmos/LUNARC
%
% configCluster; This was only needed once
c=parcluster;
c.AdditionalProperties.AccountName = 'lu2024-2-6';
c.AdditionalProperties.WallTime  = '24:00:00';
c.AdditionalProperties.QueueName ='lu48';
c.AdditionalProperties.RequireExclusiveNode=true;
c.AdditionalProperties.ProcsPerNode=48;
c.saveProfile;
end