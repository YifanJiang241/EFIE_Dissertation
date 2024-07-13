function [FSL,PL,H1,H1_60,fnTX,fnRX] = SPLAT_PL(TX,RX,a3)
% [FSL,PL,H1,H1_60,fnTX,fnRX] = SPLAT_PL(TX,RX,a3)
% function SPLAT_PL creates transmitter and receiver sites defined in TX
% and RX, respectively, then executes the SPLAT! code to calculate the path
% loss between the transmitter and receiver sites, and reads the resultant
% data files. 
% 
% The TX and RX can be either 
% -defined by pre-existing files (then the field 'FileName' must be set) or
% -by setting the fields 'SiteName', 'Lat', 'Lon', and 'AntHeightAGL'.
% 
% One can also specify additional command line arguments for Splat! using
% the char variable a3 (but Splat! does not always consider it - please 
% test your specific case prior to running a batch).
%
% One may need to set the paths to the Splat! executable (variable
% PathSplat), topology files (variable PathSDF), and TX/RX files (variable
% Path) in lines 35-38 of the code. 
% 
% Outputs:
% - FSL and PL are path loss calculations using free space and applying 
% topology, respectively
% - H1 and H1_60 are Splat!-computed recommended height to clear entire 
% first Fresnel zone or a practically necessary 60% of it. 
% - fnTX and fnRX are file names for the files created.
% 
% NB! The topology files for the region of interest must be available.
% Otherwise, Splat! makes some assumptions and does not take topology into
% accont. 
%
% A sample of usage is shown at the end of this function, in the TEST
% section.
% 
% The code was tested in Matlab 8.5 and Splat 1.3.1. 
% 
% Copyright by Albert Lysko, 2016-08-11
% Extract from Help for Splat!:
% SPLAT! imports site location information of transmitter and receiver sites analyzed by the program from
% ASCII files having a .qth extension. QTH files contain the site�s name, the site�s latitude (positive if North
% of the equator, negative if South), the site�s longitude (in degrees West, 0 to 360 degrees, or degrees East 0
% to -360 degrees), and the site�s antenna height above ground level (AGL), each separated by a single linefeed
% character. The antenna height is assumed to be specified in feet unless followed by the letter m or the
% word meters in either upper or lower case. Latitude and longitude information may be expressed in either
% decimal format (74.6864) or degree, minute, second (DMS) format (74 41 11.0).
narginchk(2,3);
%Path='c:\';
%Path='C:\Users\ALysko\Downloads\splat131\splat123\splat-1.2.3-win32\';
Path='C:\Users\ALysko\Downloads\splat131\';
%PathSplat = 'C:\Users\ALysko\Downloads\splat131\splat123\splat-1.2.3-win32\';
PathSplat = 'C:\Users\ALysko\Downloads\splat131\';
PathSDF = 'F:\Africa_hgt_terrain';
% generate TX file
if isfield(TX,'FileName'), fnTX = TX.FileName; else fnTX = fullfile(Path,'tmp_splat_tx.qth'); end;
if isfield(TX,'SiteName'), SiteName=TX.SiteName; else SiteName='TX'; end;
fid=fopen(fnTX,'wt');
if fid<=0, error('could not create TX file'); end;
fprintf(fid,'%s\n%g\n%g\n%g meters\n',SiteName,TX.Lat,TX.Lon,TX.AntHeightAGL);
fclose(fid);
% generate RX file
if isfield(RX,'FileName'), fnRX = RX.FileName; else fnRX = fullfile(Path,'tmp_splat_rx.qth'); end;
if isfield(RX,'SiteName'), SiteName=RX.SiteName; else SiteName='RX'; end;
fid=fopen(fnRX,'wt');
if fid<=0, error('could not create RX file'); end;
fprintf(fid,'%s\n%g\n%g\n%g meters\n',SiteName,RX.Lat,RX.Lon,RX.AntHeightAGL);
fclose(fid);
%
ExtraArg = '-metric';
if nargin>2 & ischar(a3),
    ExtraArg=[ExtraArg ' ' a3];
elseif nargin>2 & isstruct(a3),
    LRP = a3;
elseif nargin<3,
    % - sample -
    % 15.000 ; Earth Dielectric Constant (Relative permittivity)
    % 0.005 ; Earth Conductivity (Siemens per meter)
    % 301.000 ; Atmospheric Bending Constant (N-units)
    % 647.000 ; Frequency in MHz (20 MHz to 20 GHz)
    % 5 ;Radio Climate (5 = Continental Temperate)
    % 0 ;Polarization (0 = Horizontal, 1 = Vertical)
    % 0.50 ; Fraction of situations (50% of locations)
    % 0.90 ; Fraction of time (90% of the time)
    % 46000.0 ; ERP in Watts (optional)
    
    LRP.EarthDielectricConstant=15; % (Relative permittivity)
    LRP.EarthConductivity=0.005;    % (Siemens per meter)
    LRP.AtmosphericBendingConstant=301; % (N-units)
    LRP.FrequencyInMHz=647;         % (20 MHz to 20 GHz)
    LRP.RadioClimate=5;             % (5 = Continental Temperate)
    LRP.Polarization=0;             % (0 = Horizontal, 1 = Vertical)
    LRP.FractionOfSituations=0.5;   % (50% of locations)
    LRP.FractionOfTime=0.9;         % (90% of the time)
    LRP.ERPinWatts=0;               % (optional)
else
    error('unsupported combination of parameters'); 
end;
% generate LRP file
%if isfield(LRP,'FileName'), fnLRP = LRP.FileName; else fnLRP = fullfile(Path,'tmp_splat_lrp.lrp'); end;
[pathstr,name,ext] = fileparts( fnTX );
fnLRP1 = fullfile(pathstr,[name '.lrp']);
fnLRP = fullfile(pathstr,['splat' '.lrp']);
fid=fopen(fnLRP,'wt');
if fid<=0, error('could not create LRP file'); end;
fprintf(fid,'%.3f ; Earth Dielectric Constant (Relative permittivity)\n',LRP.EarthDielectricConstant);
fprintf(fid,'%.3f ; Earth Conductivity (Siemens per meter)\n',LRP.EarthConductivity);
fprintf(fid,'%.3f ; Atmospheric Bending Constant (N-units)\n',LRP.AtmosphericBendingConstant);
fprintf(fid,'%.3f ; Frequency in MHz (20 MHz to 20 GHz)\n',LRP.FrequencyInMHz);
fprintf(fid,'%d ;Radio Climate (5 = Continental Temperate)\n',LRP.RadioClimate);
fprintf(fid,'%d ;Polarization (0 = Horizontal, 1 = Vertical)\n',LRP.Polarization);
fprintf(fid,'%.2f ; Fraction of situations (50%% of locations)\n',LRP.FractionOfSituations);
fprintf(fid,'%.2f ; Fraction of time (90%% of the time)\n',LRP.FractionOfTime);
fprintf(fid,'%.1f ; ERP in Watts (optional)\n',LRP.ERPinWatts);
fclose(fid);
copyfile( fnLRP, fnLRP1 );
% run SPLAT!
% [status,cmdout] = dos(command,'-echo')
%splat -t sample_data\wnju-dt2.qth -r sample_data\wnju-dt.qth -metric
%ToSplat = fullfile(PathSplat,'splat.exe');
ToSplat = fullfile(PathSplat,'Splat-1-3-1-SD-mx49.exe');
% http://stackoverflow.com/questions/8055371/how-to-run-two-commands-in-one-line-in-windows-cmd
CmdLine = ['cd "' PathSplat '" && ' ...
    ToSplat ' -t "' fnTX '" -r "' fnRX '" -d ' PathSDF ' ' ExtraArg],
[status,cmdout] = dos( CmdLine ),
cmdout,
if status~=0, error('could not execute Splat'); end;
% get results
fn = [TX.SiteName '-to-' RX.SiteName '.txt'];
fid = fopen( fullfile(Path,fn), 'rt');
if fid<=0, error('could not open results file'); end;
s = fgetl(fid);
while ~feof(fid) & isempty(regexp(s,'Summary for the link between'))
    s = fgetl(fid);
end;
s = fgetl(fid);
d1 = fscanf(fid,'Free space path loss: %g dB');
s = fgetl(fid);
d2 = fscanf(fid,'Longley-Rice path loss: %g dB');
%
s1 = fgetl(fid);
while ~feof(fid) & isempty(regexp(s1,'must be raised to at least','once'))
    s2=s1; s1 = fgetl(fid);
end;
%s1,s2,
i = regexp( s1, 'raised to at least');
H1 = sscanf(s1(i:end),'raised to at least %g meters AGL');
s = fgetl(fid);
s = fgetl(fid);
s = fgetl(fid);
H1_60=nan;
% i = regexp( s, 'raised to at least');
% H1_60 = sscanf(s(i:end),'raised to at least %g meters AGL');
fclose(fid);
% NEED TO SUPPORT:
% No obstructions to LOS path due to terrain were detected by SPLAT!
% The first Fresnel zone is clear.
% 60% of the first Fresnel zone is clear.
FSL=d1;
PL=d2;
return
%% TEST
% Need to do
% 1) test the height
% 2) test over free space
% 3) test accross a hill
% compare w http://www.wirelesscommunication.nl/reference/chaptr03/pel/loss.htm
% CSIR upper: -25.742304, 28.278787
% CSIR Building Knowlenge Commons: latitude -25.7452573�, longitude 28.2781779�
% CSIR Telkom tower: -25.755417, 28.282828
% CSIR Building #44: latitude -25.757038�, longitude 28.280121�
% CSIR lower (Meraka): -25.755902, 28.278509
fclose all
clear TX RX LRP;
TX = struct('SiteName','SiteTx', 'Lat',-25.757038, 'Lon', -28.280121,'AntHeightAGL',10);    % DPSS
%RX = struct('SiteName','SiteRx', 'Lat',-25.7452573, 'Lon',-28.2781779,'AntHeightAGL',10); % KC
RX = struct('SiteName','SiteRx', 'Lat',-25.755417, 'Lon', -28.278509,'AntHeightAGL',10); % Telkom
%a3 = '-f 600';
a3 = '';
%[FSL,PL,H1,H1_60] = SPLAT_PL(TX,RX,a3),
[FSL,PL,H1,H1_60] = SPLAT_PL(TX,RX),
