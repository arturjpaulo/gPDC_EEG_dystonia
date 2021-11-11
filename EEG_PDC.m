
clear all
close all
clc
files= dir('C:\Users\artur\Documents\MATLAB\Mestrado\Selected_data');

data={};
j=1;
for i=3: length(files)
    
    if isequal(files(i).name(10),'w')
        
        
        data{j}= load([files(i).folder '\' files(i).name]);
       j=j+1;
    end
end
writting= files(i).name;

data=data';
%%
for i=1:length(data)
    
datafilt{i,1} = eegfilt (data{i,1}.resting1,250, 15, 25);

end


for i=1:length(datafilt)
    

c= datafilt{i};

Fz=c(2,:);
C3=c(8,:);
Oz=c(17,:);
Cz=c(24,:);
C4=c(25,:);


COI=vertcat(Fz,C3,Oz,Cz,C4,Cz);
COI= COI ( ~any(isnan(COI)| isinf (COI),2),:);
u=COI';
chLabels={'Fz','C3','Oz','Cz','C4','Cz'};
fs = 250;
nFreqs = 64; 
metric = 'diag';
maxIP = 30;
flgDetrend=1; 

[nSegLength,nChannels]=size(u);
if nSegLength < nChannels, error('The data might be transposed.'); end;

%===========================#
%  Channel identification   /
%===========================#

flgLabels = ~isempty(chLabels);
if flgLabels,
  if nChannels ~= max(size(chLabels))
    error('Numbers of labels and channels do not match.')
  end;
end;

%<***> Usually it's recommended to detrend the time series.

flgStandardize=0; %<***> For gPDC and iPDC, normalization has no effect.
if flgStandardize,
  disp('Be aware that the data normalization does not affect the generalized')
  disp('   PDC estimates nor its statistics results, so that data normalization')
  disp('   is not necessary.')
end;


alg = 1;
criterion = 1; 
alpha = 0.01;       
gct_signif = alpha; 
igct_signif = alpha; 
VARadequacy_signif = 0.05; 


flgColor = [0]; % Plotting option for automatic scaling for small PDC
                  % values.
                  % if flgColor = 0, y-axis scale = [0 1]
                  % elseif flgColor = 1, the pdc_xplot routine rescales 
                  % the y-axis automatically according to the following 
                  % rules:
                  %   if .01<=PDC(f) < .1 background-color = light-blue,
                  %                          so that y-axis scale = [0 .1]
                  %   elseif PDC(f) < .01 background-color = light-purple
                  %                          and y-axis = [0 .01];
                  % for flgColor=[0 1], both lay-outs are plotted.

%           [1 2 3 4 5 6 7]
flgPrinting=[1 1 1 0 0 0 2];
%            | | | | | | 7 Power Spectra (0: w/o SS; 1: Linear; 2: Log-scale)
%            | | | | | 6 Coherence
%            | | | | 5 Plot lower confidence limit
%            | | | 4 Plot upper confidence limit
%            | | 3 Significant PDC in red line
%            | 2 Patnaik threshold level in black dashed-line
%            1 PDC in green line

axis_scale = [0 0.50 -0.02 1.05];
w = fs*(0:(nFreqs-1))/2/nFreqs;
w_max = fs/2; %<***> Usually half of sampling frequency = Nyquist frequency


%==========================================================================
%==========================================================================
%        ATTENTION: BELOW THIS LINE PROBABLY YOU MIGHT NOT WANT TO EDIT,
%            UNLESS YOU WANT TO CUSTOMIZE YOUR ANALYSIS ROUTINE.
%==========================================================================

%==========================================================================
%                    Detrend and normalization options
%==========================================================================






[nChannels,nSegLength]=size(u);
if nChannels > nSegLength, u=u.'; 
   [nChannels,nSegLength]=size(u);
end;



if flgDetrend,
   for i=1:nChannels, u(i,:)=detrend(u(i,:)); end;
   disp('Time series were detrended.');
end;

flgLabels = ~isempty(chLabels);
if flgLabels,
  if nChannels ~= max(size(chLabels))
    error('Numbers of labels and channels do not match.')
  end;
end;



%==========================================================================
% Additional info for title (optional)

strTitle1 = ['PDC(' '{\alpha = ' int2str(100*alpha) '%}' ') '];
switch metric
   case 'euc'
      %NOP
   case 'diag'
      strTitle1 = ['g' strTitle1];
   case 'info'
      strTitle1 = ['i' strTitle1];
   otherwise
      error('Unknown metric.')
end;
% or set strTitle1 = [];

%==================
switch alg
  case 1
    disp('VAR estimation using Nutall-Strand algorithm.')
  case 2
    disp('VAR estimation using least-squares estimator.')
  case 3
    disp('VAR estimation using Vieira-Morf algorithm.')
  case 4
    disp('VAR estimation using QR-Arfit algorithm.')
end;

%============================#
%MAR order selection criteria/
%============================#
switch criterion
   case 1
      disp('Model order selection criteria: AIC.')
   case 2
      disp('Model order selection criteria: Hanna-Quinn.')
   case 3
      disp('Model order selection criteria: Schwartz (BIC).')
   case 4
      disp('Model order selection criteria: FPE.')
   case 5
      disp('Model order selection criteria: fixed order in maxIP.')
   otherwise
      error('Model order selection criteria: NOT IMPLEMENTED YET.')
end;

%==========================================================================
%                            VAR model estimation
%==========================================================================
[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);


disp(['Number of channels = ' int2str(nChannels) ' with ' ...
  int2str(nSegLength) ' data points; MAR model order = ' int2str(IP) '.']);

%==========================================================================
%    Testing for adequacy of MAR model fitting through Portmanteau test
%==========================================================================
   h = 20; % testing lag
   aValueVAR = 1 - VARadequacy_signif;
   flgPrintResults = 1;
[Pass,Portmanteau,st,ths]=mvarresidue(ef,nSegLength,IP,aValueVAR,h,...
                                                          flgPrintResults);

%==========================================================================
%         Granger causality test (GCT) and instantaneous GCT
%==========================================================================
   flgPrintResults = 1;
[Tr_gct, pValue_gct, Tr_igct, pValue_igct] = gct_alg(u,A,pf,gct_signif, ...
                                                         flgPrintResults);

%==========================================================================
%            PDC, threshold and confidence interval calculation.
%==========================================================================

% if alpha == 0, no asymptotic statistics is performed. ASYMP_PDC returns
% only the PDC. This option is much faster!!
 c=asymp_pdc(u,A,pf,nFreqs,metric,alpha);

% Power spectra and coherence calculation
c.SS = ss_alg(A, pf, nFreqs);
c.coh = coh_alg(c.SS);

% Statistically significant PDC on frequency scale
if alpha ~= 0,
   pdc_temp = ((abs(c.pdc)-c.th) > 0).*c.pdc + ((abs(c.pdc)-c.th) <= 0)*(-1);
   pdc_temp(ind2sub(size(pdc_temp),find(pdc_temp == -1))) = NaN;
   c.pdc_th = pdc_temp;
end;

%Adding further analysis details in the figure title.
%strTitle3 = ['[N=' int2str(nSegLength) '; IP=' int2str(c.p) ']'];
% or

strTitle3 = ['[N=' int2str(nSegLength) 'pts; IP=' int2str(c.p) '; ' ...
   datestr(now) ']'];

% or leave emptied: strTitle3=[];

%==========================================================================
%              Matrix Layout Plotting of the Analysis Results
%==========================================================================

w_max = fs/2;
%strTitle = [strTitle1 strTitle2 strTitle3];
strTitle = [strTitle1 strTitle3]
strWindowName = 'pdc Analysis Template Example';

% The following "for loop" through flgColor values, 0 and 1, and yields a
% pair of plots, one without and other with color scale rearrangement option.
% Value range of PDC and Coherence is from [0 1], but sometimes the maximum 
% peak value is small (<0.1), or even smaller, (<.01), so in these cases it
% might be interesting to have a plot with finer smaller y-axis scale. The
% white-background plot indicates full-scale [0 1] y-axis, while
% light-blue-background stands for intermediate [0 .1] scaling and
% light-purple-background shows very fine detail of small, usually not
% significant PDCs. Try flgColor = 0 or 1, or both [0 1].

for kflgColor = flgColor,
   h=figure;
   set(h,'NumberTitle','off','MenuBar','none', ...
      'Name', strWindowName )

   [hxlabel hylabel] = pdc_xplot(c,...
      flgPrinting,fs,w_max,chLabels,kflgColor);
   
% The title suplabel command should (not sure) follow the pdc_xplot routine
% In MacOS X, for flgPrinting(7) = 4 or 5, the main diagonal plotting
% gets misaligned if suplabel with 't' option is used more than once.

   [ax,hT]=suplabel( strTitle, 't' );
   set(hT,'FontSize',8) 
end;


%======================= pdc_xplot ========================================
%Plot legend:  Blue lines on the main diagonal = Power spectra;
%              Black dashed lines are Patnaik threshold for pdcn;
%              Green lines = non significant pdcn;
%              Red lines = significant pdcn;
%              Light-gray lines = coherence function.
%

%raw3(i).description

clear IP pf A pb B ef eb vaic Vaicv
 end

