clear
path_data='/Users/lyly/scripts/rbans_corrections_bootstrap/RBANS_bootstrap_adjustment_analysis';
cd (path_data)
[num,txt]=xlsread('RBANS_20170427_bootsrap_adjusted.xls');
id=num(2:end,1);
time=num(2:end,2);
NAP=num(2:end,4);
french=num(2:end,9);
pass=num(2:end,6);

version_txt=txt(2:end,7);
version=zeros(size(id));
version(strcmpi(version_txt,'B')==1)=1;
version(strcmpi(version_txt,'C')==1)=2;
version(strcmpi(version_txt,'D')==1)=3;
rbans=num(2:end,10:15);% rbans
%tbl=table(id,time,pass,NAP,version,french);
% Data filtering
c1=french==1;
c2=pass==1;
c3=id~=211391;
c4=time==0 | time==3;
c5=NAP==1;
keep=c1&c2&c3&c4&c5;

temp=[id time version french];
X=temp(keep,:);
cog=rbans(keep,:);
% testing for id having only one time point
id_test=X(:,1);
for i=1:length(id_test);ind=find(id_test==id_test(i));n(i,1)=length(ind);end
% keeping id with both BL and 3 month
X03=X(n==2,1:3);l=length(unique(X03(:,1))) % number of subject final
% to be used for analysis 
cog03adj=cog(n==2,:);
id03=X03(:,1);
time03=X03(:,2);
version03=X03(:,3);

%% loading unadjusted data
[num,txt]=xlsread('RBANS_index_scores_unadjusted.xls');
idraw=num(:,1);timeraw=num(:,2);rbans_temp=num(:,4:9);
% merging with id we kept before
cog03raw=NaN(length(id03),6);
for i =1:length(id03)
    s=idraw==id03(i);
    v=timeraw==time03(i);
    ind=s&v;
    cog03raw(i,:)=rbans_temp(ind,:);
end

%% Getting wide format
t0=time03==0;
t3=time03==3;
v3=version03(t3);% 1=B, 2=C, 3=D
id0=id03(t0);id3=id03(t3);test=prod(double(id0==id3))% will be 1 if id are in the same order
cog0adj=cog03adj(t0,:);
cog3adj=cog03adj(t3,:);
cog0raw=cog03raw(t0,:);
cog3raw=cog03raw(t3,:);

% You can calculate your formula from here :-) ...


%% BOOTSTRAP : testing for one cog 
% M=[X03 cog03adj(:,1)];
% tbl=array2table(M,'VariableNames',{'id','time','version','cog'});
% tbl.version=categorical(tbl.version);
% lme=fitlme(tbl,'cog~version+(1|id)','FitMethod','REML')
% coef=double(lme.Coefficients(2:end,2))';

% keeping id with both BL and 3 month original code
M=[X03 cog03adj];
M3=M(:,[1:4]);
tbl=array2table(M3,'VariableNames',{'id','time','version','cog'});
tbl.version=categorical(tbl.version);
lme1=fitlme(tbl,'cog~version+(1|id)','FitMethod','REML');
coef1=double(lme1.Coefficients(2:end,2))';

M3=M(:,[1:3 5]);
tbl=array2table(M3,'VariableNames',{'id','time','version','cog'});
tbl.version=categorical(tbl.version);
lme2=fitlme(tbl,'cog~version+(1|id)','FitMethod','REML');
coef2=double(lme2.Coefficients(2:end,2))';

M3=M(:,[1:3 6]);
tbl=array2table(M3,'VariableNames',{'id','time','version','cog'});
tbl.version=categorical(tbl.version);
lme3=fitlme(tbl,'cog~version+(1|id)','FitMethod','REML');
coef3=double(lme3.Coefficients(2:end,2))';

M3=M(:,[1:3 7]);
tbl=array2table(M3,'VariableNames',{'id','time','version','cog'});
tbl.version=categorical(tbl.version);
lme4=fitlme(tbl,'cog~version+(1|id)','FitMethod','REML');
coef4=double(lme4.Coefficients(2:end,2))';

M3=M(:,[1:3 8]);
tbl=array2table(M3,'VariableNames',{'id','time','version','cog'});
tbl.version=categorical(tbl.version);
lme5=fitlme(tbl,'cog~version+(1|id)','FitMethod','REML');
coef5=double(lme5.Coefficients(2:end,2))';

M3=M(:,[1:3 9]);
tbl=array2table(M3,'VariableNames',{'id','time','version','cog'});
tbl.version=categorical(tbl.version);
lme6=fitlme(tbl,'cog~version+(1|id)','FitMethod','REML');
coef6=double(lme6.Coefficients(2:end,2))';

coef_total=[coef1' coef2' coef3' coef4' coef5' coef6']

%% create a vector of even index

indices = 2:2:266; % to change only the even values (i.e 3m visit in this case )

bstr_betas = nan(5000,3,6);

for k=1:size(cog03adj,2) 
    
    M=[X03 cog03adj(:,k)]; % merging with cognition k
    
    for i = 1:5000 % N botstrp
        
        vector = datasample(indices,133); % sample 133 random indices
        
        new_table = [];
        for j = 1 : 133
            
            new_table = [new_table;M(vector(j)-1,:);M(vector(j),:)]; % j-1 is the untouch time 0
            
        end
        
        tbl2 = array2table(new_table,'VariableNames',{'id','time','version','cog'});
        tbl2.version=categorical(tbl2.version);
        lme=fitlme(tbl2,'cog~version+(1|id)','FitMethod','REML');
        bstr_betas(i,:,k)=double(lme.Coefficients(2:end,2))';
        [i k]
    end
end


save bstr5000_cog_version_betas.mat bstr_betas 
    
% %% looking at bootsrap results
% clear
path_data='/Users/lyly/scripts/rbans_corrections_bootstrap/RBANS_bootstrap_adjustment_analysis';
 figure
for i=1:6
    x1=bstr_betas(:,1,i);
    subplot(3,6,i);hist(x1);ylabel('version B');xlabel(num2str(i));xlim([-5 5])
    x2=bstr_betas(:,2,i);
    subplot(3,6,i+6);hist(x2);ylabel('version C');xlabel(num2str(i));xlim([-5 5])
    x3=bstr_betas(:,3,i);
    subplot(3,6,i+12);hist(x3);ylabel('version D');xlabel(num2str(i));xlim([-5 5])
end
% %% Get confidence interval
MEAN=NaN(3,6);
SEM=NaN(3,6);
CI=NaN(3,6,2);
for i =1:3
    for j=1:6
        x=bstr_betas(:,i,j);
        MEAN(i,j)=mean(x);
        SEM(i,j) = std(x)/sqrt(length(x));          % Standard Error
        ts = tinv([0.025  0.975],length(x)-1);      % T-Score
        ci= MEAN(i,j) + ts*SEM(i,j);
        CI(i,j,1)=ci(1);                            % Confidence Intervals 1
        CI(i,j,2)=ci(2);                            % Confidence Intervals 2
    end
end

Matrix=[id0 v3 cog0adj cog3adj]

stda=std(Matrix(:,3));
stdb=std(Matrix(Matrix(:,2)==1,[9 10 11 12 13 14]));
stdc=std(Matrix(Matrix(:,2)==2,[9 10 11 12 13 14]));
stdd=std(Matrix(Matrix(:,2)==3,[9 10 11 12 13 14]));

ratiob=stdb/stda
ratioc=stdc/stda
ratiod=stdd/stda

meana=mean(Matrix(:,3:8));
meanb=mean(Matrix(Matrix(:,2)==1,[9 10 11 12 13 14]));
meanc=mean(Matrix(Matrix(:,2)==2,[9 10 11 12 13 14]));
meand=mean(Matrix(Matrix(:,2)==3,[9 10 11 12 13 14]));

for i=1:6
    differenceb(:,i)=(Matrix(Matrix(:,2)==1,8+i))-(meanb(:,i));
    differencec(:,i)=(Matrix(Matrix(:,2)==2,8+i))-(meanc(:,i));
    differenced(:,i)=(Matrix(Matrix(:,2)==3,8+i))-(meand(:,i));
end

for i=1:6
    adjustedb(:,i)=(ratiob(:,i)*differenceb(:,i))+meana(:,i);
    adjustedc(:,i)=(ratioc(:,i)*differencec(:,i))+meana(:,i);
    adjustedd(:,i)=(ratiod(:,i)*differenced(:,i))+meana(:,i);
end

for i=1:6
    lmeadjustedb(:,i)=(Matrix(Matrix(:,2)==1,8+i));
    lmeadjustedc(:,i)=(Matrix(Matrix(:,2)==2,8+i));
    lmeadjustedd(:,i)=(Matrix(Matrix(:,2)==3,8+i));
end

