
load dataset

wave2=wave2fullclusterestimates;

inds{1}=find((wave2.Type=='Aged care') & (wave2.cluster~=82));

inds{2}=find(wave2.Type=='Healthcare');

inds{3}=find(wave2.Type=='Housing');

inds{4}=find(wave2.Type=='Packing plant/meat processing');

inds{5}=find(wave2.Type=='School');

%inds{6}=find(wave2.Type=='NA');

numtyp=length(inds);
Mu=zeros(numtyp,1);
Mu_se=zeros(numtyp,1);

for typ=1:numtyp
    
    mu=wave2.mu(inds{typ});
    mu_se=wave2.mu_se(inds{typ});
    
    Mu(typ)=sum( mu./mu_se.^2)/sum(1./mu_se.^2);
    Mu_se(typ)= sqrt(1/sum(1./mu_se.^2));
    
end


figure(1)
clf

upper=Mu+1.96*Mu_se;
lower=Mu-1.96*Mu_se;

small=0.1;

map=colormap;

for k=1:numtyp
    
    plot([k k],[lower(k) upper(k)],'.-','Color',map(k*50,:),'markersize',20,'linewidth',2)
    
    hold on
    
   plot([k-small k+small ],[Mu(k) Mu(k)],'Color',map(k*50,:),'linewidth',2)
    
end

axis([ 0 numtyp+1 min(lower)-3 max(upper)+1])

title('Estimates of mean serial interval by cluster site type','fontname','sansserif','fontsize',15)

grid on

hh=0;
text(1,hh,'Aged Care','Rotation',90,'fontname','sansserif','FontSize',15);
text(2,hh,'Healthcare','Rotation',90,'fontname','sansserif','FontSize',15);
text(3,hh,'Housing','Rotation',90,'fontname','sansserif','FontSize',15);
text(4,hh,'PP/MP','Rotation',90,'fontname','sansserif','FontSize',15);
text(5,hh,'School','Rotation',90,'fontname','sansserif','FontSize',15);
%text(6,hh,'NA','Rotation',90);


print -depsc mean_by_setting


% do we think there are different means in different settings
overall_Mu=sum(Mu./Mu_se.^2)/sum(1./Mu_se.^2);
Stat=sum((overall_Mu-Mu).^2./Mu_se.^2);
chi2cdf(Stat,numtyp-1,'Upper')

% here we delete the second one
Removal=[ 1 2 3 5];
MuR=Mu(Removal,1);
Mu_seR=Mu_se(Removal,1);
overall_MuR=sum(MuR./Mu_seR.^2)/sum(1./Mu_seR.^2);
Stat=sum((overall_MuR-MuR).^2./Mu_seR.^2);
chi2cdf(Stat,length(Removal)-1,'Upper')

