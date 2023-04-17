function [data_2D,output]=PZS(folder,N,freq_r,norder,nExpDec,datatype)
%% Pseudo-Zernike similarity (PZS) - an alternative FFCF extraction method 

% PZS can be used to evaluate FFCF from 2DIR spectra of a
% single or two-component system (overlapping spectral bands - resolved
% or unresolved).


%% Input parameters

% folder (required)  -> foldername where 2DIR matlab figures are saved.
                      % Matlab Figure Format 
                      % FigureName -'Tw_()' where () is filled with waiting time in femtoseconds
                      % FigureAxes - wtau (y-direction) and wt (x-direction)
% N      (required)  -> Different analysis ('single','overlapping') 
% freq_r (optional)  -> Frequency window parameters (cm-1)
                    % scalar input - [radius] (single component anaysis)
                    % vector input - [wtL wtU wtauL wtauU] (multi-component analysis) 
                        %wtL - lower bound for wt frequency 
                        %wtU - upper bound for wt frequency
                        %wtauL - lower bound for wtau frequency
                        %wtauU - upper bound for wtau frequency
                    % 'Default' 
                    % 'single' - full maximum width at 90% of 2DIR trace at t0. 
                    % 'overlapping' - full width maximum at 75% of 2DIR trace at t0.
% norder (optional)  -> order of pseudo-Zernike moment ( norder larger...
                        ... than 4 is recommended)
                    % 'Default' - 15
% nExpDec (optional) -> single exponential decay fit (1)...
                        ... area under the curve for multi-exponential decay fit (0)
                    % 'Default' - 1
% datatype (optional) -> 'psc' - peak shift correction
                       % 'Default' - '' (no peak shift correction)
%% Output 

    % data_2D -> spectral data extracted from 2DIR matlab figures
    % output -> struct file ('spectral_region','analysis','results','figures')

%% Extracting data (frequencies, intensity and waiting time) from 2DIR figures in a specific folder
cf=pwd;
cd(folder)
fig_list=dir('*Tw*.fig'); %% Extract the 2DIR matlab figure names (includes waiting times)

%Pre-allocation
n_list=length(fig_list);
data_cell=cell(n_list,4);
Tw=nan(n_list,1);

%% Running for loop to extract data for all matlab figures
for j=1:n_list
    figname=fig_list(j).name;
    openfig(figname,'invisible');
    fig_data=get(gca,'Children');
    all_data=get(fig_data,{'ZData','YData','XData'});
    tw = regexp(figname,'\d*','Match');
    data_cell(j,:)=[all_data(size(all_data,1),:),str2double(tw)];
    Tw(j,1)=str2double(tw);
    close;
end
data_s=sortrows(data_cell,4); % Sorting the extracted data with chronological waiting times

%% Converting spectral data from cell type to matrix 
Sc=data_s(:,1);
[n_wtau,n_wt]=size(Sc{1,1});
S_mat=nan(n_wtau,n_wt,n_list); 

for i=1:n_list
    S_mat(:,:,i)=Sc{i,1};
end

%% Save input 2DIR data
data_2D.spectra=S_mat;
data_2D.Tw_fs=sort(Tw);
data_2D.wtau=data_s{1,2}';
data_2D.wt=data_s{1,3};

cd(cf);

%% Spectral Data Analysis
    if nargin < 3
        output=spectral_region(data_2D,N);
    elseif nargin < 4
        output=spectral_region(data_2D,N,freq_r);
    elseif nargin < 5
        output=spectral_region(data_2D,N,freq_r,norder);
    elseif nargin < 6
       output=spectral_region(data_2D,N,freq_r,norder,nExpDec);
    else
       output=spectral_region(data_2D,N,freq_r,norder,nExpDec,datatype); 
    end

end




function [output]=spectral_region(data_2D,N,freq_r,norder,nExpDec,datatype)
%% Input 
    % data_2D -> data extracted from 2DIR matlab figures
    % N ->  number of transitions (spectral bands)
    % freq_r -> delta frequency (radius) from center frequencies in either wtau\wt direction
    % norder -> order of pseudo-Zernike moment
    % nExpDecay -> Single exponential decay (1) or Area under the curve for multi-exponential decay(0)

%% Assigning 2DIR data 
wt=data_2D.wt; % pump frequencies (cm-1)
wtau=data_2D.wtau; % probe frequencies (cm-1)
Tw=data_2D.Tw_fs; % waiting times (fs)
S=data_2D.spectra; %  Raw 2DIR spectra for all waiting times
Smax=S./max(S,[],[1 2]); % Max normalized 2DIR spectrum at each waiting time
n_Tw=length(Tw); % number of waiting times

if nargin < 4 || isempty(norder)
    norder=15;
end

if nargin < 5 || isempty(nExpDec)
    nExpDec=1;
end
%Figure Size
figL=500;
figW=400;
%Figure Position
figV=300;
figH=300;

%% Selecting spectral region of interest based on N
switch(N)

   
%% Analysis of single component system (non-overlapping spectral bands)       
    case 'single'
       
        if nargin < 3 || isempty(freq_r)
            
            [ymax,xmax]=find(Smax(:,:,1)==1);
            if abs(ymax-xmax)>2
                wt_shift=(xmax-ymax);
                S0_diag=diag(Smax(:,:,1),wt_shift);
            else
                S0_diag=diag(Smax(:,:,1));     
            end
            [~,max_loc]=max(S0_diag);
            peak_ratio=0.1;
            [r_max]=find(S0_diag>(min(S0_diag)+peak_ratio*range(S0_diag)));
            rcheck=ischange(r_max,"linear");
            ri=find(rcheck==1);
            if ~isempty(ri)
                if max(ri)<max_loc
                    r_max=r_max(max(ri)+1:end);
                elseif min(ri)>max_loc
                    ri=find(rcheck,1,'first');
                    r_max=r_max(1:min(ri)-1);
                elseif max(ri)>max_loc && min(ri)<max_loc
                    ru=ri(ri>max_loc);
                    rl=ri(ri<max_loc);
                    r_max=r_max(max(rl)+1:min(ru)-1);
                end
            end
            hf=ceil((max(r_max)-max_loc)/(2*cosd(45)));
            lf=ceil((max_loc-min(r_max))/(2*cosd(45)));
                if hf>=lf
                    freq_r=lf;
                elseif lf>hf
                    freq_r=hf;
                end
        end

        
        Pk_freq=nan(n_Tw,2);
        di=freq_r*4+1;
        S_p=nan(di,di,n_Tw);
        wtau_c=nan(di,n_Tw);
        wt_c=nan(n_Tw,di);
        cls=nan(n_Tw,1);    


        for ti=1:n_Tw
            Smn=Smax(:,:,ti);
            [n_r,n_c]=find(Smn==1);
            Pk_freq(ti,:)=[wtau(n_r),wt(n_c)];

                s_m=max(Smn,[],2);
                s_m(s_m<0.5)=0;
                s_l_m=s_m(islocalmax(s_m));
                if length(s_l_m)>1
                    [n_r1,n_c1]=find(Smn==s_l_m(1));
                    [n_r2,n_c2]=find(Smn==s_l_m(2));
                    if s_l_m(1)>s_l_m(2)
                        if abs(Pk_freq(1,1)-wtau(n_r1)) < (abs(wtau(n_r1)-wtau(n_r2))/2)
                            n_r=n_r1;
                            n_c=n_c1;
                        else
                            n_r=n_r2;
                            n_c=n_c2;
                        end
                    elseif s_l_m(1)<s_l_m(2)
                        if abs(Pk_freq(1,1)-wtau(n_r2)) < (abs(wtau(n_r1)-wtau(n_r2))/2)
                            n_r=n_r2;
                            n_c=n_c2;
                        else
                            n_r=n_r1;
                            n_c=n_c1;
                        end
                    end
                    Pk_freq(ti,:)=[wtau(n_r),wt(n_c)];
                end
            
            nrow_c=n_r-freq_r*2:n_r+freq_r*2;
            ncol_c=n_c-freq_r*2:n_c+freq_r*2;
            S_p(:,:,ti)=Smn(nrow_c,ncol_c);
            wtau_c(:,ti)=wtau(nrow_c);
            wt_c(ti,:)=wtau(ncol_c);



        mp=ceil(di/2);
        jstart=6; % lower limit for radius of frequency window
        nsize=freq_r-jstart+1;
        S_o=cell(nsize,1);
        Pv=cell(nsize,1);
        ffcf=nan(nsize,n_Tw);
        tau=nan(nsize,1);
        
        for j=jstart:freq_r
            ji=j-jstart+1;
            r=mp-2*j:mp+2*j;
            sp=S_p(r,r,:);
            S_o{ji,1}=sp;
            rl=length(r);
            pc=PZP_calc(rl,norder);
            
           sm=PZM(sp,pc);
           ffcf(ji,:)=similarity(sm);
           tau(ji,1)=TimeDecay(ffcf(ji,:)',Tw,nExpDec);
            Pv{ji,1}=sm;
            
        end
        output.spectral_region.Pv=Pv;
        output.spectral_region.Sp=S_o;
        output.spectral_region.center_freq_loc=Pk_freq;
        output.spectral_region.wt=wt_c;
        output.spectral_region.wtau=wtau_c;

        output.analysis.FFCF=ffcf;
        output.analysis.tau_in_ps=tau;

%% Frequency Window Selection
        tau_smooth=smooth(tau,1);
       [mind]=islocalmax(tau_smooth);

       if length(tau_smooth(mind))<1 
            [mind]=islocalmin(tau_smooth);
                if ~isempty(tau_smooth(mind))
                    mind=find(mind==1);
                    mind=mind(1);
                else
                   mind=find(tau_smooth==median(tau_smooth));
                   if ~isempty(mind)
                   mind=mind(1);
                   else
                       mind=1;
                   end
                end
       elseif length(tau_smooth(mind))==2
           mind=find(mind==1);
           mind=mind(1);
       elseif length(tau_smooth(mind))==1
           mind=find(mind==1);
       else
           mind=find(mind==1);
           mind=mind(1);
       end
        rfw=jstart:freq_r;
        output.results.FFCF=ffcf(mind,:);
        output.results.tau_in_ps=tau(mind);
        output.results.frequency_radius=rfw(mind);

            % CLS
            wtau_w=rfw(mind); % frequency width in wtau direction
            w0=2*wtau_w;
            sc=Smn(n_r-w0:n_r+w0,:);
            [~,md]=max(sc,[],2);
            mdl=length(md);
            wt_w=rfw(mind); % frequency width in wt direction 
            w1=2*wt_w;
      
            wt_n=nan(mdl,1);
            for k=1:mdl
                sc_m=sc(k,md(k)-w1:md(k)+w1);
                wt_cls=wt(md(k)-w1:md(k)+w1);
                pfit2=polyfit(wt_cls,sc_m,2);
                y_new=polyval(pfit2,wt_cls);
                [~,ymax]=max(y_new);
                wt_n(k,1)=wt_cls(ymax);
            end

            wtau_cls=wtau(n_r-w0:n_r+w0);
            p1=polyfit(wtau_cls,wt_n,1);
            cls(ti,1)=p1(1);       
        end
     %   S_p=S_p./max(S_p,[],[1 2]); 
        CLS_m=cls./max(cls);

%% Create figures for frequency window optimization and FFCF results.
        close([figure(101) figure(102) figure(103) figure(104)]);
        output.fig1=[figure(101) figure(102)];
        
        %Frequency Window Optimization plot
        figure(output.fig1(1)), line([rfw(mind) rfw(mind)],[0 output.analysis.tau_in_ps(mind)],'Color','r','LineWidth',2,'LineStyle',':'); 
        hold on, scatter(rfw,output.analysis.tau_in_ps,'k','filled');hold off
        ylabel('fitted \tau_c (ps)');
        xlabel(' radius of frequency window (cm^-^1)');
        title ('Frequency window optimization');
        set(gca,'FontSize',14);
        box on;
        axis tight;
        set(gca,'Ylim',[min(output.analysis.tau_in_ps)-0.1*range(output.analysis.tau_in_ps) max(output.analysis.tau_in_ps)+0.1*range(output.analysis.tau_in_ps)]);
        set(gca,'Xlim',[jstart-1 max(rfw)+1]);
        set(gcf, 'position',[figV figH figL figW]);


        % FFCF plot
        Tw_ps=Tw./1000;
        ffcf0=output.analysis.FFCF(mind,:);
        figure(output.fig1(2)),scatter(Tw_ps,ffcf0,'filled','Marker','s','MarkerFaceColor','r','MarkerEdgeColor','r');
        ylabel('PZS (max normalized)');
        xlabel('waiting time (ps)');
        title 'FFCF approximation';
        set(gca,'FontSize',14);
        box on;
        axis tight;
        set(gca,'Ylim',[min(ffcf0)-0.1*range(ffcf0) max(ffcf0)+0.1*range(ffcf0)]);
        set(gca,'Xlim',[min(Tw_ps)-0.5 max(Tw_ps)+0.5]);
        set(gcf, 'position',[figV figH figL figW]);

S_p=Smax;

op_r=output.results.frequency_radius;
 [~,Imax]=max(S_p,[],2);
[ix,~,iz]=size(Imax);
Imax=reshape(Imax,ix,iz);
[mp,~]=find(S_p(:,:,1)==max(S_p(:,:,1),[],'all'));
ur2=round(2*op_r);
ur1=round(2*op_r);
wtau_scan=wtau(mp-ur1:mp+ur2); 
rl=4*op_r+1;
Sp=nan(rl,rl,n_Tw,ur2+ur1+1);

n=1;
for k=mp-ur1:mp+ur2
    n_r= k-op_r*2:k+op_r*2;
    for l=1:n_Tw
        n_c= Imax(k,l)-op_r*2:Imax(k,l)+op_r*2;
        if min(n_c)>=1 && max(n_c)<=ix
        Ss=Smax(n_r,n_c,l); 
        Sp(:,:,l,n)=Ss;
        end
    end
    n=n+1;
end


            pc=PZP_calc(rl,norder);
            pz=PZM(Sp,pc);
            [~,n,p]=size(pz);
            ffcf=nan(p,n);
            tau=nan(p,1);
            for k=1:p
                cp=similarity(pz(:,:,k)');
                ffcf(k,:)=cp;
                tau(k,1)=TimeDecay(cp,Tw,nExpDec);
            end
            output.CL_PZS=ffcf;
        output.fig2=[figure(103),figure(104)];


        figure(output.fig2(1)), scatter(wtau_scan,tau,'k','filled');
        ylabel('fitted \tau_c (ps)');
        xlabel('Center frequency (cm^-^1)');
        title ('PZS along centerline');
        set(gca,'FontSize',14);
   %     text(median(wtau_scan),max(tau)+0.5,['RSD = ',rsd, ' %'],'Fontsize',14);
        box on;
        axis tight;
        set(gca,'Ylim',[min(tau)-1 max(tau)+1]);
        set(gca,'Xlim',[min(wtau_scan)-5 max(wtau_scan)+5]);
        set(gcf, 'position',[figV figH figL figW]);
        output.analysis.PZS_centerline=[wtau_scan,tau];

        figure(output.fig2(2)), scatter(CLS_m,ffcf0,60,'k','filled','s');
        ylabel('PZS');
        xlabel('CLS');
        title ('CLS vs. PZS (max normalized)');
        set(gca,'FontSize',14);
        box on;
        axis tight;
        lr=lsline;
        lr.LineWidth=2;
        lf=regstats(ffcf0,CLS_m,'linear');
        r2=lf.rsquare;
        r2=sprintf('%.4f',r2);
        text(min(CLS_m)+0.25*range(CLS_m),min(ffcf0)+0.25*range(ffcf0),['R^2 = ',r2],'Fontsize',14);
        legend off;
        set(gca,'Ylim',[min(ffcf0)-0.1*range(ffcf0) max(ffcf0)+0.1*range(ffcf0)]);
        set(gca,'Xlim',[min(CLS_m)-0.1*range(CLS_m) max(CLS_m)+0.1*range(CLS_m)]);
        set(gcf, 'position',[figV figH figL figW]);

       

%% Analysis of multi-component system (Overlapping but resolved  spectral bands)       
    case 'overlapping' 
        
        if nargin < 3 || isempty(freq_r)
            datatype='';
            [ymax,xmax]=find(Smax(:,:,1)==1);
            if abs(ymax-xmax)>2
                wt_shift=(xmax-ymax);
                S0_diag=diag(Smax(:,:,1),wt_shift);
            else
                S0_diag=diag(Smax(:,:,1));
            end
            [~,max_loc]=max(S0_diag);
            peak_ratio=0.25;
            [r_max]=find(S0_diag>(min(S0_diag)+peak_ratio*range(S0_diag)));
            rcheck=ischange(r_max,"linear");
            ri=find(rcheck==1);
            if ~isempty(ri)
                if max(ri)<max_loc
                    r_max=r_max(max(ri)+1:end);
                elseif min(ri)>max_loc
                    ri=find(rcheck,1,'first');
                    r_max=r_max(1:min(ri)-1);
                elseif max_loc<max(ri) && max_loc>max(ri)
                    ru=ri(ri>max_loc);
                    rl=ri(ri<max_loc);
                    r_max=r_max(max(rl)+1:min(ru)-1);
                end
            end
            hf=ceil((max(r_max)-max_loc)/(2*cosd(45)));
            lf=ceil((max_loc-min(r_max))/(2*cosd(45)));

            s=Smax(max_loc,:,1);
            [~,smax]=max(s);
            [~,smin]=min(s);
            f2=interp1(s(smin:smax),smin:smax,0);
            anhar1=ceil(smax-f2);
            anhar2=ceil((smax-smin)/2);
                if anhar2>=anhar1
                    anharm=anhar1;
                else
                    anharm=anhar2;
                end
            if nargin<3 ||isempty(freq_r)    
            freq_r=[lf+anharm hf];
            else
                freq_r=[lf hf];
            end

            Smn=Smax(:,:,1);
            [n_r,n_c]=find(Smn==1);
            Pk_freq=[wtau(n_r),wt(n_c)];
            wtau_lm=[2 2];
            wtc=[2 2];
            nrow_c=round(n_r-lf*wtau_lm(1)):round(n_r+hf*wtau_lm(2));
            ncol_c=round(n_c-freq_r(1)*wtc(1)):round(n_c+freq_r(2)*wtc(2));
            S_p=Smax(nrow_c,ncol_c,:); 

          %%{

    if datatype == "psc"
            p_shift=nan(n_Tw,2);
            S_p=nan(length(nrow_c),length(ncol_c),n_Tw);
        for j=1:n_Tw
          %  [nr,nc]=find(Smax(:,:,j)==1);
                Smn=Smax(:,:,j);
                s_m=max(Smn,[],2);
                s_m(s_m<0.5)=0;
                s_l_m=s_m(islocalmax(s_m));
                if length(s_l_m)>1
                    [n_r1,n_c1]=find(Smn==s_l_m(1));
                    [n_r2,n_c2]=find(Smn==s_l_m(2));
                    if s_l_m(1)>s_l_m(2)
                        if abs(Pk_freq(1,1)-wtau(n_r1)) < (abs(wtau(n_r1)-wtau(n_r2))/2)
                            nr=n_r1;
                            nc=n_c1;
                        else
                            nr=n_r2;
                            nc=n_c2;
                        end
                    elseif s_l_m(1)<s_l_m(2)
                        if abs(Pk_freq(1,1)-wtau(n_r2)) < (abs(wtau(n_r1)-wtau(n_r2))/2)
                            nr=n_r2;
                            nc=n_c2;
                        else
                            nr=n_r1;
                            nc=n_c1;
                        end
                    end
                else
                    [nr,nc]=find(Smn==s_l_m(1));
                end
               
            p_shift(j,:)=[wtau(nr),wt(nc)];
            nrow_c=round(nr-lf*wtau_lm(1)):nr+round(hf*wtau_lm(2));
            ncol_c=round(nc-freq_r(1)*wtc(1)):round(nc+freq_r(2)*wtc(2));
            S_p(:,:,j)=Smax(nrow_c,ncol_c,j);
        end
    end
        %}
         wtau_c=wtau(nrow_c);
        else
            Smn=Smax(:,:,1);
            [n_r,n_c]=find(Smn==1);
            Pk_freq=[wtau(n_r),wt(n_c)];

            if length(freq_r)==2
                nrow_c=find(wt>freq_r(1) & wt<freq_r(2));
                ncol_c=find(wt>freq_r(1) & wt<freq_r(2));
            elseif length(freq_r)==4
                nrow_c=find(wt>freq_r(3) & wt<freq_r(4));
                ncol_c=find(wt>freq_r(1) & wt<freq_r(2));
            end
            
            if isempty(nrow_c) || isempty(ncol_c)
                error('Frequency window is not correctly positioned.')
            end
            S_p=Smax(nrow_c,ncol_c,:);
            wtau_c=wtau(nrow_c);
        end

       


        %% wt shift check
    
            d1=diag(Smax(:,:,1));
            d1(d1<0.5)=0;
            d1_max=islocalmax(d1);
            d1_loc=find(d1_max==1);
            res1=length(d1_loc);
            
            d2=diag(Smax(:,:,2));
            d2(d2<0.5)=0;
            d2_max=islocalmax(d2);
            d2_loc=find(d2_max==1);
            res2=length(d2_loc);
            if res2>1
                    if abs(d2_loc(1)-d1_loc(1)) < abs(d2_loc(2)-d1_loc(1))
                        d2_loc=find(d2==d2_max(1));
                    else
                        d2_loc=find(d2==d2_max(2));
                    end
            end

              if res1==1  
                    if (d2_loc-d1_loc)>4
                        %% wt shift correction
                            wt_c=nan(n_Tw,length(ncol_c));
                            for nt=1:n_Tw
                                wtau0=max(Smax(:,:,1),[],2);
                                wtau0(wtau0<0.5)=0;
                                wtau0_max=islocalmax(wtau0);
                                    wt0=find(wtau0_max==1);
                                    wt0=wt0(1);
                                    wtau_nt=max(Smax(:,:,nt),[],2);
                                    wtau_nt(wtau_nt<0.5)=0;
                                    wtau_ntmax=islocalmax(wtau_nt);
                                    wt_nt=find(wtau_ntmax==1);
                                    wt_nt=wt_nt(1);
                                    if wt_nt==wt0
                                        S_p(:,:,nt)=Smax(nrow_c,ncol_c,nt);
                                        wt_c(nt,:)=wt(ncol_c);
                                    elseif wt_nt~=wt0
                                        S_p(:,:,nt)=Smax(nrow_c,ncol_c+(wt_nt-wt0),nt);
                                        wt_c(nt,:)=wt(ncol_c+(wt_nt-wt0));
                                    end
                    
                            end
                    else

                        wt_c=wt(ncol_c);
    
                    end
              else
                  wt_c=wt(ncol_c);
              end




        di=length(ncol_c);
        dr=length(nrow_c);
        mp=ceil(di/2);
        row_w=1:dr;
        col_w=1:di;
        S_s=S_p(row_w,col_w,:);
        wtau_s=wtau_c(row_w,:);
     
        rl=length(col_w);
        pc=PZP_calc(rl,norder);
        
        [~,y,z]=size(S_s);
        sp=zeros(y,y,z);
        S_o=nan(y,y,z,dr);

        for j=3:dr-3

            r=j-2:j+2;
            sp(mp-2:mp+2,:,:)=S_s(r,:,:);
            S_o(:,:,:,j)=sp;
        end
        ind=3:dr-3;
        Sp=S_o(:,:,:,ind);
        pz=PZM(Sp,pc);
        [~,n,p]=size(pz);
        ffcf=nan(p,n);
        tau=nan(p,1);
            for k=1:p
            cp=similarity(pz(:,:,k)');
            ffcf(k,:)=cp;
            tau(k,1)=TimeDecay(cp,Tw,nExpDec);
            end
        output.spectral_region.Pv=pz;
        output.spectral_region.Sp=Sp;
        output.spectral_region.center_freq_loc=Pk_freq;
        output.spectral_region.wt=wt_c;
        output.spectral_region.wtau=wtau_c;
    
        output.analysis.FFCF=ffcf;
        output.analysis.tau_in_ps=tau;
        wtau_s=wtau_s(ind);
        output.analysis.scan_freq=wtau_s;

%% Frequency Window Selection
        
        mri=1:length(ind);
        ffcf0=ffcf(mri,:);
        tau0=tau(mri);
        tau_smooth=smooth(tau0,1);       
        [min1]=islocalmin(tau_smooth);
        rfw=output.analysis.scan_freq(mri);
        idmin=find(min1==1);
        il=length(idmin);
        if il>=3
        tau0=tau0(idmin(1):idmin(il));
        rfw=rfw(idmin(1):idmin(il));
        ffcf0=ffcf0(idmin(1):idmin(il),:);
      
        end
        tau_smooth=smooth(tau0,1);  
        [max1]=islocalmax(tau_smooth);
        id1=find(max1==1);
        id2=find(min1==1);
        if length(id1)>1 && length(id2)<=length(id1)
        id=id1([1,length(id1)]);
        elseif length(id1)==1
            id=id1(1);
        elseif length(id1)>1 && length(id2)==length(id1)
            if id1(1)<id2(2) && id1(1)>id2(1)
            id=id1([length(id1)]);
            else
                id=id1(1);
            end
        elseif length(id2)>1 && length(id2)>=length(id1)
            id=id1;
        else
            id=1;
            warning('No local maximum found.');
        end

        output.results.FFCF=ffcf0(id,:);
        output.results.tau_in_ps=tau0(id);


   min2=id2(id2<max(id1) & id2>min(id1));
        if isempty(min2)
            
            min2=1;
        end

        output.results.peak_freq=rfw(min2);


%% Create figures for frequency resolved analysis and FFCF results.
        close([figure(101) figure(102) figure(103) figure(104)]);
        output.fig=[figure(101) figure(102)];
        
        %Frequency Resolved Analysis plot
        figure(output.fig(1)), line([rfw(id(1)) rfw(id(1))],[0 tau0(id(1))],'Color','b','LineWidth',2,'LineStyle',':'); 
        if length(id)>1
            a=length(id);
        hold on, line([rfw(id(a)) rfw(id(a))],[0 tau0(id(a))],'Color','r','LineWidth',2,'LineStyle',':'); hold off;
        end
        hold on; scatter(rfw,tau0,'k','filled');hold off
        ylabel('fitted \tau_c (ps)');
        xlabel('Center frequency (cm^-^1)');
        title ('Frequency Resolved Analysis');
        set(gca,'FontSize',14);
        box on;
        axis tight;
        set(gca,'Ylim',[min(tau0)-0.05*range(tau0) max(tau0)+0.05*range(tau0)]);
        set(gca,'Xlim',[min(rfw)-5 max(rfw)+5]);
        set(gcf, 'position',[figV figH figL figW]);


        % FFCF plot
        Tw_ps=Tw./1000;
        ffcf1=output.results.FFCF(1,:);
        figure(output.fig(2)), yyaxis left,scatter(Tw_ps,ffcf1,'filled','Marker','^','MarkerFaceColor','b','MarkerEdgeColor','b');
        ylabel('PZS (max normalized)');
        set(gca,'YColor','b');
        set(gca,'FontSize',14);
        set(gca,'Ylim',[min(ffcf1)-0.05*range(ffcf1) max(ffcf1)+0.05*range(ffcf1)]);
        set(gca,'Ylim',[min(ffcf1)-0.05*range(ffcf1) max(ffcf1)+0.05*range(ffcf1)]);
        if length(id)>1
            ffcf2=output.results.FFCF(length(id),:);
                    yyaxis right, scatter(Tw_ps,ffcf2,'filled','Marker','s','MarkerFaceColor','r','MarkerEdgeColor','r');
                    set(gca,'YColor','r');
                    ylabel('PZS (max normalized)');
                    
        xlabel('waiting time (ps)');
        title 'FFCF approximation';
        set(gca,'FontSize',14);
        box on;
        axis tight;
        legend('Low Frequency','High Frequency');
        yyaxis left, set(gca,'Ylim',[min(ffcf1)-0.05*range(ffcf1) max(ffcf1)+0.05*range(ffcf1)]);
        yyaxis right, set(gca,'Ylim',[min(ffcf2)-0.05*range(ffcf2) max(ffcf2)+0.05*range(ffcf2)]);
        set(gca,'Xlim',[min(Tw_ps)-0.25 max(Tw_ps)+0.25]);
        set(gcf, 'position',[figV figH figL figW]);
        else

        xlabel('waiting time (ps)');
        title 'FFCF approximation';
        set(gca,'FontSize',14);
        box on;
        axis tight;

        set(gca,'Ylim',[min(ffcf1)-0.05*range(ffcf1) max(ffcf1)+0.05*range(ffcf1)]);
        set(gca,'Xlim',[min(Tw_ps)-0.25 max(Tw_ps)+0.25]);
        set(gcf, 'position',[figV figH figL figW]);     
        
        end


npf=length(output.results.peak_freq);
taus=output.results.tau_in_ps;
tl=length(taus);

if tl>1
    diff_percent=(max(taus)-min(taus))/min(taus);
    if diff_percent==0
        diff_percent=1;
    end
else
    diff_percent=1;
end

    if npf == 1 && diff_percent > 0.25
    
        %% High Frequency Analysis
        dr=size(Sp,4);
            ind=13:dr-3;
            if length(taus)>1
                if taus(2)>taus(1)
                    df=0; % delta frequency from center (extreme overlapping features of spectral bands)
                else
                    df=round(anharm/3);
                end
            else
                df=0;
            end
            pd=2*df;
 
            rl_hf=ceil(length(col_w)/2)-pd;
            if rl_hf<=0
                pd=0;
                rl_hf=ceil(length(col_w)/2)-pd;
            end
            pc_hf=PZP_calc(rl_hf,norder);
            
            Sp_hf=Sp(mp-floor(mp/2)+1+pd/2:mp+ceil(mp/2)-pd/2,mp+pd:end,:,ind);
            wt_hf=wt_c(mp+pd:end);
            pz_hf=PZM(Sp_hf,pc_hf);
            [~,n,p]=size(pz_hf);
            ffcf_hf=nan(p,n);
            tau_hf=nan(p,1);
                for k=1:p
                cp=similarity(pz_hf(:,:,k)');
                ffcf_hf(k,:)=cp;
                tau_hf(k,1)=TimeDecay(cp,Tw,nExpDec);
                end
           output.spectral_region.Pv_hf=pz_hf;
           output.spectral_region.Sp_hf=Sp_hf;
           output.spectral_region.wt_hf=wt_hf;
           output.analysis.FFCF_hf=ffcf_hf;
           output.analysis.tau_in_ps_hf=tau_hf;
           output.analysis.scan_freq_hf=wtau_s(ind);
    
           
      
    %% Low Frequency Analysis
           ind=1:(dr-13);
           if length(taus)>1
                if taus(2)>taus(1)
                    df=round(anharm/3); % delta frequency from center (extreme overlapping features of spectral bands)
                else
                    df=0;
                end
            else
                df=0;
           end
            pd=2*df;
            rl_lf=ceil(length(col_w)/2)-pd;
                
            if rl_lf<=0
                pd=0;
                rl_lf=ceil(length(col_w)/2)-pd;
            end
            pc_lf=PZP_calc(rl_lf,norder);
            Sp_lf=Sp(mp-floor(mp/2)+1+pd/2:mp+ceil(mp/2)-pd/2,1:mp-pd,:,ind);
            wt_lf=wt_c(1:mp-pd);
            pz_lf=PZM(Sp_lf,pc_lf);
            [~,n,p]=size(pz_lf);
            ffcf_lf=nan(p,n);
            tau_lf=nan(p,1);
                for k=1:p
                cp=similarity(pz_lf(:,:,k)');
                ffcf_lf(k,:)=cp;
                tau_lf(k,1)=TimeDecay(cp,Tw,nExpDec);
                end
           output.spectral_region.Pv_lf=pz_lf;
           output.spectral_region.Sp_lf=Sp_lf;
           output.spectral_region.wt_lf=wt_lf;
           output.analysis.FFCF_lf=ffcf_lf;
           output.analysis.tau_in_ps_lf=tau_lf;
           output.analysis.scan_freq_lf=wtau_s(ind);
    
    %% Frequency Window Selection
            % Low
            taul_smooth=smooth(tau_lf,1);
            [maxl]=islocalmax(taul_smooth);
            [minl]=islocalmin(taul_smooth);
            rfwl=output.analysis.scan_freq_lf;
            idl=find(maxl==1);
            if isempty(idl)
                idl=1;
            end
           
            output.results.FFCF_lf=ffcf_lf(idl(1),:);
            output.results.tau_in_ps_lf=tau_lf(idl(1));
            output.results.peak_freq_lf=rfwl(minl);
    
            
            
            % High
                    tauh_smooth=smooth(tau_hf,1);
            [maxl]=islocalmax(tauh_smooth);
            [minl]=islocalmin(tauh_smooth);
            rfwh=output.analysis.scan_freq_hf;
            id=find(maxl==1);
            if length(id)>=2
            idh=id(length(id));
            output.results.FFCF_hf=ffcf_hf(idh,:);
            output.results.tau_in_ps_hf=tau_hf(idh);
            else
            idh=id(1);
            output.results.FFCF_hf=ffcf_hf(idh,:);
            output.results.tau_in_ps_hf=tau_hf(idh);
            end
            output.results.peak_freq_hf=rfwh(minl);
    %% Create figures for in-depth frequency resolved analysis and FFCF results.
            close([figure(103) figure(104)]);
            output.fig2=[figure(103) figure(104)];
            
            %Frequency Resolved Analysis plot
            tau_sh=output.analysis.tau_in_ps_hf;
            tau_sl=output.analysis.tau_in_ps_lf;
            figure(output.fig2(1)), line([rfwl(idl(1)) rfwl(idl(1))],[-1 tau_sl(idl(1))],'Color','b','LineWidth',2,'LineStyle',':'); 
            hold on; line([rfwh(idh(1)) rfwh(idh(1))],[-1 tau_sh(idh)],'Color','r','LineWidth',2,'LineStyle',':'); hold off;
            hold on; scatter(rfwl,tau_sl,'b','filled','Marker','^');hold off
            hold on; scatter(rfwh,tau_sh,'r','filled','Marker','s');hold off
            ylabel('fitted \tau_c (ps)');
            xlabel('Center frequency (cm^-^1)');
            title ('Frequency Resolved Analysis');
            set(gca,'FontSize',14);
            box on;
            axis tight;
            set(gca,'Ylim',[min([tau_sh;tau_sl])-0.05*range([tau_sh;tau_sl]) max([tau_sh;tau_sl])+0.05*range([tau_sh;tau_sl])]);
            set(gca,'Xlim',[min([rfwl;rfwh])-1 max([rfwl;rfwh])+1]);
            set(gcf, 'position',[figV figH figL figW]);
    
    
            % FFCF plot
            Tw_ps=Tw./1000;
            ffcf1=output.results.FFCF_lf;
            figure(output.fig2(2)),yyaxis left, scatter(Tw_ps,ffcf1,'filled','Marker','^','MarkerFaceColor','b','MarkerEdgeColor','b');
            ylabel('PZS (max normalized)');
            set(gca,'YColor','b')
            set(gca,'Xlim',[min(Tw_ps)-0.5 max(Tw_ps)+0.5]);
            set(gca,'FontSize',14);
    
            ffcf2=output.results.FFCF_hf;
            yyaxis right, scatter(Tw_ps,ffcf2,'filled','Marker','s','MarkerFaceColor','r','MarkerEdgeColor','r');
            set(gca,'YColor','r')
            ylabel('PZS (max normalized)');
            xlabel('waiting time (ps)');
            title 'FFCF approximation';
            set(gca,'FontSize',14);
            box on;
            axis tight;
            legend('Low Frequency','High Frequency');
      
            yyaxis left, set(gca,'Ylim',[min(ffcf1)-0.05*range(ffcf1) max(ffcf1)+0.05*range(ffcf1)]);
            yyaxis right, set(gca,'Ylim',[min(ffcf2)-0.05*range(ffcf2) max(ffcf2)+0.05*range(ffcf2)]);
            set(gca,'Xlim',[min(Tw_ps)-0.25 max(Tw_ps)+0.25]);
            set(gcf, 'position',[figV figH figL figW]);
            
    end
       
    otherwise
        error('Analysis method is not selected correctly.')
  
end

end

%% Function to calculate pseudo-Zernike moment
function [M]= PZM(spectra,PZcoeff)
        fl=size(spectra,1);
        N=fl-1;
        x = (-N:2:N)/N;
        [X,Y] = meshgrid(x);
        [~,r] = cart2pol(X,Y);
  
        is_in_circle = r <= 1;
        n_Tw=size(spectra,3);
        iic=reshape(is_in_circle,fl^2,1);
        if ndims(spectra)==3
        
        
            F=reshape(spectra,fl^2,n_Tw);
            F(~iic,:) = nan;
            Fn=F;   
            Fc=Fn(iic,:);
            M=Fc'*PZcoeff;
        elseif ndims(spectra)==4
            n_scan=size(spectra,4);
            F=reshape(spectra,fl^2,n_Tw,n_scan);
            F(~iic,:,:) = nan;
            Fn=F;   
            Fc=Fn(iic,:,:);
            Pc=repmat(PZcoeff',[1 1 n_scan]);
            M=pagemtimes(Pc,Fc);
        end
       
end

%% Function to calculate pseudo-Zernike polynomials
function [PZcoeff]=PZP_calc(Nsize,n_max)
    N=Nsize-1;
    x = (-N:2:N)/N;
    [X,Y] = meshgrid(x);
    [theta,r] = cart2pol(X,Y);
    is_in_circle = r <= 1;
    r = r(is_in_circle);
    theta = theta(is_in_circle);
        n = zeros(1,(n_max+1)^2);
        m = zeros(1,(n_max+1)^2);
        for k = 0:n_max
            n(k^2+1:(k+1)^2) = repmat(k,1,2*k+1);
            m(k^2+1:(k+1)^2) = -k:k;
        end
        PZcoeff = pzernfun(n,m,r,theta);
end

%% Function to calculate pseudo-Zernike polynomial at each pixel
function p = pzernfun(n,m,r,theta)
if ( ~any(size(n)==1) ) || ( ~any(size(m)==1) )
    error('pzernfun:NMvectors','N and M must be vectors.')
end

if length(n)~=length(m)
    error('pzernfun:NMlength','N and M must be the same length.')
end

n = n(:);
m = m(:);

if any(m>n)
    error('pzernfun:MlessthanN', ...
          'Each M must be less than or equal to its corresponding N.')
end

if any( r>1 | r<0 )
    error('pzernfun:Rlessthan1','All R must be between 0 and 1.')
end

if ( ~any(size(r)==1) ) || ( ~any(size(theta)==1) )
    error('pzernfun:RTHvector','R and THETA must be vectors.')
end

r = r(:);
theta = theta(:);
length_r = length(r);
if length_r~=length(theta)
    error('pzernfun:RTHlength', ...
          'The number of R- and THETA-values must be equal.')
end

m_abs = abs(m);
rpowers = [];
for j = 1:length(n)
    rpowers = [rpowers m_abs(j):n(j)]; %#ok<AGROW>
end
rpowers = unique(rpowers);

if rpowers(1)==0
    rpowern = arrayfun(@(p)r.^p,rpowers(2:end),'UniformOutput',false);
    rpowern = cat(2,rpowern{:});
    rpowern = [ones(length_r,1) rpowern];
else
    rpowern = arrayfun(@(p)r.^p,rpowers,'UniformOutput',false);
    rpowern = cat(2,rpowern{:});
end

p = zeros(length_r,length(n));
for j = 1:length(n)
    s = 0:(n(j)-m_abs(j));
    pows = n(j)-s;
    for k = length(s):-1:1
        c = (1-2*mod(s(k),2))*                 ...
               prod(2:2*n(j)+1-s(k))/          ...
               prod(2:s(k))/                   ...
               prod(2:(n(j)+m_abs(j)+1-s(k)))/ ...
               prod(2:(n(j)-m_abs(j)  -s(k)));
        idx = (pows(k)==rpowers);
        p(:,j) = p(:,j) + c*rpowern(:,idx);
    end
    p(:,j) = p(:,j)*sqrt((n(j)+1)*(1+(m_abs(j)~=0))/pi);
end

idx_pos = m>0;
idx_neg = m<0;

if any(idx_pos)
    p(:,idx_pos) = p(:,idx_pos).*cos(theta*m_abs(idx_pos)');
end
if any(idx_neg)
    p(:,idx_neg) = p(:,idx_neg).*sin(theta*m_abs(idx_neg)');
end

end

%% Calculate pseudo-Zernike similarity (cosine)
function [Wavg]=similarity(C)
wp=[cosd(45) 1-cosd(45)];
Cn=C./vecnorm(C,2,2);
nT=size(C,1);
mi=Cn(1,:);
mf=Cn(nT,:);
CS=nan(nT,1);
CD=nan(nT,1);
        for j=1:nT       
        mt=Cn(j,:);
        CS(j,1)=dot(mi,mt);
        CD(j,1)=(1-dot(mt,mf));
        end
W=wp(1)*CS+wp(2)*CD;
Wavg=W./(max(W));
end

%% Function to calcualte time decay constant with exponential decay fit.
function [tau_c]=TimeDecay(SM,t,model)
ExpDec1=fittype('A*exp(-t/tau)+y0','independent','t');
    if ~isnan(SM)
        y=SM;    
        t=t./1000;
                    ye=min(y);
            ys=max(y);
        y_r=(y-ye)/(ys-ye); 
        if model==1
        
            try
            fx = fit(t,y,ExpDec1,'Algorithm','Levenberg-Marquardt','StartPoint',[max(y)-min(y) 1 min(y)]);         
            catch
                warning off backtrace;
                warning('Exponential Decay Fit did not converge at all frequencies while scanning.');
                
                fx.tau=0;
            end
            tau_c=fx.tau;
            if tau_c >2*max(t) || tau_c<0
                tau_c=nan;
            end
        
        else
 
        tau_c=trapz(t,y_r);
        end
    else
        tau_c=nan;
    end
end
