function txt2fig_2DIR(fig_folder)
%%
    % Convert 2DIR .txt file into matlab figures at each waiting time
    % Save the generated matlab figures in a new folder (fig_folder) in the same directory as .txt file
%%
cd(uigetdir)
wtau=readmatrix('wtau_freq.txt');
wt=readmatrix('wt_freq.txt');
txt_list=dir('*Tw*.txt');

figure(222)
sp=mkdir(fig_folder);
    for n=1:length(txt_list)
        spectra=readmatrix(txt_list(n).name);
        contourf(wt,wtau,spectra);
        cf=pwd;
        cd(fig_folder);
        savefig([erase(txt_list(n).name,'.txt'),'.fig'])
        cd(cf);
    end
end
