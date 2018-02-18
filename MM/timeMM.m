diary('timing_results')
t = linspace(0,20000,11)';
t1 = cputime;
num = 10000;
for i = 1:num
    %MichMent( timepoints, species_init, parameters, suppress_plot )
    %function description
    [~,~, species_out, ~]=MichMent(t, [],[1,-2.77,-1,-2],1 );
end
t2 = cputime;
t2-t1
diary off
%%
dlmwrite('mex_results.txt',species_out,'delimiter','\t','precision',16);