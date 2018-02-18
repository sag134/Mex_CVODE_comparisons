gcc = importdata('output_MichMent_gcc');
mex = importdata('mex_results.txt');
bng = importdata('MichMent.gdat');
bng = bng.data;
bng_rearranged = [bng(:,1),bng(:,2),bng(:,3),bng(:,5),bng(:,4)];

bng_mex_error = bng_rearranged(:,2:end)-mex;
bng_gcc_error = bng_rearranged-gcc;
mex_gcc_error = mex-gcc(:,2:end);
dlmwrite('bng_mex_error.txt',bng_mex_error,'\t')
dlmwrite('bng_gcc_error.txt',bng_gcc_error,'\t')
dlmwrite('mex_gcc_error.txt',bng_gcc_error,'\t')