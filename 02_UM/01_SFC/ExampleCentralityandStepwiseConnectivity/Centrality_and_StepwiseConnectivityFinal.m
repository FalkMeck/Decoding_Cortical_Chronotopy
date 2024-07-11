% Based on Sepulcre et al, Journal of Neuroscience 2012
% Other works based on SFC from the Lab: Sepulcre, Cerebral Cortex 2015; Ortiz-Teran et al, PNAS 2017; Ibai & Sepulcre, Nat Commun 2018; Bueichekú et al, PNAS 2020.

fmri_file='100110_HU67UU_sub4.nii.gz';
mask_file='maskGMtwohemis_sub4.nii.gz';
seed_file='V1_sub4.nii.gz';

q_fdr=0.001; % q for the fdr
num_steps=7; % number of steps for the SFC

%%Read brain mask
mask=MRIread(mask_file);
mask=mask.vol;
mask_indx=find(reshape(mask,[],1)>0);
num_rois=length(mask_indx);

%%Read brain seed
seed=MRIread(seed_file);
seed=seed.vol;
seed=reshape(seed(mask_indx),[],1);
seed_indx=find(seed>0);

%%Read fmri file, convert to 2d matrix and apply brain mask
fmri=MRIread(fmri_file);
fmri=fmri.vol;
fmri=reshape(fmri,[],size(fmri,4)); %convert to 2d
fmri=fmri(mask_indx,:); %mask the data

%%Compute functional matrix
[fc,pvals]=corr(fmri');
fc(1:(1+num_rois):end)=0; %diagonal to zero
fc= 0.5 * log( ( 1 + fc)./( 1 - fc) ); %fisher (optional)

%%Take positive links and filter using fdr (all optional but recommended)
fc(find(fc<0))=0;
indx=find(fc>0);
[n_signif,indx_signif]=fdr(pvals(indx),q_fdr,'original','mean');
fdr_matrix=zeros(size(pvals));
fdr_matrix(indx(indx_signif))=1;
fdr_matrix=uint8(fdr_matrix);
adj=fc.*double(fdr_matrix);

%%Centrality & Stepwise
adj=(adj - min(adj(:)))./ (max(adj(:)) - min(adj(:))); %normalize (optional)
[m n]=size(adj);
adj(1:(m+1):end)=0; %diagonal to zero
step1=adj;
degree_centrality=sum(step1,2);
adjp=step1;
step1_seed=adjp(seed_indx,:);
for s=2:num_steps
  adj=adjp*adjp; % for exponential increase, multiply by output; for arithmetic increase multiply always by step1
  adj=(adj - min(adj(:)))./ (max(adj(:)) - min(adj(:))); %normalize (optional)
  adj(1:(m+1):end)=0; %diagonal to zero
  adjp=adj;
  eval(sprintf('step%s=adjp;',num2str(s)));
  eval(sprintf('step%s_seed=adjp(seed_indx,:);',num2str(s)));
end

hdr0=MRIread(mask_file); 
img=zeros(size(mask));
img(mask_indx)=degree_centrality;
hdr0.vol=img;
hdr0.fspec='degree_centrality.nii.gz';
MRIwrite(hdr0,hdr0.fspec);  
clear hdr0 img 

hdr0=MRIread(mask_file); 
img=zeros(size(mask));
img(mask_indx)=step1_seed;
hdr0.vol=img;
hdr0.fspec='SFC_1_V1seed.nii.gz';
MRIwrite(hdr0,hdr0.fspec); 
clear hdr0 img 

hdr0=MRIread(mask_file); 
img=zeros(size(mask));
img(mask_indx)=step2_seed;
hdr0.vol=img;
hdr0.fspec='SFC_2_V1seed.nii.gz';
MRIwrite(hdr0,hdr0.fspec); 
clear hdr0 img 

hdr0=MRIread(mask_file); 
img=zeros(size(mask));
img(mask_indx)=step3_seed;
hdr0.vol=img;
hdr0.fspec='SFC_3_V1seed.nii.gz';
MRIwrite(hdr0,hdr0.fspec);
clear hdr0 img  

hdr0=MRIread(mask_file); 
img=zeros(size(mask));
img(mask_indx)=step4_seed;
hdr0.vol=img;
hdr0.fspec='SFC_4_V1seed.nii.gz';
MRIwrite(hdr0,hdr0.fspec); 
clear hdr0 img 

hdr0=MRIread(mask_file); 
img=zeros(size(mask));
img(mask_indx)=step5_seed;
hdr0.vol=img;
hdr0.fspec='SFC_5_V1seed.nii.gz';
MRIwrite(hdr0,hdr0.fspec); 
clear hdr0 img 

hdr0=MRIread(mask_file); 
img=zeros(size(mask));
img(mask_indx)=step6_seed;
hdr0.vol=img;
hdr0.fspec='SFC_6_V1seed.nii.gz';
MRIwrite(hdr0,hdr0.fspec);
clear hdr0 img  

hdr0=MRIread(mask_file); 
img=zeros(size(mask));
img(mask_indx)=step7_seed;
hdr0.vol=img;
hdr0.fspec='SFC_7_V1seed.nii.gz';
MRIwrite(hdr0,hdr0.fspec);
clear hdr0 img  

