% Sparse Coding with Anomaly Detection
% Amir Adler · Michael Elad · Yacov Hel-Or ·Ehud Rivlin

clear all
close all

load('data/MITDB.mat');

sig109  = data{signal_numbers == 109}(:,1);
Fs109   = Fs{signal_numbers == 109};
ann109  = ann{signal_numbers == 109};
type109 = type{signal_numbers == 109};
sptype109 = [repmat('  ', size(type109,1),1) type109];

[b,a] = butter(5,[1 100]/Fs109*2, 'bandpass');
x = filtfilt(b, a, sig109);
plot([sig109, x])

plot(sig109);hold on;grid on
plot(ann109, sig109(ann109), 'ro');
text(ann109, sig109(ann109), sptype109);

% Take first 5 minutes
seglen = Fs109*300;
segnum = 3;
segidx = (segnum-1)*seglen+1:segnum*seglen;  % indices in x
%
xseg   = x(segidx);                          % restricted x
annidx = ann109>segidx(1)& ann109 <= segidx(end);
annseg = ann109(annidx);
annseg = annseg - segidx(1);
typeseg   = type109(annidx);
sptypeseg = sptype109(annidx);
%
plot(xseg);hold on;grid on
%plot(annseg, xseg(annseg), 'ro');
%text(annseg, xseg(annseg),sptypeseg);
%
annsegP    = annseg(typeseg~='L');
sptypesegP = sptypeseg(typeseg~='L');
plot(annsegP, xseg(annsegP), 'ro');
text(annsegP, xseg(annsegP),sptypesegP);

% Get windows of segment
winlen = 256;
wincount = length(xseg) - winlen + 1;
win = zeros(winlen, wincount);
for i = 1:wincount
    win(:,i) = xseg(i : i+winlen-1);
end
P = pca(win');
PCAsize = 32;
Pr = P(:,1:PCAsize);

%plot(Pr*Pr'*win(:,1:10))
Cwin = Pr'*win;
% 
% %% Dictionary learning
% dltrainset = Cwin;
% 
% dlparam.mode   = 2;
% dlparam.K      = 128; % codebook_size;
% dlparam.lambda = 0.5;
% dlparam.iter   = 50;
% %dlparam.pos    = true;
% dlparam.verbose = false;
% dlparam.modeD   = 3;
% %dlparam.posD    = true;
% dlparam.gamma1  = 0.5;
% dlparam.batchsize = 512;
% dlparam.clean  = true;
% %D0 = [dctmtx(NFEAT)];
% %dlparam.D = D0;
% fprintf('%s --- --- Using %d vectors for dictionary learning\n',datetime('now'), size(dltrainset,2));
% D=mexTrainDL(dltrainset,dlparam); 
% if any(any(isnan(D)))
%     error('Error in Dict Learning!')
% end
% 
% %Normalize dictionary
% for i = 1:size(D,2)
%     D(:,i) = D(:,i) / norm(D(:,i),2);
% end
%     
% %% Sparse coding of all data
% fprintf('%s --- --- Encoding feature vectors\n',datetime('now'));        
% % Lasso
% % param.mode = 2;
% % %param.lambda=3.9;
% % param.lambda = 0.5;
% % %param.pos = true;
% % SCcoef = mexLasso(Cwin, D, param);
% % OMP
% clear param
% param.L  = 4;
% SCcoef = mexOMP(Cwin, D, param);
% 
% % AA = py.numpy.array([[2,3] [2 3]])
% % BB = 
% % CC =
% % SCcorf = py.GLS.greedy_least_squares([7 7 7; 1 2 3],7,7)
% SCcorf = py.GLS.greedy_least_squares(toNdarray(Cwin), toNdarray(D), 14);

%% Dictionary learning
dltrainset = Cwin;
K = 128;
iter = 50;
L = 4;
%rp = randperm(size(dltrainset,2));
%D0 = dltrainset(:,rp(1:K));
D0 = dltrainset*randn(size(dltrainset,2), K);
for i = 1:size(D0,2), D0(:,i) = D0(:,i)/norm(D0(:,i)); end
% K-SVD with GLSP
[DGLSP, residGLSP, DGLSP_hist, coefGLSP_hist, GLSP_coef_hist, GLSP_supp_hist] = K_SVD(true, dltrainset, D0, iter, L);
np_result = cell(py.GLS.greedy_least_squares(mat2np(Cwin), mat2np(DGLSP), L));
GLSP_coef = np2mat(np_result{1});
GLSP_supp = np2mat(np_result{2});
GLSP_debias = np2mat(np_result{3});
pca_recerror2   = np2mat(np_result{4});
% K-SVD with OMP
[DOMP, residOMP, DOMP_hist, coefOMP_hist, ~, ~] = K_SVD(false, dltrainset, D0, iter, L);
clear param
param.L  = 4;
SCcoef = mexOMP(Cwin, DOMP, param);
pca_recerror1 = sqrt(sum((Cwin - DOMP*SCcoef).^2));
     
        
figure
imshow(full(SCcoef(:, 1:236:100000)))
figure
imshow(full(GLSP_coef(:, 1:236:100000)))
figure
imshow(win(:, 1:256:100000))

figure, plot(recerror1)
figure, plot(recerror2)

% Check reconstruction error
winRec = Pr*D*SCcoef;
recerror = (sum((win - winRec).^2));
figure
%plot(recerror(1:256:end))
plot(recerror(1:end))