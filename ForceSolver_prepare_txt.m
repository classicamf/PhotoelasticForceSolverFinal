function ForceSolver_prepare_txt(file_name)
% This function based on the results of balance_reduce solutions, stored in the
% txt file.
tic;
% load the mat file for the neural network prediction
Imax=0.61;
Imin=0.09;

load('C:/Personal/Dukework/Linear/2020Data/ML4-m3-ff.mat');
PL = imread(['Processing/Adj_PL_3/' file_name '_3_PLonly.png']);
PL = im2double(PL);
if exist('Processing/Inverse','dir')~=7
    mkdir Processing/Inverse
    mkdir Processing/Inverse/SolvedMat
end
contact_thresh = 0;
[Cmatrix,Cneighbor,~,~,~]=linear_contact_new_3(file_name,contact_thresh,0);

%% Read in the information
P = csvread(['Processing/Position_corrected/' file_name '.txt']);
P(P(:,1)==0,:) = [];
%Np is the number of particles
Np=length(P);
% Create the matrix to store the data.
Mat = zeros(5*Np,10);
fsigma = 100; %photoelastic stress coefficient
%The radius of the particles
R_b=61.5;
R_s=R_b*6.35/7.95;
Meter_per_pixel = 7.95*10^(-3)/R_b ;% Real meter per pixel. Currently not used.
% Create two square masks for big and small particles respectively.
[x,y]=meshgrid(-R_b:R_b,-R_b:R_b);
% D records the distance to center from each pixels
Dis=sqrt(x.^2+y.^2);
% Ds and Db are masks for small and big particles.
Ds=Dis;
Ds(Dis>R_s)=0;
Ds(Dis<=R_s)=1;
Db=Dis;
Db(Dis>R_b)=0;
Db(Dis<=R_b)=1;
% add the particle indexes
for i=1:Np
    Mat(i,1) = i;
    Mat(Np+i,1) = i;
    Mat(2*Np+i,1) = i;
    Mat(3*Np+i,1) = i;
    Mat(4*Np+i,1) = i;
    Mat(i,1) = i;
    Mat(i,2) = P(i,1);
    Mat(i,3) = P(i,2);
    Mat(i,4) = P(i,3);
    if P(i,3)==0
        %This is a small particle
        Mat(i,5) = R_s;
        Mat(i,6) = R_s * Meter_per_pixel;
        mask = Ds;
    else
        Mat(i,5) = R_b;
        Mat(i,6) = R_b * Meter_per_pixel;
        mask = Db;
    end
    Mat(i,7) = fsigma;
    Mat(i,8) = Cmatrix(i,1);
    Mat(i,9) = Cmatrix(i,2);
    for k=1:Mat(i,8)
        %Loop all the contacts for ith particle.
        neighbor_index=Cneighbor(i,k);
        Mat(4*Np+i,k+1) = neighbor_index;
        xk=P(neighbor_index,1);
        yk=P(neighbor_index,2);
        beta_k=getangle(Mat(i,2),Mat(i,3),xk,yk);
        Mat(3*Np+i,k+1) = beta_k;
    end
    %% Give the initial guess.
    %Calculating the particle image for this particle.
    cropXstart=round(Mat(i,2)-R_b);
    cropXstop=round(Mat(i,2)-R_b)+2*R_b;
    cropYstart=round(Mat(i,3)-R_b);
    cropYstop=round(Mat(i,3)-R_b)+2*R_b;
    cimg=PL(cropYstart:cropYstop,cropXstart:cropXstop);
    forceImage=cimg;
    forceImage = forceImage.*mask;
    %% Calculate radial_g2
    if Mat(i,8)>1
        radial_g2 = boundaryg2_new(forceImage,Mat(i,5),Mat(3*Np+i,2:Mat(i,8)+1),1,0);
        % radial_g2 is zx1 array
        for k = 1:Mat(i,8)
            if radial_g2(k)>0.5*10^(-3)
                %forces(i)=pres(n).radial_g2(i)*63.82-0.1;
                %% Use neural network to guess this force in particular
                % Note, mask1, k1, etc parameters are
                % loaded from the Mat file.
                [force_new,alpha_new] = neural_guess(forceImage,Mat(3*Np+i,2:Mat(i,8)+1),k,Imax,Imin,mask1,k1,b1,k2,b2,k3,b3);
                Mat(Np+i,k+1) = force_new;
                Mat(2*Np+i,k+1) = alpha_new;
            else
                Mat(Np+i,k+1) = 0.01;
                Mat(2*Np+i,k+1) = mod(Mat(3*Np+i,k+1)+pi,2*pi);
            end
        end
    end
end

csvwrite(['Processing/Inverse/SolvedMat/' file_name '_prepare.txt'],Mat);
toc;
