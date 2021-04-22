function [C,Cneighbor,Znr,Z,fnr]=linear_contact_new_3(file_name,thresh,draw)
% This function calculate the contact for particles away from the
% boundary...
tic;
% thresh is used to identify stressed contacts.
% File name should be like: 0000_Forward for files like 0000_Forward_PL.jpg
List = find_boundary_particle(file_name,0);
% List records whether this particle is a boundary particle.
% First check whether should we create a new directory
if exist('Processing\Contact','dir')~=7
    mkdir Processing\Contact
end
%% PART -1: Identify the ROI within which we calculate the contact network information.
%c1 = 656;
%c2 = 5206;
%r1 = 786;
%r2 = 2913;

%% PART 0: Read in and pre-process the PL image.
disp('Initialization');
%temp = matfile(['Processing/Adj_PL_3/' file_name '_3_PLonly.mat']);
%PL = temp.PL_corrected2;

PL = imread(['Processing/Adj_PL_3/' file_name '_3_PLonly.png']);
PL = im2double(PL);

%% PART I: Initializing parameters
% P contains the positions for corrected images.
P = csvread(['Processing/Position_corrected/' file_name '.txt']);
P(P(:,1)==0,:) = [];
% thresh is the intensity thresh used to identify contacts.
%thresh = 0.5;
%thresh = 0;
%R_s and R_b are used to determine the boundary pairs.
%coef is determined by a test
%coef=1.0825;
%coef=1.0825;
coef=1.03;
R_b=61.5*coef;
R_s=R_b*6.35/7.95;
%dR is the radius of the ROI, determined by visual estimation.
dR=10;
%n is the number of particles
n=size(P,1);
%Rparticle should be exactly the particle radius.
%Rextend is the radius of the particle used to determine geometric
%neighbors.

Rparticle=zeros(1,n);
Rextend=zeros(1,n);
for i=1:n
    if P(i,3)==0
        %Now this particle is a small particle.
        Rparticle(i)=round(R_s/coef);
        Rextend(i)=round(R_s);
    else
        Rparticle(i)=round(R_b/coef);
        Rextend(i)=round(R_b);
    end
end

%G=im2double(PL(620:3004+50,410:5245+100+100,2));
[~,f] = linear_g2particle_new(PL,P,1);
if draw == 1
    figure,imagesc(PL,[0,1.5]);
    %figure,imagesc(log(f),[-15,-5]);
    axis image;
    colormap(jet);
    hold on;
end
%%%%%%%%
%C is a matrix contain contact number information. C(i,1) is the number
%of mechanical contact of the ith particle. C(i,2) is 0 or 1. C(i,2)=0
%means it is a boundary particle. C(i,2)=1 means it is a particle of
%interest.
C=zeros(n,2);

for i=1:n
    %if P(i,1)>c1 && P(i,1)<c2 && P(i,2)>r1 && P(i,2)<r2
        %C(i,2) = 1;
    %end
    C(i,2) = List(i);
    if draw == 1
        %text(P(i,1),P(i,2),num2str(i));
        plot(P(i,1)+Rextend(i)*sin(0:0.01:2*pi),P(i,2)+Rextend(i)*cos(0:0.01:2*pi),'k-');
        plot(P(i,1)+Rparticle(i)*sin(0:0.01:2*pi),P(i,2)+Rparticle(i)*cos(0:0.01:2*pi),'k-');
        if C(i,2)==0
            %plot(P(i,1)+Rparticle(i)*sin(0:0.01:2*pi),P(i,2)+Rparticle(i)*cos(0:0.01:2*pi),'w-','LineWidth',1);
        else
            %plot(P(i,1)+Rparticle(i)*sin(0:0.01:2*pi),P(i,2)+Rparticle(i)*cos(0:0.01:2*pi),'g-','LineWidth',1);
        end
    end

end

%Cneighbor is a matrix contain neighbor information for each particle.
Cneighbor=zeros(n,10);
%Int_ROI is used to store the intensity inside the ROI
%G2_ROI is used to store the G2 inside the ROI
%Int_ROI=zeros(n,10);
G2_ROI=zeros(n,10);

%Dist_pair contains all pair distance.
Dist_pair=zeros(n,n);
for i=1:n
    Dist_pair(i,:)=sqrt((P(i,1)-P(:,1)).^2+(P(i,2)-P(:,2)).^2);
end

disp('Working on calculating contacts');

if draw == 1
    ROI_info = [];
end

%Begin to find possible contacts.
for i=1:n-1
    for j=i+1:n
        %For the i-j pair
        if Dist_pair(i,j)<Rextend(i)+Rextend(j)
            
            %(xi,yi) is the center of i ROI, (xj,yj) is the center of j
            %ROI.
            thi=getangle(P(i,1),P(i,2),P(j,1),P(j,2));
            thj=getangle(P(j,1),P(j,2),P(i,1),P(i,2));
            %Get the calculated
            %Int_i=getroi(P(i,1),P(i,2),G,thi,Rparticle(i),dR);
            %Int_j=getroi(P(j,1),P(j,2),G,thj,Rparticle(j),dR);
            Int_i=getroi_20(P(i,1),P(i,2),PL,thi,Rparticle(i),dR,0);
            Int_j=getroi_20(P(j,1),P(j,2),PL,thj,Rparticle(j),dR,0);
            
            if draw == 1
                ROI_info = [ROI_info Int_i];
                ROI_info = [ROI_info Int_j];
            end
            %If the center distance is smaller than the condition.
            %C_i and C_j are index for i,j ROI
            C_i=0;
            C_j=0;
            %if Int_i-ROI_initial(i)>thresh
            if Int_i>thresh
                C_i=1;
            end
            %if Int_j-ROI_initial(j)>thresh
            if Int_j>thresh
                C_j=1;
            end
            if C_i*C_j==1
                %Calculate the G2 of ROI for each contact.
                g2_i=getroi_20(P(i,1),P(i,2),f,thi,Rparticle(i),dR,0);
                g2_j=getroi_20(P(j,1),P(j,2),f,thj,Rparticle(j),dR,0);
                %Add contact number.
                C(i,1)=C(i,1)+1;
                C(j,1)=C(j,1)+1;
                %Record the neighbor
                Cneighbor(i,C(i,1))=j;
                Cneighbor(j,C(j,1))=i;
                if draw == 1
                    plot([P(i,1) P(j,1)],[P(i,2) P(j,2)],'r--','LineWidth',1);
                end
                
                %Record ROI int
                %Int_ROI(i,C(i,1))=Int_i;
                %Int_ROI(j,C(j,1))=Int_j;
                %Record ROI g2
                G2_ROI(i,C(i,1))=g2_i;
                G2_ROI(j,C(j,1))=g2_j;
            end
        end
    end
end


CC=C(C(:,2)==1);
Znr=mean(CC(CC>=1));
Z=mean(CC(:,1));
fnr=numel(CC(CC>=1))/numel(CC);

disp(['Znr = ' num2str(Znr)]);
disp(['fnr = ' num2str(fnr)]);

if draw == 1
    [x,pdf,pdfuncertain,mx,hist]=ringpdf(ROI_info,20,max(ROI_info),min(ROI_info));
    figure,plot(x,pdf,'o');
end

toc;
% Save the information for individual particles...
%csvwrite(['Processing\Contact\' file_name '_Cn.txt'],Cneighbor);
%csvwrite(['Processing\Contact\' file_name '_C.txt'],C);
%csvwrite(['Processing\Contact\' file_name '_G2_ROI.txt'],G2_ROI);
