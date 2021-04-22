function [List] = find_boundary_particle(file_name,draw)
% This function use the bw image of Wh image to identify the boundary
% regions and find the particles that are closest to the boundary region.
% List is a 1xN array with 1 (0) if the particle is (or not) a boundary
% particle.

%% PART I. Read in the particle position.
P = csvread(['Processing/Position_corrected/' file_name '.txt']);
Im = imread(['Picture/' file_name '_Wh.jpg']);
Im = 1/2*im2double(Im(:,:,1))+1/2*im2double(Im(:,:,2));
% Initialization of particle radius.
R_b=61.5;
R_s=R_b*6.35/7.95;
List = ones(1,length(P));
%% PART I. Creat the boundary bw image.
% B = im2double(Im);
% B = imadjust(B);
% Walls = bwareaopen(imdilate(~imbinarize(B, 0.2),strel('disk',20)), 10^5);
% Walls = ~bwareaopen(imcomplement(Walls), 10^7);
% A_pin(Walls)=0;
[A_pin, Walls] = Correction_yz(Im,2);
if draw == 1
figure,imshow(A_pin);
end
%% PART II. Identify the boundary particles.

% First for the left and right side, just use the c position is enough.
x_limit = [635 5235];
% For the upper and lower side, need to
bw_dist = bwdist(Walls);
for i = 1:length(P)
    if P(i,1)~=0
        temp = bw_dist(round(P(i,2)),round(P(i,1)));
        if P(i,1)<x_limit(1) || P(i,1)>x_limit(2) || temp<2*R_b
            List(i) = 0;
            if draw == 1
                if P(i,3) == 1
                    Rp = R_b;
                else
                    Rp = R_s;
                end
                hold on;
                plot(P(i,1)+Rp*cos(0:0.01:2*pi),P(i,2)+Rp*sin(0:0.01:2*pi),'g-','LineWidth',2);
            end
        end
    end
end
