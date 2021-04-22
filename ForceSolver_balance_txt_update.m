function ForceSolver_balance_txt_update(file_name,verbose,file_name_previous,display_setting, iteration_ratio,update_based_on_manual)
% This function based on the results of balance_reduce solutions, stored in the
% txt file.

% It uses previous step result (as well as reaction forces) as initial guesses, if the current result is
% not good.


% INPUTS:
% 1. file_name: for example, it should be in the form of '0000_Forward'
% 2. verbose should be an array like [0 0 0 0 1 0 0 0 0]. verbose(i) == 1
% means we will show the fitting process visually for particles with
% contact number i.
% 3. file_name_previous are in the form of '0000_Forward'
% 4. display_setting = 1 means many details associated with the fitting
% process will be displayed on the GUI. Otherwise it will not display
% anything.
% OUTPUTS:
% 1. Ave_f is the averaged contact forces for non-rattler particles
% 2. A mat file is  saved.
% INVOKED EXTERNAL FUNCTIONS
% 1. remove_contact : defined inside this function
% 2. fitforce_nobalance_new : external
% 3. get_corresponding_balance : external
% 4. fitforce_balance_new : external

% Thresholds used to find the particles need to update.
Err_thresh = 0.12;
% F_thresh is used to evaluate the bad action-reaction forces.
F_thresh = 0.07;
% In Newton, all contacts below this threshold will be droped.
F_th = 0.005 * 0.6368;
% Display the current working step file name
disp(['Work Load =  Step = ' file_name]);

%% May subject to change parameters

%The experimental light environment, Imax is the maximum intensity and Imin
%is the minimum balckground light intensity.
Imax=0.61;
Imin=0.09;
% mu is the interparticle friction coefficient, we set it to be much larger
% than the physical value to avoid artificial accumulation of solved forces
% to lie on the preset mu boundary.
mu=3;
%% The following line should NOT be here.
% printfig is used to print out the step-wise optimization process for
% demonstration process
printfig = 0;
%% Readin the data from a former step...
% it must be built on the unbalance solved results.
%particle_old=load([pwd '/Processing/Inverse/SolvedMat/' file_name_previous '_solved_3.mat'],'pres');
% if exist(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce.txt'],'file') == 2
%     Mat = csvread(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce.txt']);
% else
%     Mat = csvread(['Processing/Inverse/SolvedMat/' file_name '_solved_3.txt']);
% end
%
totally_refit = 0;

if  exist(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce_update_final.txt'],'file') == 2
    disp('use current update final file');
    Mat = csvread(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce_update_final.txt']);
     read_name = [file_name '_balance_reduce_update_final.txt'];
elseif  exist(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce_final.txt'],'file') == 2
    Mat = csvread(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce_final.txt']);
    read_name = [file_name '_balance_reduce_final.txt'];
elseif exist(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce_manual.txt'],'file') == 2
    disp('use current balance reduce manual');
    Mat = csvread(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce_manual.txt']);
    read_name = [file_name '_balance_reduce_manual.txt'];
elseif exist(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce_update.txt'],'file') == 2
    Mat = csvread(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce_update.txt']);
    read_name = [file_name '_balance_reduce_update.txt'];
elseif exist(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce.txt'],'file') == 2
    Mat = csvread(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce.txt']);
    read_name = [file_name '_balance_reduce.txt'];
elseif exist(['Processing/Inverse/SolvedMat/' file_name '_solved_3.txt'],'file') == 2
    Mat = csvread(['Processing/Inverse/SolvedMat/' file_name '_solved_3.txt']);
    read_name = [file_name '_solved_3.txt'];
else
    totally_refit = 1;
    Mat = csvread(['Processing/Inverse/SolvedMat/' file_name '_prepare.txt']);
    read_name = [file_name '_prepare.txt'];
end
disp(['current file is ' read_name]);
% Mat_old is the forces from the previous run

if  exist(['Processing/Inverse/SolvedMat/' file_name_previous '_balance_reduce_update_final.txt'],'file') == 2
    disp('use previous update final file');
    Mat_old = csvread(['Processing/Inverse/SolvedMat/' file_name_previous '_balance_reduce_update_final.txt']);
    read_previous_name = [file_name_previous '_balance_reduce_update_final.txt'];
elseif  exist(['Processing/Inverse/SolvedMat/' file_name_previous '_balance_reduce_final.txt'],'file') == 2
    Mat_old = csvread(['Processing/Inverse/SolvedMat/' file_name_previous '_balance_reduce_final.txt']);
    read_previous_name = [file_name_previous '_balance_reduce_final.txt'];
elseif exist(['Processing/Inverse/SolvedMat/' file_name_previous '_balance_reduce_manual.txt'],'file') == 2
    disp('use previous manual');
    read_previous_name = [file_name_previous '_balance_reduce_manual.txt'];
    Mat_old = csvread(['Processing/Inverse/SolvedMat/' file_name_previous '_balance_reduce_manual.txt']);
elseif exist(['Processing/Inverse/SolvedMat/' file_name_previous '_balance_reduce_update.txt'],'file') == 2
    read_previous_name = [file_name_previous '_balance_reduce_update.txt'];
    Mat_old = csvread(['Processing/Inverse/SolvedMat/' file_name_previous '_balance_reduce_update.txt']);
elseif exist(['Processing/Inverse/SolvedMat/' file_name_previous '_balance_reduce.txt'],'file') == 2
    read_previous_name = [file_name_previous '_balance_reduce.txt'];
    Mat_old = csvread(['Processing/Inverse/SolvedMat/' file_name_previous '_balance_reduce.txt']);
elseif exist(['Processing/Inverse/SolvedMat/' file_name_previous '_solved_3_balance_reduce.txt'],'file') == 2
    read_previous_name = [file_name_previous '_solved_3_balance_reduce.txt'];
    Mat_old = csvread(['Processing/Inverse/SolvedMat/' file_name_previous '_solved_3_balance_reduce.txt']);
    disp('use previous 3_balance_reduce');
elseif exist(['Processing/Inverse/SolvedMat/' file_name_previous '_solved_3.txt'],'file') == 2
    read_previous_name = [file_name_previous '_solved_3.txt'];
    Mat_old = csvread(['Processing/Inverse/SolvedMat/' file_name_previous '_solved_3.txt']);
else
    read_previous_name = [file_name_previous '_prepare.txt'];
    Mat_old = csvread(['Processing/Inverse/SolvedMat/' file_name_previous '_prepare.txt']);
end
disp(['previous file is ' read_previous_name]);
% This is the current image.
original_img = imread(['Processing/Adj_PL_3/' file_name '_3_PLonly.png']);
original_img = im2double(original_img);
% original_old = imread(['Processing/Adj_PL_3/' file_name_previous '_3_PLonly.png']);
% original_old = im2double(original_old);

% Read out the stored information
% N is the number of particles.
N = length(Mat)/5;
% ratio between physical length and the pixel length.
meterperpixel = Mat(1,6) / Mat(1,5);
% read the stress-optic constant
fsigma = Mat(1,7);
% R_b is the pixel radius for the large particle.
R_b = max(Mat(1:N,5));
R_s = min(Mat(1:N,5));

% calculate c_mask that would be used later
% Create two square masks for big and small particles respectively.
[x,y]=meshgrid(-R_b:R_b,-R_b:R_b);
% D records the distance to center from each pixels
Dis=sqrt(x.^2+y.^2);
% Ds and Db are masks for small and big particles.
Ds=Dis;
Ds(Dis>R_s)=0;
Ds(Dis<=R_s)=1;
c_mask_small=Ds;
Db=Dis;
Db(Dis>R_b)=0;
Db(Dis<=R_b)=1;
c_mask_big=Db;


%% Create a list to track the particles, we assume a particle do not move longer than a small radius in one step
N_old = length(Mat_old)/5;
P_old = [Mat_old(1:N_old,2),Mat_old(1:N_old,3)];
P = [Mat(1:N,2) Mat(1:N,3)];
% The ith particle is the old_index(i) th particle in P_old.
old_index = zeros(1,N);
for k = 1:N
    Dis = sqrt((P(k,1) - P_old(:,1)).^2+(P(k,2)-P_old(:,2)).^2);
    temp = find(Dis==min(Dis));
    if Dis(temp(1))<R_s
        old_index(k) = temp(1);
    end
end


%% Find current bad particles


ForceImgs = cell(1,N);
Original_Errors = zeros(1,N);
Original_SynthImgs = cell(1,N);

num_bad = 0 ;
if totally_refit == 1
    refit_index = ones(1,N);
else
    refit_index = zeros(1,N);
end

previous_manual = zeros(1,N);

for n = 1:N % Note N is the number of particles at this step, in case some were missed by the center finding algorithm.
    % calculate the original image, namely forceImage.
    x = Mat(n,2); y = Mat(n,3);
    cropXstart=round(x-R_b);
    cropXstop=round(x-R_b)+2*R_b;
    cropYstart=round(y-R_b);
    cropYstop=round(y-R_b)+2*R_b;
    forceImage = original_img(cropYstart:cropYstop,cropXstart:cropXstop);
    ForceImgs{n} = forceImage;
    
    z = Mat(n,8);
    % calculate the fit image, namely synthImage
    betas = Mat(3*N+n,2:z+1);
    pixel_r = Mat(n,5);
    if pixel_r == R_b
        c_mask = c_mask_big;
    else
        c_mask = c_mask_small;
    end
    forces = Mat(N+n,2:z+1);
    alphas = Mat(2*N+n,2:z+1);
    synthImg = Force_Image_new(betas,forces,alphas,size(forceImage,1),pixel_r,Imax,Imin,meterperpixel,fsigma);
    Original_SynthImgs{n} = synthImg;
    temp_diff = ( (synthImg - forceImage).*c_mask ).^2;
    original_fit_error = sqrt(sum(sum(temp_diff))/sum(sum(c_mask)));
    Original_Errors(n) = original_fit_error;
    
    % find the differences between reaction-action forces. % here we
    % consider all reaction forces, not only the good ones.
    [reaction_forces,~,original_forces,~] = get_reaction_forces(Mat,n,N,1,Original_Errors);
    
    delta_f = abs(reaction_forces - original_forces);
    if ~isempty(delta_f)
        if ( (original_fit_error>Err_thresh) || (max(delta_f)>F_thresh) ) && Mat(n,2)<5100 && Mat(n,9) == 1
            num_bad = num_bad + 1;
            refit_index(n) = 1;
        end
    end
end

% track errors
current_fit_error = Original_Errors;

%% Calculate the old force images


%
% ForceImgs_old = cell(1,N_old);
%
% for n = 1:N_old % Note N is the number of particles at this step, in case some were missed by the center finding algorithm.
%     % calculate the original image, namely forceImage.
%     if Mat_old(n,9)==1
%     x = Mat_old(n,2); y = Mat_old(n,3);
%     cropXstart=round(x-R_b);
%     cropXstop=round(x-R_b)+2*R_b;
%     cropYstart=round(y-R_b);
%     cropYstop=round(y-R_b)+2*R_b;
%     %disp(cropXstart);
%     %disp(cropXstop);
%     forceImage_old = original_old(cropYstart:cropYstop,cropXstart:cropXstop);
%     ForceImgs_old{n} = forceImage_old;
%     end
% end


%%


tic;
% remove_contact function removes all the contacts with F<F_th from pres
% structure. However, after such process, we may want to re-fit the
% particles that was changed. But at this point, we want to refit all the
% particles, in order to improve accuracy...

% I wish to remove a contact only when forces below F_th on both sides of
% the contact!


% set whether or not to display the details of the fitting results.
if display_setting == 1
    opt_display_setting = 'final-detailed';
else
    opt_display_setting = 'off';
end


%% Use previous forces as initial guesses
% Construct a C matrix for the convinience of future calculations.
% Some parameters for display.
C = zeros(N,2);
for i=1:N
    % C(i,1) is z.
    C(i,1) = Mat(i,8);
    % C(i,2) is the boundary index.
    C(i,2) = Mat(i,9);
end
% % Maximum number of contact
zmax=max(C(:,1));
num_particle_with_z = zeros(1,zmax);
CC = C(C(:,2)==1,1);
for z = 1:zmax
    num_particle_with_z(z) = sum(CC == z);
end



%% The actual fitting part.
for nz=2:zmax % This is the wave fitting method, fit particles with same contact number first. We do not fit z = 1 particle.
    if nz>=6
        fitoptions = optimoptions('lsqnonlin','MaxIter',6000*iteration_ratio,'MaxFunEvals',6000*iteration_ratio,'TolFun',10^(-15),'Display',opt_display_setting);
    elseif nz == 5
        fitoptions = optimoptions('lsqnonlin','MaxIter',5000*iteration_ratio,'MaxFunEvals',5000*iteration_ratio,'TolFun',10^(-15),'Display',opt_display_setting);
    elseif nz == 4
        fitoptions = optimoptions('lsqnonlin','MaxIter',4000*iteration_ratio,'MaxFunEvals',4000*iteration_ratio,'TolFun',10^(-15),'Display',opt_display_setting);
    elseif nz <=3
        fitoptions = optimoptions('lsqnonlin','MaxIter',3000*iteration_ratio,'MaxFunEvals',3000*iteration_ratio,'TolFun',10^(-15),'Display',opt_display_setting);
    end
    z_particle_index = 0; % z_particle_index is a index to record the number of particles with nz contacts we are processing
    for n=1:N % Loop each particle to perform the fit.
        z = Mat(n,8);
        
        if z==nz  && Mat(n,9) == 1 && refit_index(n) == 1 %&& n ==1030  %% IMPORTANT: for particle need to be refitted.
            verbose_temp = verbose(z);
            %Display the number of particles we are processing
            z_particle_index = z_particle_index + 1;
            if display_setting == 1
                disp('---------------------------------------------------------------------------------------');
            end
            %% Here for each particle, we fit it multiple time, until such a particle's stress does not change much.
            display(['id ' num2str(n) ' Old enhance: Current Processing: ' num2str(z) ' contact fittings: fitting force(s) to particle ',num2str(z_particle_index) ' / ' num2str(num_particle_with_z(z)) ]);
            %% make sure to check whether this is correct.
            % calculate some basic information
            [forces,alphas,betas,pixel_r,c_mask,forceImage] = read_force_information(Mat,n,N,R_b,ForceImgs,c_mask_small,c_mask_big);
            synthImg = Force_Image_new(betas,forces,alphas,size(forceImage,1),pixel_r,Imax,Imin,meterperpixel,fsigma); % reconstructed image.
            temp_diff = ( (synthImg - forceImage).*c_mask ).^2;
            original_fit_error = sqrt(sum(sum(temp_diff))/sum(sum(c_mask))); % original error is the fit error before refit.
            
            
            %
            % Use old forces as initial guesses.
            [old_forces,old_alphas,old_manual] = find_previous_force(Mat,n,Mat_old,old_index);
            go_on = 1;
            if update_based_on_manual == 1 && old_manual == 0
                go_on = 0;
            end
            if go_on == 1
            
            %             figure,imagesc(ForceImgs{n});
            %             title('new');
            %             % Plot contacts on old...
            %             axis image;
            %             hold on;
            %             xcenter = 1/2*length(ForceImgs{n});
            %             ycenter = xcenter;
            %             for k = 1:Mat(n,8)
            %                 beta_c = Mat(3*N+n,k+1);
            %                 %alpha_c = Mat(2*N+n,k+1);
            %                 alpha_c = old_alphas(k);
            %                 xc = xcenter + pixel_r*cos(beta_c);
            %                 yc = ycenter - pixel_r*sin(beta_c);
            %                 plot([xc xc+10*cos(alpha_c)],[yc yc-10*sin(alpha_c)],'LineWidth',2,'Color','k');
            %             end
            %
            %
            %             figure,imagesc(ForceImgs_old{old_index(n)});
            %             title('old');
            %             axis image;
            %             hold on;
            %             n_old = old_index(n);
            %             xcenter = 1/2*length(ForceImgs{n});
            %             ycenter = xcenter;
            %             for k = 1:Mat_old(n_old,8)
            %                 beta_c = Mat_old(3*N_old+n_old,k+1);
            %                 alpha_c = Mat_old(2*N_old+n_old,k+1);
            %                 xc = xcenter + pixel_r*cos(beta_c);
            %                 yc = ycenter - pixel_r*sin(beta_c);
            %                 plot([xc xc+10*cos(alpha_c)],[yc yc-10*sin(alpha_c)],'LineWidth',2,'Color','k');
            %             end
            
            
            if old_manual~=0
                %[forces,alphas] = get_corresponding_balance(betas,old_alphas,old_forces,mu); % forces and alphas are balanced solutions based on unbalanced results.
                %if z >= 3
                %    [forces,alphas] = get_corresponding_balance(betas,old_alphas,old_forces,mu); % forces and alphas are balanced solutions based on unbalanced results.
                %else
                forces = old_forces;
                alphas = old_alphas; % forces and alphas are balanced solutions based on unbalanced results.
                %end
                Mat(n,10) = 0.5;
                
                disp('initial gauess forces =');
                disp(forces);
                disp('initial gauess alphas =');
                disp(alphas/pi*180);
                
                verbose_nobalance = 0;
                [synthImg,forces,alphas] = fitforce_nobalance_new(forceImage,size(forceImage,1),betas,pixel_r,forces,alphas,verbose_nobalance,z,c_mask,printfig,Imax,Imin,meterperpixel,fsigma,mu,fitoptions);
                      %  input('input anything here = ');
                      %  close();
                      % figure;imagesc(nobalance_Img,[0,0.6]); axis image;
                      % input('input anything here = ');
                      % close();
                
            else
                verbose_nobalance = 0;
                [~,forces_nobalance,alphas_nobalance] = fitforce_nobalance_new(forceImage,size(forceImage,1),betas,pixel_r,old_forces,old_alphas,verbose_nobalance,z,c_mask,printfig,Imax,Imin,meterperpixel,fsigma,mu,fitoptions);
                if z>=3
                    [forces,alphas] = get_corresponding_balance(betas,alphas_nobalance,forces_nobalance,mu); % forces and alphas are balanced solutions based on unbalanced results.
                else
                    forces = forces_nobalance;
                    alphas = alphas_nobalance;
                end
            end
            
            % The following code performs the balance fit part.
            
            if update_based_on_manual == 0
                [synthImg,forces,alphas] = fitforce_balance_new(forceImage,size(forceImage,1),betas,pixel_r,forces,alphas,verbose_temp,z,c_mask,printfig,Imax,Imin,meterperpixel,fsigma,mu,fitoptions);
            end
            disp('fitted forces =');
            disp(forces);
            disp('fitted alphas =');
            disp(alphas/pi*180);
            
            
            
            % calculate the new error
            temp_diff = ( (synthImg - forceImage).*c_mask ).^2;
            new_fit_error = sqrt(sum(sum(temp_diff))/sum(sum(c_mask)));
            
            if (new_fit_error < original_fit_error) || old_manual ~=0
                Mat(N+n,2:z+1) = forces;
                Mat(2*N+n,2:z+1) = alphas;
                current_fit_error(n) = new_fit_error;
                if old_manual == 0
                    disp(['old forces improve fitting quality -- old error = ' num2str(original_fit_error,4) ' ; new error = ' num2str(new_fit_error,4)]);
                else
                    disp(['record data based on old manual -- old error = ' num2str(original_fit_error,4) ' ; new error = ' num2str(new_fit_error,4)]);
                end
            end
            end
        end
    end
end



%% Update using the reaction forces.

if update_based_on_manual ~= 1
    
    for nz=2:zmax % This is the wave fitting method, fit particles with same contact number first. We do not fit z = 1 particle.
        if nz>=6
            fitoptions = optimoptions('lsqnonlin','MaxIter',6000*iteration_ratio,'MaxFunEvals',6000*iteration_ratio,'TolFun',10^(-15),'Display',opt_display_setting);
        elseif nz == 5
            fitoptions = optimoptions('lsqnonlin','MaxIter',5000*iteration_ratio,'MaxFunEvals',5000*iteration_ratio,'TolFun',10^(-15),'Display',opt_display_setting);
        elseif nz == 4
            fitoptions = optimoptions('lsqnonlin','MaxIter',4000*iteration_ratio,'MaxFunEvals',4000*iteration_ratio,'TolFun',10^(-15),'Display',opt_display_setting);
        elseif nz <=3
            fitoptions = optimoptions('lsqnonlin','MaxIter',3000*iteration_ratio,'MaxFunEvals',3000*iteration_ratio,'TolFun',10^(-15),'Display',opt_display_setting);
        end
        z_particle_index = 0; % z_particle_index is a index to record the number of particles with nz contacts we are processing
        for n=1:N % Loop each particle to perform the fit.
            z = Mat(n,8);
            
            % here we only consider reaction forces came from good neighbours.
            [reaction_forces,reaction_alphas,original_forces,~] = get_reaction_forces(Mat,n,N,0,current_fit_error);
            delta_f = abs(reaction_forces - original_forces);
            
            if z==nz  && Mat(n,9) == 1 && max(delta_f)>F_thresh %% IMPORTANT: for particle need to be refitted.
                verbose_temp = verbose(z);
                %Display the number of particles we are processing
                z_particle_index = z_particle_index + 1;
                if display_setting == 1
                    disp('---------------------------------------------------------------------------------------');
                end
                %% Here for each particle, we fit it multiple time, until such a particle's stress does not change much.
                display(['Reaction force enhance: Current Processing: ' num2str(z) ' contact fittings: fitting force(s) to particle ',num2str(z_particle_index) ' / ' num2str(num_particle_with_z(z))]);
                %% make sure to check whether this is correct.
                % calculate some basic information
                [forces,alphas,betas,pixel_r,c_mask,forceImage] = read_force_information(Mat,n,N,R_b,ForceImgs,c_mask_small,c_mask_big);
                synthImg = Force_Image_new(betas,forces,alphas,size(forceImage,1),pixel_r,Imax,Imin,meterperpixel,fsigma); % reconstructed image.
                temp_diff = ( (synthImg - forceImage).*c_mask ).^2;
                original_fit_error = sqrt(sum(sum(temp_diff))/sum(sum(c_mask))); % original error is the fit error before refit.
                
                % Use reaction forces as initial guesses.
                verbose_nobalance = 0;
                [~,forces_nobalance,alphas_nobalance] = fitforce_nobalance_new(forceImage,size(forceImage,1),betas,pixel_r,reaction_forces,reaction_alphas,verbose_nobalance,z,c_mask,printfig,Imax,Imin,meterperpixel,fsigma,mu,fitoptions);
                if z>=3
                    [forces,alphas] = get_corresponding_balance(betas,alphas_nobalance,forces_nobalance,mu); % forces and alphas are balanced solutions based on unbalanced results.
                else
                    forces = forces_nobalance;
                    alphas = alphas_nobalance;
                end
                
                % The following code performs the balance fit part.
                [synthImg,forces,alphas] = fitforce_balance_new(forceImage,size(forceImage,1),betas,pixel_r,forces,alphas,verbose_temp,z,c_mask,printfig,Imax,Imin,meterperpixel,fsigma,mu,fitoptions);
                % calculate the new error
                temp_diff = ( (synthImg - forceImage).*c_mask ).^2;
                new_fit_error = sqrt(sum(sum(temp_diff))/sum(sum(c_mask)));
                
                if new_fit_error < original_fit_error
                    Mat(N+n,2:z+1) = forces;
                    Mat(2*N+n,2:z+1) = alphas;
                    current_fit_error(n) = new_fit_error;
                    disp('reaction forces improve fitting quality');
                end
            end
        end
    end
    
    if iteration_ratio == 1
        %csvwrite(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce_update_final.txt'],Mat);
    else
        %csvwrite(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce_update.txt'],Mat);
    end
    %% Remove the small forces.
    [Mat, num_weak_contact,List_reduced] = remove_contact(Mat,F_th);
    
    % Make sure to update the fitting for all particles at least once.
    if totally_refit ~= 1
        List_reduced = ones(1,N);
    end
    while sum(List_reduced)>0
        % Construct a C matrix for the convinience of future calculations.
        %% The actual fitting part.
        for nz=2:zmax % This is the wave fitting method, fit particles with same contact number first. We do not fit z = 1 particle.
            if nz>=6
                fitoptions = optimoptions('lsqnonlin','MaxIter',6000*iteration_ratio,'MaxFunEvals',6000*iteration_ratio,'TolFun',10^(-15),'Display',opt_display_setting);
            elseif nz == 5
                fitoptions = optimoptions('lsqnonlin','MaxIter',5000*iteration_ratio,'MaxFunEvals',5000*iteration_ratio,'TolFun',10^(-15),'Display',opt_display_setting);
            elseif nz == 4
                fitoptions = optimoptions('lsqnonlin','MaxIter',4000*iteration_ratio,'MaxFunEvals',4000*iteration_ratio,'TolFun',10^(-15),'Display',opt_display_setting);
            elseif nz <=3
                fitoptions = optimoptions('lsqnonlin','MaxIter',3000*iteration_ratio,'MaxFunEvals',3000*iteration_ratio,'TolFun',10^(-15),'Display',opt_display_setting);
            end
            % z_particle_index is a index to record the number of particles with nz contacts we
            % are processing
            z_particle_index = 0;
            for n=1:N % Loop each particle to perform the fit.
                z = Mat(n,8);
                if z==nz  && C(n,2) == 1 && List_reduced(n) == 1 % && n==500%Only fit particles that are with particular number of contacts (wave fitting). And we must make sure the particle is NOT a boundary particle C(n,2)=1
                    %Display the number of particles we are processing
                    z_particle_index = z_particle_index + 1;
                    if display_setting == 1
                        disp('---------------------------------------------------------------------------------------');
                    end
                    %% Here for each particle, we fit it multiple time, until such a particle's stress does not change much.
                    display(['Id '  num2str(n) ' : Rmove small forces: Current Processing: ' num2str(z) ' contact fittings: fitting force(s) to particle ',num2str(z_particle_index) ' / ' num2str(num_particle_with_z(z)) ]);
                    % read force information
                    [forces,alphas,betas,pixel_r,c_mask,forceImage] = read_force_information(Mat,n,N,R_b,ForceImgs,c_mask_small,c_mask_big);
                    synthImg = Force_Image_new(betas,forces,alphas,size(forceImage,1),pixel_r,Imax,Imin,meterperpixel,fsigma);
                    temp_diff = ( (synthImg - forceImage).*c_mask ).^2;
                    original_fit_error = sqrt(sum(sum(temp_diff))/sum(sum(c_mask)));
                    if z>=3
                        verbose_nobalance = 0;
                        [~,forces_nobalance,alphas_nobalance] = fitforce_nobalance_new(forceImage,size(forceImage,1),betas,pixel_r,forces,alphas,verbose_nobalance,z,c_mask,printfig,Imax,Imin,meterperpixel,fsigma,mu,fitoptions);
                        [forces,alphas] = get_corresponding_balance(betas,alphas_nobalance,forces_nobalance,mu);
                    end
                    verbose_temp = verbose(z);
                    % The following code performs the balance fit part.
                    [synthImg,forces,alphas] = fitforce_balance_new(forceImage,size(forceImage,1),betas,pixel_r,forces,alphas,verbose_temp,z,c_mask,printfig,Imax,Imin,meterperpixel,fsigma,mu,fitoptions);
                    % calculate the new error
                    temp_diff = ( (synthImg - forceImage).*c_mask ).^2;
                    new_fit_error = sqrt(sum(sum(temp_diff))/sum(sum(c_mask)));
                    %% Here we evaluate the change of the forces as well as the change of angles after such an optimization. If the forces and the angles changed a lot, we want to do another fitting to make sure the optimization was achieved.
                    if new_fit_error < original_fit_error
                        % assign new forces and alphas
                        Mat(N+n,2:z+1) = forces;
                        Mat(2*N+n,2:z+1) = alphas;
                        current_fit_error(n) = new_fit_error;
                        disp('remove forces improve fitting quality');
                    end
                end
            end
        end
        [Mat,num_weak_contact,List_reduced] = remove_contact(Mat,F_th);
    end
    
end

toc;


%% Save the outputs.


if iteration_ratio == 1
    if update_based_on_manual == 1
        csvwrite(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce_update_final_fixed.txt'],Mat);
    else
        csvwrite(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce_update_final.txt'],Mat);    
    end
else
    csvwrite(['Processing/Inverse/SolvedMat/' file_name '_balance_reduce_update.txt'],Mat);
end

%% Helper 1: remove_contact
function [Mat,num_weak_contact,List_reduced] = remove_contact(Mat,F_th)
% This supporter function is used to remove all the contacts with force
% smaller than F_th. Note, a contact is removed, only when the fitted force
% magnitude on both of the contacting particles is below F_th.



num_weak_contact = 0;
N = length(Mat)/5;
% if List_reduced(i) = 1 the ith particle experience a contact change, and will be used to fit again in the
% next while loop.
List_reduced = zeros(1,N);
for n = 1:N
    if Mat(n,9) == 1 % this is a non-boundary contact.
        z = Mat(n,8);
        if z>= 2
            rm_index = zeros(1,z);
            for k = 1:z
                current_force = Mat(N+n,k+1);
                if current_force < F_th
                    % Now on this side the contact is a weak contact, how
                    % about the other side?
                    %% Find the neighbor particle and check if the reaction force is larger than the threshold
                    % temp_index is the index of the other particle.
                    temp_index = Mat(4*N+n,k+1);
                    other_z = Mat(temp_index,8);
                    rm_index_temp = zeros(1,other_z);
                    % we note the other particle could have only 1 contact,
                    % in which case we defenitely remove the contact. But,
                    % if the other particle has more than 1 contact, we
                    % need to consider the fitted force magnitude on that
                    % particle.
                    % initialization to a checker.
                    reaction_force_index = 0;
                    other_boundary = Mat(temp_index,9);
                    if other_z <= 1 || other_boundary == 0
                        reaction_force_index = 1;
                    else
                        other_neighbours = Mat(4*N+temp_index,2:other_z+1);
                        other_forces = Mat(N+temp_index,2:other_z+1);
                        for j = 1 : other_z
                            if other_neighbours(j) == n && other_forces(j) <F_th
                                % the reaction force is also small, need to
                                % remove.
                                reaction_force_index = 1;
                                rm_index_temp(j) = 1;
                            end
                        end
                    end
                    
                    % if the reaction force is small as well.
                    if reaction_force_index == 1
                        List_reduced(n) = 1;
                        List_reduced(temp_index) = 1;
                        % NOW we ACTUALLY need to remove weak contact.
                        num_weak_contact = num_weak_contact + 1;
                        % remove this contact only under such case.
                        rm_index(k) = 1;
                        % This part we also need to remove the corresponding
                        % particle's information.
                        
                        other_z_new = sum(~rm_index_temp);
                        % update neighbours
                        other_neighbours = Mat(4*N+temp_index,2:other_z+1);
                        Mat(4*N+temp_index,2:other_z+1) = 0;
                        other_neighbours = other_neighbours(~rm_index_temp);
                        Mat(4*N+temp_index,2:other_z_new+1) = other_neighbours;
                        
                        % update betas
                        other_betas = Mat(3*N+temp_index,2:other_z+1);
                        Mat(3*N+temp_index,2:other_z+1) = 0;
                        other_betas = other_betas(~rm_index_temp);
                        Mat(3*N+temp_index,2:other_z_new+1) = other_betas;
                        
                        % update other forces and alphas
                        other_forces = Mat(N+temp_index,2:other_z+1);
                        other_alphas = Mat(2*N+temp_index,2:other_z+1);
                        if sum(other_forces)>0
                            % update other forces
                            Mat(N+temp_index,2:other_z+1) = 0;
                            other_forces = other_forces(~rm_index_temp);
                            Mat(N+temp_index,2:other_z_new+1) = other_forces;
                            
                            % update other alphas
                            Mat(2*N+temp_index,2:other_z+1) = 0;
                            other_alphas = other_alphas(~rm_index_temp);
                            Mat(2*N+temp_index,2:other_z_new+1) = other_alphas;
                        end
                        % update the contact number.
                        Mat(temp_index,8) = other_z_new;
                    end
                end
            end
            
            %% remove the current contact
            
            % update current forces
            z_new = sum(~rm_index);
            current_forces = Mat(N+n,2:z+1);
            Mat(N+n,2:z+1) = 0;
            current_forces = current_forces(~rm_index);
            Mat(N+n,2:z_new+1) = current_forces;
            
            % update current neighbours
            current_neighbours = Mat(4*N+n,2:z+1);
            Mat(4*N+n,2:z+1) = 0;
            current_neighbours = current_neighbours(~rm_index);
            Mat(4*N+n,2:z_new+1) = current_neighbours;
            
            % update current betas
            current_betas = Mat(3*N+n,2:z+1);
            Mat(3*N+n,2:z+1) = 0;
            current_betas = current_betas(~rm_index);
            Mat(3*N+n,2:z_new+1) = current_betas;
            
            % update current alphas
            current_alphas = Mat(2*N+n,2:z+1);
            Mat(2*N+n,2:z+1) = 0;
            current_alphas = current_alphas(~rm_index);
            Mat(2*N+n,2:z_new+1) = current_alphas;
            
            % update the current z
            Mat(n,8) = z_new;
            
        end
    end
end


disp(['Reduced ' num2str(num_weak_contact) ' number of contacts']);

%% Helper function 2: get_reaction_forces

function [reaction_forces,reaction_alphas,original_forces,original_alphas] = get_reaction_forces(Mat,n,N,type,current_fit_error)
% Mat is the overal Mat, n is the index of the particle
% out put would be the magnitude and alphas given by the reaction forces.
% Note, if a reaction force is from a non-fitted particle then the original
% force and alpha are used in the output array.

% if type is 1 : get all reaction forces
% if type is 0 : get only good reaction forces.

z = Mat(n,8);

original_forces = Mat(N+n,2:z+1);
original_alphas = Mat(2*N+n,2:z+1);
reaction_forces = original_forces;
reaction_alphas = original_alphas;

if Mat(n,9) == 1 && z>=2 % this is a non-boundary contact and has enough number of contacts.
    
    for k = 1:z
        %% Find the neighbor particle and check if the reaction force is larger than the threshold
        % temp_index is the index of the other particle.
        temp_index = Mat(4*N+n,k+1);
        other_z = Mat(temp_index,8);
        other_boundary = Mat(temp_index,9);
        if other_z <= 1 || other_boundary == 0
            % this is not a fitted particle % do nothing.
        else
            other_neighbours = Mat(4*N+temp_index,2:other_z+1);
            other_forces = Mat(N+temp_index,2:other_z+1);
            other_alphas = Mat(2*N+temp_index,2:other_z+1);
            for j = 1 : other_z
                if other_neighbours(j) == n %&& Original_Errors(temp_index) <0.15 % this is the reaction forces. The error threshold is chosen empirically.
                    if (type == 1) || (type==0 && (current_fit_error(temp_index) < current_fit_error(n)))
                        reaction_forces(k) = other_forces(j);
                        reaction_alphas(k) = mod(other_alphas(j)+pi,2*pi);
                    end
                end
            end
        end
    end
end

%% Helper function 3: find previous forces
function [old_forces,old_alphas,old_manual] = find_previous_force(Mat,n,Mat_old,old_index)
% Mat is the information matrix
% n is the particle index

N = length(Mat)/5;
N_old = length(Mat_old)/5;
z = Mat(n,8);
old_forces = Mat(N+n,2:z+1);
old_alphas = Mat(2*N+n,2:z+1);

old_n = old_index(n);
old_manual = Mat_old(old_n,10);

for k=1:z
    neighbour = Mat(4*N+n,k+1);
    % find previous n and previous neighbour
    old_neighbour = old_index(neighbour);
    
    
    % find if previous_neighbour and previous_n has a contact
    % disp(4*N_old+old_n);
    if old_n~=0
        old_n_neighbourlist = Mat_old(4*N_old+old_n,2:Mat_old(old_n,8)+1);
        temp = find(old_n_neighbourlist==old_neighbour);
        if isempty(temp)
            % if there is no old contact
            old_forces(k) = 0.005;
            old_alphas(k) = mod(Mat(3*N+n,k+1) + pi , 2*pi);
        else
            % there is a old contact
            if Mat_old(old_n,9) == 1
                old_forces(k) = Mat_old(N_old + old_n , temp+1);
                old_alphas(k) = Mat_old(2*N_old + old_n , temp+1);
            end
        end
    end
end

% calculate the angle between old_alphas and old_betas
%  angles = zeros(1,z);
%  new_betas = Mat(3*N+n,2:z+1);
%  for k = 1:z
%     vec_r = [cos(new_betas(k)),-sin(new_betas(k))];
%     vec_old_f = [cos(old_alphas(k)),-sin(old_alphas(k))];
%     cosphi = sum(vec_r.*vec_old_f);
%     angles(k) = acos(cosphi);
%  end
%  disp(angles/pi*180);




%% helper function, read information from Mat file for a given particle.
function [forces,alphas,betas,pixel_r,c_mask,forceImage] = read_force_information(Mat,n,N,R_b,ForceImgs,c_mask_small,c_mask_big)
% this function read the force information needed for fitting from the Mat
% matrix.
z = Mat(n,8);
forceImage = ForceImgs{n};
% calculate the fit image, namely synthImage
betas = Mat(3*N+n,2:z+1);
pixel_r = Mat(n,5);
forces = Mat(N+n,2:z+1);
alphas = Mat(2*N+n,2:z+1);
if pixel_r == R_b
    c_mask = c_mask_big;
else
    c_mask = c_mask_small;
end






