%%
% Wasserstein Modulation Index (wMI)
%
% This function file ("Wasserstein_MI.m") calculates phase-amplitude coupling via the Wasserstein distance. 
% See  Ohki, 2022, J Neurosci Methods, DOI: 10.1016/j.jneumeth.2022.109578
% 
% How to use
% [wMI, optpath, MeanAmp] = Wasserstein_MI(phase, amplitude, position)
% NOTE: The current setting for H0 is defined as the uniform distribution (U).
%             I would recommend trying modifying this setting in a data-driven manner instead of using U.
% Inputs:
% Phase = phase time series
% Amp = amplitude time series
% position = phase bins
%
% Outputs:
% wMI = Wasserstein Modulation Index
% optpath = coupling phase matrix
% MeanAmp = amplitude distribution over phase bins (non-normalized)
%
% Questions -> takefumi2ohki@gmail.com


function [wMI, optpath, MeanAmp] = Wasserstein_MI(phase, amplitude, position)
       
    nbin = length(position);  
    winsize = 2*pi/nbin;

    % DataBox 
    MeanAmp=zeros(1, nbin);
    
    % now we compute the mean amplitude in each phase:
    for bin_idx = 1 : nbin
          I = find(phase <  position(bin_idx)+winsize & phase >=  position(bin_idx));
         MeanAmp(bin_idx)=mean(amplitude(I));
    end


    % transform amplitude values in each bin into probability
    MeanAmp = MeanAmp/sum(MeanAmp);

    % the linear equality constraint1
    Aeq = zeros(2*nbin, nbin^2);
    
    for bin_idx = 1 : nbin
         Aeq(bin_idx, (bin_idx-1) * nbin +1 : bin_idx * nbin) = 1;
         Aeq(nbin+bin_idx, bin_idx : nbin : end) = 1;
    end
    
    % the linear equality constraints: H0 and vectroize Mean amplitude Distribution
    Beq = [(ones(nbin,1)/nbin); MeanAmp']; % 1st dim = uniform dist via the built-in function 'one'; 2nd dim = empirial data
                                                                        % you can modify 1st term based on your data, your working hypothesis or your control data                
    
    %  the objective or cost function (i.e., the distance in the phase plane)
    c = zeros(nbin^2,1);
    
    for bin1_idx =1 : nbin
        for bin2_idx =1 : nbin
            
            %  L1 wasserstein distance
            c((bin1_idx-1)*nbin + bin2_idx) = min(abs(bin1_idx - bin2_idx), nbin-abs(bin1_idx- bin2_idx));
            
            % L2 wasserstein distance
            % c((bin1_idx-1)*nbin + bin2_idx) = min(abs(bin1_idx- bin2_idx), nbin-abs(bin1_idx - bin2_idx))^2;
        end
    end
    
    % the constrain term (i.e., the lower bound)
    lb = zeros(nbin^2, 1);
    
    % wMI calculation via linear programming; i.e., minimize c'*x SUCH THAT A*x = b, x >= 0
    % Set undisplaying a message: 'optimal solution found.' in order to
    % enhance the calculation speed.
    options = optimoptions('linprog','Display','none');
    
    % calculate the optimal trransport matrix and wMI
    [x, wMI] = linprog(c, [],[], Aeq, Beq, lb, [], options);

    % coupling phase matrix (phase bin * phase bin)
    optpath = reshape(x, nbin, length(x)/nbin);

     
end


