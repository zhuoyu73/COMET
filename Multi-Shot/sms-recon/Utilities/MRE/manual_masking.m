%% Manual Masking (Human Brain)
close all;
load('t2mask_bet.mat')
load('t2stack.mat')
tmp = (t2stack.*mask)./(max(t2stack)); % GE EPI Data or Siemens Spiral Data
%tmp = (t2stack.*mask)./(1000); %Siemens EPI Data
%for remask
%maskx = load('maskx_old.mat');
%maskx = maskx.maskx;

maskx = zeros(size(mask));

% look at the whole brain
figure;im(tmp);caxis([0 1])
while true
    % Display the figure to allow user to view the mask and slices
    figure;
    im(tmp .* mask .* abs(1 - maskx));
    caxis([0 1]);
    % Prompt user for the start and end slice
    start_slice = input('Enter the start slice: ');
    end_slice = input('Enter the end slice: ');
    disp(['Start slice: ', num2str(start_slice)]);
    disp(['End slice: ', num2str(end_slice)]);
    last_ss = start_slice;  % Initialize last_ss for error handling
    % Loop through slices in reverse order
    for ss = last_ss:-1:end_slice
        try
            disp(['Processing slice: ', num2str(ss)]);
            % Process the mask for the current slice
            maskx(:,:,ss) = double(roipoly(tmp(:,:,ss)));
            %maskx(:,:,ss) = double(roipoly(maskx(:,:,ss)));
        catch ME
            % Handle the error, restart the loop from the failed slice
            disp(['Error occurred at slice ', num2str(ss)]);
            disp('Restarting from the last successful slice...');
            % Display error message for debugging
            disp(ME.message);
            % Set last_ss to the current slice so that the loop restarts from here
            last_ss = ss;
            break;  % Break out of the for loop to restart from the error slice
        end
    end
    % Save the mask after processing all slices
    save('maskx.mat', 'maskx'); 
    % Display the final mask to check if it's correct
    figure;
    im(tmp .* mask .* abs(1 - maskx));
    caxis([0 1]);
    % Ask the user whether they are satisfied with the mask
    continue_input = input('Please enter Y if you are satisfied with your mask or N if you would like to make changes: ', 's');
    % If user is satisfied, break out of the while loop
    if strcmpi(continue_input, 'Y' )
        disp('User is satisfied with the mask.');
        break;
    % If user wants to make changes, restart the slice selection process
    elseif strcmpi(continue_input, 'N')
        disp('User wants to make changes. Re-enter slices.');
    % Handle invalid inputs
    else
        disp('Invalid input. Please enter Y or N.');
    end
end