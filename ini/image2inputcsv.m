function image2inputcsv(imagefile, csvname)
    %IMAGE2INPUTCSV convert an image to a csv as an input for our final project
    %   Turns a image (.png or .jpeg file) into a .csv file
    %   where 0 denotes a fluid cell and 1 an obstacle cell. The program
    %   checks for the two cell criteria and adjusts the obstacle if
    %   needed. 
    %   The input image can either be an rgb or greyscale image. Obstacles
    %   are assumbed to be pixels which greyvalues are <= 256/2. 
    % 
    %   Warning! This does not check the two cell criteria at the boundary
    %   cells as this depends on the boundary conditions! Set them
    %   accordingly in the .txt files. 
    % 
    %   Input: 
    %       - imagefile (string): e.g. "test-case.png"
    %       - csvname (string): name for output csv file e.g. "test-case"
    %
    %   by Kim Kröner, Magnus Ostertag, Janis Wissinger
    %   Last modified 8.02.2021
    %   MATLAB R2021a  
    
    
    %% read in image as matrix
    A = imread(imagefile);
    
    % turn it into greyscale 
    A = im2gray(A);
    
    % imshow(A) % as a sanity check
    
    % convert it so that a 1 indicates an obstacle cell and 0 a fluid cell
    A_obs = 1 - round(A./256);
    
    % imshow(A_obs*256);
    % heatmap(A_obs) % as an alternative santity check
    
    %% check two cells criterium for obstacle cells
    
    % iterate until no cells are changed to fluid anymore
    modifiedCells = true;
    while modifiedCells
        
        modifiedCells = false; % init for this iteration
        
        % go over the inner domain, outer domain
        for i=2:size(A_obs, 1)-1
            for j=2:size(A_obs, 2)-1
                
                if A_obs(i,j) == 1 % is obstacle cell
                    
                    % if two opposing sides of an obstacle cell are fluid, turn
                    % obstacle in fluid cell
                    if A_obs(i-1,j) + A_obs(i+1, j) == 0
                        A_obs(i,j) = 0;
                        modifiedCells = true;
                    end
                    if A_obs(i,j-1) + A_obs(i,j+1) == 0
                        A_obs(i,j) = 0;
                        modifiedCells = true;
                    end
                    
                end
            end
        end
    end
    
    %     heatmap(A_obs);
    
    %% write to csv
    name = strcat(csvname, ".csv");
    writematrix(A_obs, name)
end