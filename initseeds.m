function [segind, Labels, final_alpha_gcost, seeds] = initseeds(Image, L, f, conn, camera_noise, varargin)
%%init seed & hard constraint(INTERACTIVE GRAPH CUTS FOR OPTIMAL BOUNDARY & REGION SEGMENTATION...
%%OF OBJ IN N-D IMG bY YURI BOYKOV 2001)
%%change labeling based on alpha expansion move minimizes the current energy via... 
%%alpha-expansion until convergence(FAST APPROXIMATE ENERGY MINIMIZATION VIA GRAPH CUT by BOYKOV 2001)

%% Explanation of parameters:
%%Image: grayscaled and double-classed image to be segmented(if NOT we make it)
%%conn: 4-connected neighbourhood
%L: centroids of clusters in Kmeans clustering
%%f: current labels of all pixels created by Kmeans
%%camera_noise: estimated by mean of intensity diffrences between adjacenet pixels
%%prior_likelihood or varagin:custom prabability distribiution function fitted to a training sample data

I = Image;
if size(I, 3) == 3
    I = rgb2gray(I);%Image should be gray scale,if not,convert it
end
if ~isa(I, 'double')
    I = double(I);%Image should be double,if not,convert it
end
[n, m] = size(I);
I_data = I(:);
num_L = length(L);

prior_likelihood = {};
loaded_seeds = [];
switch nargin
    case 7
        prior_likelihood = varargin{1};
        loaded_seeds = varargin{2};
        if ~any(loaded_seeds)
            disp('no seeds exist for this slice in the init_seeds matrix you loaded')
            return
        end
    case 6
        if iscell(varargin{1})
            prior_likelihood = varargin{1};
        elseif ismatrix(varargin{1})
            loaded_seeds = varargin{1};
            if ~any(loaded_seeds)
                disp('no seeds exist for this slice in the init_seeds matrix you loaded')
                return
            end
        else
            disp('last input argumnet is not recognized as a valid input')
            return
        end
end

format compact
str = 'n';
if isempty(loaded_seeds)
    prompt = 'Do you want to mark seeds? (y/n) ';
    str = input(prompt, 's');
    if isempty(str)
        str = 'n';
    end
end


seed_dist = zeros(n*m, num_L); % density distribution of intensity of each segment, derived from the seeds of that segment
seeds = zeros(n*m, num_L);
if strcmp(str,'y') || ~isempty(loaded_seeds)
    fprintf('Order of marking seeds, must be the same as order of labels in the input labels vector.\n')
    color = {'m', 'yellow', 'green', 'blue', 'red'}; % assuming at most 5 regions
    
    for k = 1:num_L
        if strcmp(str,'y')
            fprintf('\nBrush on the %dth segment''s seeds:\n', k)
            fprintf('Start by clicking mouse and continue brushing and press any key to finish selecting seeds for the desired segment.')
            str_engulf = 'n';
            if num_L==2
                prompt = '\nDoes this segment engulf the other one? (y/n) ';
                str_engulf = input(prompt, 's');
            end
            figure, imshow(I, [], 'initialMagnification', 'fit'), hold on
            set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1])
            x = zeros(n*m, 1);%the position of rows of Entered seeds
            y = zeros(n*m, 1);%the position of columns of Entered seeds

            while true
                h = imfreehand('Closed', false);%you can draw desired line
                setColor(h, color{k})
                if strcmp(str_engulf,'y')
                    mask = createMask(h);%create mask with imfreehand
                    %%mask contains 1s inside the ROI and 0s everywhere else
                    boundary = bwboundaries(mask, conn);%traces the exterior boundaries of mask with 4 connected neighbouring
                    boundary = boundary{1};
                    x( nnz(x)+1:nnz(x)+size(boundary,1) ) = boundary(:,2);
                    y( nnz(y)+1:nnz(y)+size(boundary,1) ) = boundary(:,1);
                else
                    accepted_pos = getPosition(h);
                    x( nnz(x)+1:nnz(x)+size(accepted_pos,1) ) = round( accepted_pos(:,1) );
                    y( nnz(y)+1:nnz(y)+size(accepted_pos,1) ) = round( accepted_pos(:,2) );
                end
                seed_ending_flag = waitforbuttonpress;%Wait for key press or mouse-button click
                %%0 if it detects a mouse button click
                %%1 if it detects a key press then function break
                
                if seed_ending_flag==1
                    break
                end
            end
            title(['location of ', num2str(k), 'th segment''s seeds, as drawn by the user'])
            close
%             set(gcf,'Visible','off')
            x = nonzeros(x);%final position of rows of seeds
            y = nonzeros(y);%final position of columns of seeds
            seeds(1:length(x),k) = sub2ind([n,m], y, x);
            
        elseif ~isempty(loaded_seeds)
            seeds(:,k) = loaded_seeds(:,k);
            [y,x] = ind2sub( [n,m], nonzeros(seeds(:,k)) );
        end
        
     
        observed_data_pts = diag( I(y,x) );%find the location of pixels on the main image
                [prob,~] = ksdensity(observed_data_pts, I_data, 'support', 'positive');  
                %%this function find 
                %%positive= Restrict the density to positive values
                %%suport=log function for observed data
                %%default parameters: PDF,normal
        seed_dist(:,k) = prob;        
        f( sub2ind([n, m], y, x) ) = L(k);
    end    
end

regional_bias = zeros(n*m, num_L);
if ~isempty(prior_likelihood)
    switch class(prior_likelihood{1})
        case 'ProbDistUnivParam' % for the case when prior_likelihood is a cell array of ProbDist Objects
            for i = 1:num_L % num_L should also be the length of the cell array prior_likelihood
                regional_bias(:,i) = pdf( prior_likelihood{i}, I_data );%I_data:A numeric array of values where you want to evaluate the PDF.
            end
        case 'function_handle' % for the case when prior_likelihood is a cell array of custom pdf functions that I have fitted to a training sample data
            for i = 1:num_L % num_L should also be the length of the cell array prior_likelihood
                dummy_pdf = prior_likelihood{i};
                regional_bias(:,i) = dummy_pdf(I_data);
            end
    end
end


%% change labeling based on alpha expansion move at(FAST APPROXIMATE ENERGY MINIMIZATION VIA GRAPH CUT by BOYKOV 2001)

tic%starts a stopwatch timer to measure performance
min_cut_capacity = inf;%Infinity values result from operations like division by zero and overflow
cycles = 1;
iterations = 1;
success = 0;
while ~success
    fprintf( '\starting %dth cycle of Alpha Expansion algorithm.\n', cycles )
    random_order = randperm(num_L);%=[1,2](random)a certain order that can be random or fixed(page3 of fast approximate by boykov)
    random_labels = L( random_order );
    random_Data_term = seed_dist( :,random_order );
    random_seeds = seeds( :,random_order );
    random_regional_bias = regional_bias( :,random_order );
    for i = 1:num_L
        fprintf( '\starting %dth iteration of Alpha Expansion algorithm.\n', iterations )
        alpha = random_labels(i);
        fprintf( '\constructing the alpha graph for the %dth randomly chosen label.\n', i )
        tAlphaConstructStart = tic;
        [alpha_gcost, num_nodes] = AlphaGraph(I, f, camera_noise, random_labels,...
                                           random_Data_term, i, random_regional_bias);
        if any(seeds) % applying changes to the graph, to incorporate hard constraints (based on Funka-lea, Boykov 2006)
            alpha_gcost( nonzeros(seeds), end ) = 0;
            alpha_gcost( end-1, nonzeros(seeds) ) = inf;
            alpha_gcost( nonzeros( random_seeds(:,i) ) , end ) = inf;
            alpha_gcost( end-1 , nonzeros( random_seeds(:,i) ) ) = 0;
        end
        fprintf( '-->alpha graph constructed. it took %f seconds. applying a max-flow algorithm to it.\n',toc(tAlphaConstructStart) )
        tMaxFlowStart = tic;
%         [~, states, new_min_cut_capacity] = fDinic(num_nodes, alpha_gcost);
        [~, states, new_min_cut_capacity] = maxflow_Dinic(num_nodes, alpha_gcost);
%         [~, states, new_min_cut_capacity, ~] = flabeling_method_alpha_graph(num_nodes, alpha_gcost); % "states" is a vector with values 1 and 2.
                                                                                                     % these must be converted to the labels in "f"
        fprintf('-->alpha graph cut. it took %f seconds.\n', toc(tMaxFlowStart) )
        
        if new_min_cut_capacity < min_cut_capacity
            fprintf( '--->%dth iteration yielded a local minimum within one alpha-expansion of current labeling.\n', iterations )
            states = states(1:m*n); % Because only the first m*n rows of alpha_gcost are pixels and the rest of rows correspond to labels of auxillary nodes
            f( states==1 ) = alpha; % Cf. the aforementioned 2001 paper: nodes cut from, and not connected to, alpha will be so-labeled
                                    % Though f is a matrix and states is a vecotor, this assignment "does" work
            success = 1;
            min_cut_capacity = new_min_cut_capacity;
        end
        iterations = iterations+1;
    end
    cycles = cycles+1;
    iterations = 1;
    success = ~success;
end

fprintf( '\n%d cycles were run before convergance.\n', cycles-1 )
f( isnan(I_data) ) = -1; % -1 is just a dummy label for NaN-valued pixels
Labels = f;
final_alpha_gcost = alpha_gcost;
toc

segind = Labels; % Labels has the elements of the input L vector as its elements. segind changes them to a simple number from 1 to num_L for every segment
% hsubplot_seed_segres = figure; subplot(1,2,1)
% figure, imshow(I, [], 'initialMagnification', 'fit'), title('location of all segments'' seeds, as digitized by the computer')
% set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1])
% hold on
for i = 1:num_L
    segind(segind==L(i)) = i;
    
%     dummy = nonzeros( seeds(:,i) );
%     [dummy_row, dummy_col] = ind2sub( size(I), dummy);
%     plot(dummy_col, dummy_row, 'Marker','o','MarkerFaceColor',color{i},'MarkerSize',4,'Color','none')
end
% hold off
segind(segind==-1) = 0;
set( findobj('Visible','off'),'Visible','on' )
% figure, imshow(segind, []), title('segmentation result')
% set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1])
% figure(hsubplot_seed_segres), subplot(1,2,2)
% figure, imshow( label2rgb(segind, [0 1 1; 1 0 1; 1 1 0], 'k'),'initialMagnification','fit' ), title('segmentation result (color-coded)')
% set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1])

format loose