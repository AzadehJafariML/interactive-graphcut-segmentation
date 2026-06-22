function [alpha_gcost, num_nodes] = AlphaGraph(Image, f, camera_noise, randlabels, randSeedDist, currlabind, randRegionalBias)
%%alpha graph(FAST APPROXIMATE ENERGY MINIMIZATION VIA GRAPH CUT by BOYKOV 2001 page5 part5)

%% Explanation of parameters:
%%Image: grayscaled and double-classed image to be segmented
%%f: current labels of all pixels created by Kmeans
%%camera_noise: estimated by mean of intensity diffrences between adjacenet pixels
%%randlabels:set of labels, randomly ordered for this cycle of alpha_expansion move
%%randseedDist:the Data_term (the same as the old fD function), with its columns ordered accroding to randlabels
%%currentlabind:the index of the label in randlabels that is the ALPHA terminal for this iteration of alpha_expansion
%%randregionalBias:regional_bias (prior likelihoods of pixels belonging to certain labels), with its columns ordered accroding to randlabels
alpha = randlabels(currlabind);
[n,m] = size(Image);
A = 1;
dist_flag = 1;
gamma = .5;
lambda =2;%constant that is determined by changes and runing code for each value
if any(randRegionalBias) %if any array elements are nonzero
    for i = 1:size(randRegionalBias,2)
        if unique(randRegionalBias(:,i),'legacy')==1 % this makes sure that if prior knowledge for a segment is not known, (in which case we use p=1 for all pixels so that log(p)=0], regional cost...
            randRegionalBias(:,i) = -lambda*(1-gamma)*log(randSeedDist(:,i)); % for that segment is not 0 for all pixels. rather, regional cost is chosen to be the probability deriven from seeds
        else
            % randRegionalBias = uplim.*randRegionalBias;
            randRegionalBias(:,i) = -lambda*(1-gamma)*log(randRegionalBias(:,i)); % if because of not entering seeds, randDterm was zeros, this command would make it inf, which would impair our algorithm
            % randRegionalBias = -log10(randRegionalBias);
        end
    end
    if 1-gamma==0  % to solve the problem of 0*inf = NaN, which happens when 1-gamma=0 and p=0 so that log0=inf
        randRegionalBias = zeros(size(randRegionalBias));
    end
else
    gamma = 1; % this is for when we don't provide any regional bias term (which means it is zeros(n*m,num_L)). In that case gamma should be 1 so that A alone controls data term vs neighborhood term
end
if any(randSeedDist) % this makes sure that randDterm is not zeros because of the user "not marking" seeds
    % randDterm = uplim.*randDterm;
    randSeedDist = -lambda*gamma*log(randSeedDist); % if because of not entering seeds, randDterm was zeros, this command would make it inf, which would impair our algorithm
    % randDterm = -log10(randDterm);
    if gamma==0 % to solve the problem of 0*inf = NaN, which happens when gamma=0 and p=0 so that log0=inf
        randSeedDist = zeros(size(randSeedDist));
    end
end

num_pixels = n*m;
num_AN = 4*num_pixels - 3*(n+m) + 2; %AN = Auxillary Node
num_nodes = num_pixels + num_AN + 2; % this includes pixels, auxillary nodes, and terminals(of which there are always 2, alpha and alpha_bar)
srow = zeros(1,26*num_nodes); scol = zeros(1,26*num_nodes); sz = zeros(1,26*num_nodes);

I_data = Image(:);
f = f(:);
mask_pixels = [-n-1 -1 n-1; ...
               -n    0 n; ...
               -n+1 +1 n+1];
dist_coef1 = [1/sqrt(2) 1 1/sqrt(2); ... % as suggested by page 5 of boykov,jolly2001 paper and page10...
             1         1 1; ...         %  of boykov,funka-lea2006 paper, boundary penalty should...
             1/sqrt(2) 1 1/sqrt(2)];%be a function of distance between the two pixels and we have...
                                     %sqrt(2) for diagnol neighbouring
dist_coef2 = ones(3);
dist_coef = {dist_coef1, dist_coef2};
dist_coef = dist_coef{dist_flag};
  
if any( isnan(I_data) )
    f( isnan(I_data) ) = -1; % this is just a dummy label for NaN-valued pixels, because unfortunately NaN itself return zero for isequal(NaN,NaN)
    randlabels(end+1) = -1; % this is just a dummy label for NaN-valued pixels, because unfortunately NaN itself return zero for isequal(NaN,NaN)
    randSeedDist(:,end+1) = NaN; % because while you want the label for NaN-valued pixels to be "-1", you want their edge costs to be NaN
    randRegionalBias(:,end+1) = NaN;
end

scol( 1:num_pixels ) = num_nodes; % because we are going to consider the connection of "all" the pixels, to the "sink" first.
alpha_pixels = find( f==alpha );
flag = length(alpha_pixels);
if flag~=0
    srow(1:flag) = alpha_pixels;
      sz(1:flag) = inf;
end

if flag~=num_pixels
    dummy1 = randRegionalBias';
    dummy2 = randSeedDist';
    dummy3 = dummy2( (repmat(f,1,length(randlabels))==repmat(randlabels,num_pixels,1))' );%returns copies of f
    dummy4 = dummy1( (repmat(f,1,length(randlabels))==repmat(randlabels,num_pixels,1))' );
    srow( flag+1:num_pixels ) = find( f~=alpha );
      sz( flag+1:num_pixels ) = dummy3( f~=alpha ) + dummy4( f~=alpha ); % + fD( f(f~=alpha),I_data(f~=alpha) );
end

flag = 2*num_pixels;
srow( num_pixels+1:flag ) = num_nodes-1;
scol( num_pixels+1:flag ) = 1:num_pixels;
  sz( num_pixels+1:flag ) = randSeedDist(:,currlabind) + randRegionalBias(:,currlabind); % + fD(alpha, I_data);

counter1 = 0; % number of pairs of neighboring pixels which fp = fq (page 7 of the paper)
counter2 = 0; % number of pairs of neighboring pixels which fp ~= fq (page 7 of the paper)
nonNaN_pixels = find( ~isnan(I_data) )'; % we want any link weight from/to NaN pixels to be zero. so we "not consider" them in the next for loop altogether
for p = nonNaN_pixels
    row = mod(p, n); % matricial_row_of_p
    row(row==0) = n; % this is necessary for pixels in the bottom row of the image, because row=mod(n,n)==0 is meaningless
    col = ceil(p/n); % matricial_col_of_p
    nghbr_AN = fmask_AN(n, row, col);
    
    if p == 1
        maskp = mask_pixels(2:3,2:3);
        nghbr_ANp= nghbr_AN(2:3,2:3);
        dist_coefp = dist_coef(2:3,2:3);
        
    elseif p == n
        maskp = mask_pixels(1:2,2:3);
        nghbr_ANp= nghbr_AN(1:2,2:3);
        dist_coefp = dist_coef(1:2,2:3);
        
    elseif p == num_pixels-n+1
        maskp = mask_pixels(2:3,1:2);
        nghbr_ANp= nghbr_AN(2:3,1:2);
        dist_coefp = dist_coef(2:3,1:2);
        
    elseif p == num_pixels
        maskp = mask_pixels(1:2,1:2);
        nghbr_ANp= nghbr_AN(1:2,1:2);
        dist_coefp = dist_coef(1:2,1:2);
        
    elseif (p>=2) && (p<=n-1)
        maskp = mask_pixels(:,2:3);
        nghbr_ANp= nghbr_AN(:,2:3);
        dist_coefp = dist_coef(:,2:3);
        
    elseif (p>=num_pixels-n+2) && (p<=num_pixels-1)
        maskp = mask_pixels(:,1:2);
        nghbr_ANp= nghbr_AN(:,1:2);
        dist_coefp = dist_coef(:,1:2);
        
    elseif (p>=n+1) && (p<=num_pixels-n)
        if mod(p,n) == 1
            maskp = mask_pixels(2:3,:);
            nghbr_ANp= nghbr_AN(2:3,:);
            dist_coefp = dist_coef(2:3,:);
        elseif mod(p,n) == 0
            maskp = mask_pixels(1:2,:);
            nghbr_ANp= nghbr_AN(1:2,:);
            dist_coefp = dist_coef(1:2,:);
        else
            maskp = mask_pixels;
            nghbr_ANp= nghbr_AN;
            dist_coefp = dist_coef;
        end
        
    end
    
    maskp = maskp(:)+p;
    nghbr_ind = maskp( maskp~=p ); % indices of the 8 pixels in the neighborhood of pixel p
    
    dist_coefp = dist_coefp(:);
    dist_coefp = dist_coefp( maskp~=p );
    
    nghbr_ANp = nghbr_ANp(:);
    nghbr_ANp_ind = nghbr_ANp( nghbr_ANp~=0 ); % indices of the 8 auxillary nodes in neighborhood of pixel p
    
    nghbr_ANp_ind = nghbr_ANp_ind( ~isnan(I_data(nghbr_ind)) ); % removing ANs corresponding to NaN-valued pixels from the set of neighbors
    dist_coefp = dist_coefp( ~isnan(I_data(nghbr_ind)) );
    nghbr_ind = nghbr_ind( ~isnan(I_data(nghbr_ind)) ); % removing NaN-valued pixels from the set of neighbors
    
    same_label_nghbrs = nghbr_ind( f(p)==f(nghbr_ind) );
    same_label_nghbrs_dist_coef = dist_coefp( f(p)==f(nghbr_ind) );
    if ~isempty(same_label_nghbrs)
        dummy = length(same_label_nghbrs);
        srow( flag+1:flag+dummy ) = p;
        scol( flag+1:flag+dummy ) = same_label_nghbrs;
          sz( flag+1:flag+dummy ) = same_label_nghbrs_dist_coef.*fV(f(p), alpha, I_data(p), I_data(same_label_nghbrs), camera_noise, A);
        flag = flag+dummy;
        counter1 = counter1+dummy;
    end
    
    diff_label_nghbrs = nghbr_ind( f(p)~=f(nghbr_ind) );
    diff_label_nghbrs_dist_coef = dist_coefp( f(p)~=f(nghbr_ind) );
    if ~isempty(diff_label_nghbrs)
        diff_label_ANps = nghbr_ANp_ind( f(p)~=f(nghbr_ind) );
        dummy = length(diff_label_nghbrs);
        srow(flag+1:flag+2*dummy) = [repmat(p,dummy,1); num_pixels+diff_label_ANps];
        scol(flag+1:flag+2*dummy) = [num_pixels+diff_label_ANps; repmat(p,dummy,1)];
          sz(flag+1:flag+2*dummy) = repmat( diff_label_nghbrs_dist_coef.*fV(f(p), alpha, I_data(p), I_data(diff_label_nghbrs), camera_noise, A),2,1 );
        flag = flag+2*dummy;
        
        srow(flag+1:flag+dummy) = num_pixels+diff_label_ANps;
        scol(flag+1:flag+dummy) = num_nodes;
          sz(flag+1:flag+dummy) = 0.5*diff_label_nghbrs_dist_coef.*fV(f(p), f(diff_label_nghbrs), I_data(p), I_data(diff_label_nghbrs), camera_noise, A); % multiplied...
          % by 0.5 because this line will execute again when p becomes the neighbor between p (the current p) and whom this AN is being created
        flag = flag+dummy;
        counter2 = counter2+dummy;
    end

end

srow = nonzeros(srow);
scol = nonzeros(scol);
sz = sz(1:nnz(srow));
sz( isnan(sz) ) = 0;
alpha_gcost = sparse(srow, scol, sz, num_nodes, num_nodes);
% alpha_gcost = floor(alpha_gcost); % floor must happen after creating the matrix using "sparse" command, so that summing of sz-values for the connection...
% of AN nodes to sink, according to line 153, occurs first and so the rounding effect of floor doesn't reduce their final values by 1
fprintf('++>%d pairs of neighboring pixels had similar labels.\n', counter1)
fprintf('++>%d pairs of neighboring pixels didn''t have similar labels, so %d auxillary nodes were created.\n', counter2, counter2)
end

%%%%%%%%%%%%%%% Local Funcitons %%%%%%%%%%%%%%%%%%%%
function nghbr_AN = fmask_AN(n, row, col)
% row and col are the row and column of the current pixel, p, in the image
left_cons = (col-1)*(n-1) + (col-2)*(3*n-2) + (row-3)*3+2+2;
midd_cons = (col-1)*(n-1) + (col-1)*(3*n-2) + row-2;
rght_cons =     col*(n-1) + (col-1)*(3*n-2) + (row-1)*3-1;
nghbr_AN = [left_cons+1, midd_cons+1, rght_cons+1;...
            left_cons+3, 0,           rght_cons+2;...
            left_cons+5, midd_cons+2, rght_cons+3];
end

% function D = fD(fp, Ip)
% % p is a vector containing indices of the pixels
% % fp is a vector containing labels of the pixels in p
% % Ip is the intensity of the pixels in p
% % D = min( (fp - Ip).^2, K );
% D = abs( fp - Ip );
% % D = 255 - abs(fp - Ip);
% end

function V = fV(fp, fq, Ip, Iq, camera_noise, A)
%%fp is the label of the pixel p
%%fq is the label of the set of pixels q neighbor to p
%%Ip = intensity of the pixel p
%%Iq = intensity of the set of pixels q neighbor to p
%%V = K.*min(3, abs(fp-fq));
%%V = 255.*exp( -abs(fp-fq) );
%%V = abs(fp-fq);
%%camera_noise: this is a variance based on boykov and jolly, 2001, page 5

% Potts model of interaction penalty:
if length(fq) == 1
    fq = repmat(fq, length(Iq), 1);
end
V = zeros(length(fq),1);
V(fp~=fq) = A*exp( -(Ip-Iq(fp~=fq)).^2/(2*camera_noise^2) );
end