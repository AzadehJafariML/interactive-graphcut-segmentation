function [residual_graph, tree, min_cut_capacity] = maxflow_Dinic(num_nodes, alpha_graph) 
%%Dinic algorithm based on (DINITZ ALGORITHM:THE ORIGINAL VERSION AND EVEN'S VERSION)
%%1TLN stands for "1-Terminal Layered Network"


tStartDinic = tic;

t = num_nodes;
residual_graph = alpha_graph;
flag_stlen = 1; % this flag indicates that all paths of a certain length are finished, and that s-t distance has increased by 1
tree = ones(t, 1); % like the "status" variable of flabeling_method code, tree=1 means connected to sink, while tree=2 means connected to source

phase = 1;%do all steps of dinic algorithm
iteration = 1;
while true
    
    % Constructing Layered Network Data Structure
    if flag_stlen == 1
        tStartPhase = tic;
        fprintf('@@>starting Phase %d\n',phase)
        node_layer = zeros(t,1);
        edge_layer = spalloc(t,t,t);%Allocate space for sparse matrix
        
        node_layer(t-1) = 1; % actually dist of soure to source is zero, but because of matlab's practicalities we use a dist higher by 1 for all nodes
        i = 2;%start with layer with length i=2 based on BFS and then increse it
        while true
            [~,dummy] = find( residual_graph( node_layer==i-1,: )>0 );
            dummy = dummy( node_layer(dummy)==0 );
            node_layer(dummy) = i;
            
            prev_layer = find( node_layer==i-1 );
            curr_layer = find( node_layer==i );
            for j = 1:length(prev_layer)
                edge_layer( prev_layer(j) , curr_layer(residual_graph(prev_layer(j),curr_layer)>0) ) = i;
            end
            
            if ismember(t,dummy,'legacy')
                stlen = i; % stlen is the length (number of arcs+1) of the shortest path between s (source) and t (sink)
                break
            elseif isempty(dummy) % maximum flow is reached
                tree( node_layer~=0 ) = 2; % these are the nodes that are connected to the source after pushing maximum flow
                min_cut_capacity = sum( residual_graph(:,t-1) ); % remember, this is residual graph with residual capacity values, not flow values. So...
                % instead of residual_graph(t-1,:) we should use residual_graph(:,t-1)
                fprintf('@@>Maximum Flow is reached. It took %f seconds\n', toc(tStartDinic))
                return
            end
            i = i+1;
        end
    end
    
    % Finding and Pushing Flow along an Augmenting Path
    path = zeros(stlen,1);
    arch_flow = zeros(stlen,1);
    
    path(1) = t-1; path(end) = t;
    arch_flow(end) = inf;
    for j = stlen-1:-1:1
        dummy = find( edge_layer(:,path(j+1)) ); % possibly, i can remove the ==j condition
        path(j) = datasample( dummy,1 );
        arch_flow(j) = min( arch_flow(j+1),residual_graph(path(j),path(j+1)) );
    end
    
    path_arc = sub2ind( [t,t],path(1:end-1),path(2:end) );
    oppo_arc = sub2ind( [t,t],path(2:end),path(1:end-1) ); % opposite arcs of the arcs in the path
    residual_graph(path_arc) = residual_graph(path_arc) - arch_flow(1); % this might saturate because its residual capacity is being decreased
    residual_graph(oppo_arc) = residual_graph(oppo_arc) + arch_flow(1); % this can't saturate because its residual capacity is being increased
    % to cure inf-inf problem that returens NaN instead of 0:
    dummy = residual_graph(path_arc);
    dummy( isnan(dummy) ) = 0;
    residual_graph(path_arc) = dummy;
    
    % Cleaning
    flag_stlen = 0;
    
    Sat = path_arc( residual_graph(path_arc)==0 ); % the set of saturated edges
    edge_layer(Sat) = 0;
    [~,Satr] = ind2sub( [t,t],Sat );
    
    Qr = zeros(t,1); % right queue of egdes, as defined by Dinic 2006 paper
    Qr(1:length(Sat)) = Satr;
    i = 1; % Right Pass
    while true
        if ~any( edge_layer(:,Qr(i)) )
            dummy = find( edge_layer(Qr(i),:) );
            Qr(nnz(Qr)+1: nnz(Qr)+length(dummy)) = dummy;
            edge_layer(Qr(i),:) = 0;
        end
        if i == nnz(Qr)
            break
        end
        i = i+1;
    end
    
    if ~any( edge_layer(:,t) )
        fprintf('@@>Phase %d consisted of %d iterations, cumulatively taking %f seconds\n', phase, iteration, toc(tStartPhase))
        phase = phase+1;
        flag_stlen = 1;
        iteration = 1;
    end
    
    iteration = iteration+1;
end