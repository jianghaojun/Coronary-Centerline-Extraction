function [] = coronary_refine( rpath )
% This function processes each probability map of coronary arteries under
% 'rpath'. The processing steps include but not limited to thresholding, 
% filling holes, thinning, detecting bifurcation or end points, 
% reconnecting disconnected branches, removing isolated branches, and 
% obtains a coronary artery tree finally.
% 
% Examples:
%   coronary_refine('path_of_parent_directory_containg_volumes')

% create output directory
% rpath = 'C:\Users\hjjiang\Desktop\大三上\数字图像处理\作业\综合作业2\code';
wpath = fullfile(rpath, 'coronary');
if ~exist(wpath, 'dir'), mkdir(wpath); end

% processing each volume under rpath
img_list = dir(fullfile(rpath, '*.mha'));
binary_threshold = [0.5, 0.70];
for ii = 1:length(img_list)
    %% read mha volume
    img_path = fullfile(rpath, img_list(ii).name);
    [img_prop, img_info] = mha_read_volume(img_path); %img_prop 是指图像每个位置属于冠脉的可能性
    
    %% thresholdingCoro
    img_bin = img_prop >= (binary_threshold(ii)*intmax('uint16'));
    % check the binary image obtained by considering it as a volume data,
    % and you can also store the binary volume into a single file (.mha file)
    % volumeViewer(img_bin);
    w_info = img_info; % header information of volume to be written
    w_info.DataType = 'uchar'; % change the data type to uchar (uint8)
    mha_write(img_bin, w_info, 'path'); %将图片信息写入到名为path的文件中去
    
    %% filling holes
    img_bin = imfill(img_bin, 'holes');
    
    %% thinning
    img_thin = bwskel(img_bin);
    img_thin = bwareaopen(img_thin, 10); 
    
    %% detecting bifucation and end points
    img_branchpoints = bwmorph3(img_thin, 'branchpoints'); % 找到分叉点
    [branchpoints_x, branchpoints_y, branchpoints_z] = ind2sub(size(img_branchpoints), find(img_branchpoints));
    
    img_endpoints = bwmorph3(img_thin, 'endpoints');% 找到端点
    [endpoints_x, endpoints_y, endpoints_z] = ind2sub(size(img_endpoints), find(img_endpoints));
    
    %% reconnecting disconnected branches & removing isolate points or branches
    % connectivity label 如果属于一个连通对象，对应的标签也是相同的
    connectivity_label = zeros(size(img_thin));
    endpoints_connectivity_label = zeros(size(endpoints_x));
    CC = bwconncomp(img_thin);
    breakbranch_label = zeros(size(CC.PixelIdxList, 2), 1);
    for i = 1:size(CC.PixelIdxList, 2)
        [X,Y,Z] = ind2sub(size(img_thin), CC.PixelIdxList{i});
        for j = 1:size(X, 1)
            connectivity_label(X(j),Y(j),Z(j)) = i;
        end
        if(size(X, 1) < 100)
            breakbranch_label(i) = 1; % 判断断开的连通对象可能是断裂分支，一般断裂分支长度短于100
        end
    end
    
    for i = 1:size(endpoints_x, 1)
        tmp_x = endpoints_x(i);
        tmp_y = endpoints_y(i);
        tmp_z = endpoints_z(i);
        endpoints_connectivity_label(i) = connectivity_label(tmp_x,tmp_y,tmp_z);
    end
    
    % 断连分支连接，先用端点和分叉点的信息进行修复；有的断连分支明显是主干的分支，但是可能没有对应的分叉点，选取主干上离分支端点最近的一点进行直线连接
    branchendpoints_x = [endpoints_x' branchpoints_x']; %分叉点和端点合并
    branchendpoints_y = [endpoints_y' branchpoints_y'];
    branchendpoints_z = [endpoints_z' branchpoints_z'];
    radius = [35, 25]; % 考虑距离端点半径范围内的端点和分叉点
    tmp_x = 0;
    tmp_y = 0;
    tmp_z = 0;
    for i = 1:size(endpoints_x, 1)
        tmp_prop = img_prop;
        x = endpoints_x(i);
        y = endpoints_y(i);
        z = endpoints_z(i);
        correspondingpoint_index = 0; % 如果该端点是断连点，它对应的连接点
        min_distance = 1000000;
        distance = 0;
        if(breakbranch_label(connectivity_label(x,y,z)) == 1) % 该端点属于断裂分支的端点
            for j = 1:size(branchendpoints_x, 2)
                if(j ~= i)
                    tmp_x = branchendpoints_x(j);
                    tmp_y = branchendpoints_y(j);
                    tmp_z = branchendpoints_z(j);
                    if(connectivity_label(x,y,z) ~= connectivity_label(tmp_x,tmp_y,tmp_z))% 利用connectivity_label 判断两个点是否连通 如果不连通才可能是两个断连的点 
                        if(branchendpoints_x(j) - x < radius(ii))
                            distance = sqrt((x - tmp_x)^2 + (y - tmp_y)^2 + (z - tmp_z)^2);
                            if(distance < radius(ii) && distance < min_distance)
                                min_distance = distance;
                                correspondingpoint_index = j;
                            end
                        end
                    end
                end
            end
        end
        
        % 判断断连分支另一端是否离correspondingpoint更近，如果更近则直接跳过该端点的连接，在后面会用更近的端点连接
        if(correspondingpoint_index ~= i && correspondingpoint_index ~= 0)
            index = find(endpoints_connectivity_label == connectivity_label(x,y,z));
            for k = 1:size(index,1)
                if(index(k) ~= i)
                    distance = sqrt((endpoints_x(index(k)) - branchendpoints_x(correspondingpoint_index))^2 + (endpoints_y(index(k)) - branchendpoints_y(correspondingpoint_index))^2 + (endpoints_z(index(k)) - branchendpoints_z(correspondingpoint_index))^2);
                    if(distance < min_distance)
                        correspondingpoint_index = 0;
                        break;
                    end
                end
            end
        end
        
        % 断支连接
        if(correspondingpoint_index ~= i && correspondingpoint_index ~= 0)
            start_point = [endpoints_x(i) endpoints_y(i) endpoints_z(i)];
            target_point = [branchendpoints_x(correspondingpoint_index) branchendpoints_y(correspondingpoint_index) branchendpoints_z(correspondingpoint_index)];
            index = 1;
            min_distance = sqrt((start_point(1) - target_point(1))^2 + (start_point(2) - target_point(2))^2 + (start_point(3) - target_point(3))^2);
            neighbor_coordinate = zeros(27,3);
            while(min_distance > sqrt(3))
                % 找到 tmp_point 周围3*3*3块中距离目标点最近的一点
                for k = -1:1
                    for m = -1:1
                        for n = -1:1
                            neighbor_coordinate(((k+1)*9+(m+1)*3+(n+2)),1) = start_point(1) + k; 
                            neighbor_coordinate(((k+1)*9+(m+1)*3+(n+2)),2) = start_point(2) + m;
                            neighbor_coordinate(((k+1)*9+(m+1)*3+(n+2)),3) = start_point(3) + n;
                        end
                    end
                end
                
                alldistance = sqrt(sum(bsxfun(@minus, neighbor_coordinate, target_point).^2, 2));
                index = find(alldistance == min(alldistance));
                
                start_point(1) = neighbor_coordinate(index(1),1);
                start_point(2) = neighbor_coordinate(index(1),2);
                start_point(3) = neighbor_coordinate(index(1),3);
                img_thin(start_point(1),start_point(2),start_point(3)) = 1;
                if(alldistance(index(1)) < min_distance)
                    min_distance = alldistance(index(1));
                end
            end
            
            % 更新connectivity_label和endpoints_connectivity_label
            CC = bwconncomp(img_thin);
            for p = 1:size(CC.PixelIdxList, 2)
                [X,Y,Z] = ind2sub(size(img_thin), CC.PixelIdxList{p});
                for q = 1:size(X, 1)
                    connectivity_label(X(q),Y(q),Z(q)) = p;
                end
                if(size(X, 1) < 100)
                    breakbranch_label(p) = 1; % 判断断开的连通对象可能是断裂分支，一般断裂分支长度短于100
                end
            end
            for p = 1:size(endpoints_x, 1)
                tmp_x = endpoints_x(p);
                tmp_y = endpoints_y(p);
                tmp_z = endpoints_z(p);
                endpoints_connectivity_label(p) = connectivity_label(tmp_x,tmp_y,tmp_z);
            end
        end
    end
    
    % 进一步修复
    img_endpoints = bwmorph3(img_thin, 'endpoints');% 更新端点
    [endpoints_x, endpoints_y, endpoints_z] = ind2sub(size(img_endpoints), find(img_endpoints));
    endpoints_connectivity_label = zeros(size(endpoints_x));
    CC = bwconncomp(img_thin);
    breakbranch_label = zeros(size(CC.PixelIdxList, 2), 1);
    for i = 1:size(CC.PixelIdxList, 2)
        [X,Y,Z] = ind2sub(size(img_thin), CC.PixelIdxList{i});
        for j = 1:size(X, 1)
            connectivity_label(X(j),Y(j),Z(j)) = i;
        end
        if(size(X, 1) < 100)
            breakbranch_label(i) = 1; % 判断断开的连通对象可能是断裂分支，一般断裂分支长度短于100
        end
    end
    
    for i = 1:size(endpoints_x, 1)
        tmp_x = endpoints_x(i);
        tmp_y = endpoints_y(i);
        tmp_z = endpoints_z(i);
        endpoints_connectivity_label(i) = connectivity_label(tmp_x,tmp_y,tmp_z);
    end
    
    for i = 1:size(endpoints_x, 1)
        min_distance = 1000000;
        distance = 0;
        find_targetpoint_flag = 0;
        if(breakbranch_label(endpoints_connectivity_label(i)) == 1) % 如果端点属于断裂分支
            for dx = max(round(endpoints_x(i)- radius(ii)/2), 1):min(round(endpoints_x(i) + radius(ii)/2), size(img_thin, 1))
                for dy = max(round(endpoints_y(i) - radius(ii)/2), 1):min(round(endpoints_y(i) + radius(ii)/2), size(img_thin, 2))
                    for dz = max(round(endpoints_z(i) - radius(ii)/2), 1):min(round(endpoints_z(i) + radius(ii)/2), size(img_thin, 3))
                        if(img_thin(dx, dy, dz) == 1 && connectivity_label(dx, dy, dz) ~= endpoints_connectivity_label(i))
                            distance = sqrt((dx - endpoints_x(i))^2 + (dy - endpoints_y(i))^2 + (dz - endpoints_z(i))^2);
                            if(distance < min_distance && distance < radius(ii))
                                min_distance = distance;
                                target_point(1) = dx;
                                target_point(2) = dy;
                                target_point(3) = dz;
                                find_targetpoint_flag = 1;
                            end
                        end
                    end
                end
            end
            
            if(find_targetpoint_flag == 1)
                start_point = [endpoints_x(i) endpoints_y(i) endpoints_z(i)];
                index = 1;
                min_distance = sqrt((start_point(1) - target_point(1))^2 + (start_point(2) - target_point(2))^2 + (start_point(3) - target_point(3))^2);
                neighbor_coordinate = zeros(27,3);
                while(min_distance > sqrt(3))
                    % 找到 tmp_point 周围3*3*3块中距离目标点最近的一点
                    for k = -1:1
                        for m = -1:1
                            for n = -1:1
                                neighbor_coordinate(((k+1)*9+(m+1)*3+(n+2)),1) = start_point(1) + k; 
                                neighbor_coordinate(((k+1)*9+(m+1)*3+(n+2)),2) = start_point(2) + m;
                                neighbor_coordinate(((k+1)*9+(m+1)*3+(n+2)),3) = start_point(3) + n;
                            end
                        end
                    end

                    alldistance = sqrt(sum(bsxfun(@minus, neighbor_coordinate, target_point).^2, 2));
                    index = find(alldistance == min(alldistance));

                    start_point(1) = neighbor_coordinate(index(1),1);
                    start_point(2) = neighbor_coordinate(index(1),2);
                    start_point(3) = neighbor_coordinate(index(1),3);
                    img_thin(start_point(1),start_point(2),start_point(3)) = 1;
                    if(alldistance(index(1)) < min_distance)
                        min_distance = alldistance(index(1));
                    end
                end
            end
        end
    end
    
    % 清除孤立点和分支
    img_thin = bwareaopen(img_thin, 100);
    
    %% obtain coronary artery tree (by tracking or other methods)
    % 更新分支点和端点
    img_branchpoints = bwmorph3(img_thin, 'branchpoints'); % 找到分叉点
    [branchpoints_x, branchpoints_y, branchpoints_z] = ind2sub(size(img_branchpoints), find(img_branchpoints));
    
    img_endpoints = bwmorph3(img_thin, 'endpoints');% 找到端点
    [endpoints_x, endpoints_y, endpoints_z] = ind2sub(size(img_endpoints), find(img_endpoints));
    
    % 寻找所有的分支并对分支的点进行排序(空间相邻，存储也要相邻)
    coronary_tree = {};
    tmp_img_thin = img_thin;
    while(isempty(find(tmp_img_thin)) == 0) % 如果tmp_img_thin中还有点没有被分到相应的分支
        % 要更新端点和分叉点，因此算法中原先的端点会消失，而部分分叉点会变成端点
        img_branchpoints = bwmorph3(tmp_img_thin, 'branchpoints'); % 找到分叉点
        [branchpoints_x, branchpoints_y, branchpoints_z] = ind2sub(size(img_branchpoints), find(img_branchpoints));

        img_endpoints = bwmorph3(tmp_img_thin, 'endpoints');% 找到端点
        [endpoints_x, endpoints_y, endpoints_z] = ind2sub(size(img_endpoints), find(img_endpoints));
        
        [endpoints_z, index] = sort(endpoints_z);
        endpoints_y = endpoints_y(index);
        endpoints_x = endpoints_x(index);
        
        [branchpoints_z, index] = sort(branchpoints_z);
        branchpoints_y = branchpoints_y(index);
        branchpoints_x = branchpoints_x(index);
        
        for i = 1:size(endpoints_x, 1)
            branch_coordinate = cell(1);
            find_branchpoint_flag = 0;
            %find_all_nearby_branchpoint_flag = 0;
            tmp_x = endpoints_x(i);
            tmp_y = endpoints_y(i);
            tmp_z = endpoints_z(i);
            while(find_branchpoint_flag ~= 1)
                branch_coordinate{1} = [branch_coordinate{1};[tmp_x, tmp_y, tmp_z]];
                tmp_img_thin(tmp_x, tmp_y, tmp_z) = 0;
                neighbor = tmp_img_thin(tmp_x-1:tmp_x+1, tmp_y-1:tmp_y+1, tmp_z-1:tmp_z+1);
                [neighbor_x, neighbor_y, neighbor_z] = ind2sub(size(neighbor), find(neighbor));
                if(size(neighbor_x, 1) > 1) % 该点3*3*3周围有至少2个点，可能分支点，将坐标与分支点对比
                    if(isempty(find(branchpoints_x == tmp_x)) == 0)
                        if(isempty(find(branchpoints_y == tmp_y)) == 0)
                            if(isempty(find(branchpoints_z == tmp_z)) == 0)
                                find_branchpoint_flag = 1; % tmp_x 这点就是分支点
                                % 将该分支点附近的分支点全部纳入该分支中，即原始模型从此处断开
                                neighbor = tmp_img_thin(tmp_x-1:tmp_x+1, tmp_y-1:tmp_y+1, tmp_z-1:tmp_z+1);
                                [neighbor_x, neighbor_y, neighbor_z] = ind2sub(size(neighbor), find(neighbor));
                                for j = 1:size(neighbor_x, 2)
                                    if(isempty(find(branchpoints_x == tmp_x + neighbor_x(j) - 2)) == 0)
                                        if(isempty(find(branchpoints_y == tmp_y + neighbor_y(j) - 2)) == 0)
                                            if(isempty(find(branchpoints_z == tmp_z + neighbor_z(j) - 2)) == 0)
                                                tmp_img_thin(tmp_x + neighbor_x(j) - 2, tmp_y + neighbor_y(j) - 2, tmp_z + neighbor_z(j) - 2) = 0;
                                                branch_coordinate{1} = [branch_coordinate{1};[tmp_x + neighbor_x(j) - 2, tmp_y + neighbor_y(j) - 2, tmp_z + neighbor_z(j) - 2]];
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if(find_branchpoint_flag == 0) % 可能分支中间有点聚集
                        tmp_x = tmp_x + neighbor_x(1) - 2;
                        tmp_y = tmp_y + neighbor_y(1) - 2;
                        tmp_z = tmp_z + neighbor_z(1) - 2;
                    end
                end
                if(size(neighbor_x, 1) == 1) % 该点周围只有一个点，不可能是分支点，继续追踪下去
                    tmp_x = tmp_x + neighbor_x(1) - 2;
                    tmp_y = tmp_y + neighbor_y(1) - 2;
                    tmp_z = tmp_z + neighbor_z(1) - 2;
                end
                if(size(neighbor_x, 1) == 0) % 该点周围没有一个点，那么该段终止
                    find_branchpoint_flag = 1;
                end
            end
            if(size(branch_coordinate{1}, 1) ~= 1)
                coronary_tree{end+1} = branch_coordinate{1};
            end
        end
    end
    
    % plot the complete coronary artery tree in different color according 
    % to the ids of branches, e.g. denote 'coro_tree', a cell array, as the 
    % coronary artery tree obtained, and each element is a coordinate array
    % of a single branch
    coronary_show(coronary_tree);
    
    %% save the tree obtained into a mat file (.mat)
    tree_name = split(img_list(ii).name, '.');
    tree_name = [tree_name{1}, '.mat'];
    tree_path = fullfile(wpath, tree_name);
    save(tree_path, 'coronary_tree');
end
end