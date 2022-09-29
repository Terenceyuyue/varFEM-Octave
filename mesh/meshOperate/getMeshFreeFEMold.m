function [node,elem] = getMeshFreeFEMold(meshName)

% e.g. meshName = 'meshdata_heatex.msh'
data = importdata(meshName);

%% Numbers
N = data(1,1); NT = data(1,2);
data(1,:) = []; % delete fist line

%% node
node = data(1:N,1:2);
data(1:N,:) = [];

%% elem
% �ò��������� 4 �У�matlab ��� 4 ����� NaN ����� 6 ��Ԫ�ط�Ϊ���У��Ա�֤
% ��ǰ�����ж�Ӧ��Ϊ�˸ò���ÿ���ж�Ӧһ����Ԫ
elem = data(1:2:2*NT,:);
% data(1:2*NT,:) = [];

% %% bdEdge
% ind = data(:,3);
% label = unique(ind);
% bdEdgeIdxType = cell(length(label),1);
% for i = 1:length(label)
%     bdEdgeIdxType{i} = data(ind==i,1:2);
% end

% meshOutputName = meshName(1:end-4);
% save(meshOutputName,'node','elem','bdEdgeIdxType')
% showmesh(node,elem)