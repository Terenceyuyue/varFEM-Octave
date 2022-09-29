% Ref: �ϰ�Բ������
% dline: ���ҵ����� (�ߵ���������������ڲ�����ఴ��ʱ�����)

function [fd,fh,BdBox,pfix] = Upper_Circle_Circle_Domain
BdBox = [-1 1 0 1];
fd = @DistFnc;
fh = @huniform;
pfix = [-1 0; -0.5 0; 0.5 0; 1 0];

    function Dist = DistFnc(p)
        d1 = dcircle(p,0,0,1);
        d2 = dcircle(p,0,0,0.5);
        d3 = dline(p,0,0,1,0);
        Dist = dintersect(d3,ddiff(d1,d2));
    end

end