function [pmarks,rads,Is_neg,Is_pos,depths,dtIs_neg,dtIs_pos,vlats_neg,vlats_pos,vlongs_neg,vlongs_pos]=pos_marks(params,dx,dy,srh,srv)
try 
dx/srv;
 catch Me
       fclose('all');
       disp(Me)
       disp('Horz dims and sampling rate are indiviseable!')
end
try 
dy/srh;
 catch Me
       fclose('all');
       disp(Me)
       disp('Vert dims and sampling rate are indiviseable!')
end

pmarks=zeros(dx/srv,dy/srh);
rads=zeros(dx/srv,dy/srh);
Is_neg=zeros(dx/srv,dy/srh);
Is_pos=zeros(dx/srv,dy/srh);
depths=zeros(dx/srv,dy/srh);
dtIs_neg=zeros(dx/srv,dy/srh);
dtIs_pos=zeros(dx/srv,dy/srh);
vlats_neg=zeros(dx/srv,dy/srh);
vlats_pos=zeros(dx/srv,dy/srh);
vlongs_neg=zeros(dx/srv,dy/srh);
vlongs_pos=zeros(dx/srv,dy/srh);

for ii=1:length(params(1,:));
    c1=round(params(1,ii)/srv);
    c2=round(params(2,ii)/srh);
    if c1==0
        c1=1;
    end
    if c2==0
        c2=1;
    end
    pmarks(c1,c2)=1; 
    rads(c1,c2)=params(3,ii);
    if params(4,ii)<0
        Is_neg(c1,c2)=params(4,ii);
    else
        Is_pos(c1,c2)=params(4,ii);
    end
    depths(c1,c2)=params(7,ii);
    if params(8,ii)<0
        dtIs_neg(c1,c2)=params(8,ii);
    else
        dtIs_pos(c1,c2)=params(8,ii);
    end
    if params(9,ii)<0
        vlats_neg(c1,c2)=params(9,ii);
    else
        vlats_pos(c1,c2)=params(9,ii);
    end
    if params(10,ii)<0
        vlongs_neg(c1,c2)=params(10,ii);
    else
        vlongs_pos(c1,c2)=params(10,ii);
    end
end
end

