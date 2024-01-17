function [partable]=params_loop(table)
%table='range_14.txt';
%table
fileID=fopen(num2str(table),'r');
for n=1:2
line=fgetl(fileID);
if n==1
    linetmp=line;
end
continue
end

colheaders = numel(regexp(linetmp, '[darRiIhHtTPp]'));

H.data = linetmp;

header = cell(1,colheaders);
disp(header)

%strsplit(A(1).data)

%strsplit(A(n).data)
header(:) = strsplit(H.data);
headref=['d';'a';'r';'R';'i';'I';'h';'H';'t';'T';'P';'p'];
    
for hh=1:length(header)
        for rr=1:length(headref)
            if header{hh}==headref(rr)
                header{hh}=num2str(rr);
            end
        end
end
%header = cellfun(@str2double,header);
cols = numel(regexp(line, '\w*0\w*'));
m=1;
while ischar(line)
  line = fgetl(fileID);
  A(m).data = line;
  m = m+1;

end

fclose(fileID);

B = cell(length(A),cols);

for n = 1:length(A)-1
    A(n).data
    B(n,:) = strsplit(A(n).data);
end
B=[header;B];
B = cellfun(@str2double,B);
B=B';
B=sortrows(B);
B=B';

j=length(B(:,1))-1;
B=B(2:j,:);

% %Range definitions
 pmaxv=zeros(j,1);
 pminv=zeros(j,1);

 jj=1; kk=1;
for ii=1:j-1
    if mod(ii,2)==0
    pmaxv(jj:jj+cols-1)=B(ii,:)';
    jj=jj+cols;
    end
    if mod(ii,2)==1
    pminv(kk:kk+cols-1)=B(ii,:)';
    kk=kk+cols;
    end
end



partable=cat(2,pmaxv,pminv);


end