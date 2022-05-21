nstart=1;
nend=9;
maxv=0;

nmats=nend-nstart+1;
rs=0;

nrows=ceil(sqrt(nmats));
colormap(gray);

load spindiff;
for n=nstart:nend
name=['density_' num2str(n) '_1;'];
eval(['data=' name ';']);


if (rs==0)
   rs=1;
   cs=1;
   for j=1:length(data)
      rs=rs+size(data{j},1);
      cs=cs+size(data{j},2);
   end
   disp(['Rows: ' num2str(rs)]);
   disp(['Cols: ' num2str(cs)]);
end

A=zeros(rs,cs);
currow=1;
curcol=2;
for j=1:length(data)
   cr=size(data{j},1);
   cc=size(data{j},2);
   rowr=[currow:currow+cr-1];
   colr=[curcol:curcol+cc-1];
   A(rowr,colr)=data{j};
   A(colr,rowr)=data{j}';
   currow=currow+cr;
   curcol=curcol+cc;
end
%disp(A);
%norm(A)
A=abs(A);
if maxv==0
   maxv=max(max(A));
end
curp=n-nstart+1;
%subplot(nrows,nrows,curp);
%axis square;
imagesc(A,[0 maxv]);
M(curp)=getframe;
end




   



