function [name]=findname(filename)
index=strfind(filename,'/');
sindex=index(size(index,2));
index=strfind(filename,'.');
eindex=index(size(index,2));
name=filename(sindex+1:eindex-1);
end