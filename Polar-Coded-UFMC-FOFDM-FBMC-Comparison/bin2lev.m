function out=bin2lev(data_in,M)
DM=log(M)/log(2);
nob=length(data_in);
a=reshape(data_in,nob/DM,DM);
out=bi2de(a);

if M==2
    out=data_in;
end