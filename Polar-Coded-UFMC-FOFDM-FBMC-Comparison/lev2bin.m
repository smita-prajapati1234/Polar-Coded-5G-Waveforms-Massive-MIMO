function out=lev2bin(data_in,M)
DM=log(M)/log(2);
nob=length(data_in);
a=de2bi(data_in,DM);
out=reshape(a,nob*DM,1)';
if M==2
    out=data_in;
end