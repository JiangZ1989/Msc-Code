function output = FIR(input,h_wind)
global buffer;

buffer(2:end)=buffer(1:end-1);
buffer(1)=input;
output =sum(buffer.*h_wind);

end