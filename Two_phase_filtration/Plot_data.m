clc;
close all;
clear variables;

A = load('output_copy.txt');
x = A(1,:);

fd = figure(1);
v = VideoWriter('Sw.avi','Motion JPEG AVI');
v.Quality = 95;
v.FrameRate = 20;
open(v);

for i=2:length(A(:,1))
    plot(x,A(i,:));
    grid on;
    xlabel('x,m');
    ylabel('Sw');
    title(ceil(i*1000/60/60));
    
    F=getframe(fd);
     writeVideo(v,F);
    
    pause(1e-6);
end

close(v);