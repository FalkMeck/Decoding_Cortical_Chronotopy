files=dir('H:\04.06\Slideshow\stimulus');

filename= sprintf('pics7.txt');
myFileID=fopen(filename,'a');

for k=3:247
    name=string(files(k).name);       
    fprintf(myFileID,name);
            fprintf(myFileID,'\n');
end
fclose(myFileID);