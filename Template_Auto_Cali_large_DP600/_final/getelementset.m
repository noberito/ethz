function ELEMENT_N = getelements(name,inp);
% file-import - Abaqus *.inp file
fid = fopen(inp, 'r');
if (fid<0)
    warndlg('File broken or not available!',...
        'File ERROR');
end

% Read Input file and count number of lines (N_LINES)
node_area = 0;
N_LINES=0;
tline = fgetl(fid);
while(node_area<2)
    if (tline == -1) % end of file
        break;
    elseif (node_area == 0 && ~strcmp(tline, name))
        tline = fgetl(fid);
    elseif (strcmp(tline, name))
        node_area = node_area + 1;
        tline = fgetl(fid);
    elseif (node_area == 1 && ~strncmp(tline, '*', 1))
        N_LINES = N_LINES + 1;
        tline = fgetl(fid);
    elseif (node_area == 1 && strncmp(tline, '*', 1))
        node_area = 2;
        break;
    end
end
frewind(fid);

% Read Input file and write ELEMENT NO 
ELEMENT_N = zeros(N_LINES*16,1);
node_area = 0;
C_line = 1;
tline = fgetl(fid);
eval(['generate = ''' num2str(name) ', generate'';']);
while(node_area<2)
    if (tline == -1) % end of file
        break;
    elseif (node_area == 0 && ~strcmp(tline, name) && ~strcmp(tline, generate))
        tline = fgetl(fid);
    elseif (strcmp(tline, name) )
        node_area = node_area + 1;
        tline = fgetl(fid);
    elseif (strcmp(tline, generate) )
        node_area = 1.5;
        tline = fgetl(fid);    
    elseif (node_area == 1.5 && ~strncmp(tline, '*', 1))
        C = sscanf(tline, '%f,%f,%f');
        for i=C(1):C(3):C(2)
          ELEMENT_N(i,1) = C(1)+(i-1)*C(3);
        end
        node_area = 2;
        break;    
    elseif (node_area == 1 && ~strncmp(tline, '*', 1))
        C = sscanf(tline, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f');
        for i=1:numel(C)
          ELEMENT_N(i+((C_line-1)*16),1) = C(i);
        end
        C_line = C_line + 1;
        tline = fgetl(fid);
    elseif (node_area == 1 && strncmp(tline, '*', 1))
        node_area = 2;
        break;
    end
end
fclose(fid); 
bb = any(ELEMENT_N(:,1)~=0,2);
ELEMENT_N = ELEMENT_N(bb,:);