function [NODES] = getnodes(inp);
% file-import - Abaqus *.inp file
fid = fopen(inp, 'r');
if (fid<0)
    warndlg('File broken or not available!',...
        'File ERROR');
end

% Read Input file and count node-number (N_NODES)
node_area = 0;
N_NODES=0;
tline = fgetl(fid);
while(node_area<2)
    if (tline == -1) % end of file
        break;
    elseif (node_area == 0 && ~strcmp(tline, '*Node, nset= NALL'))
        tline = fgetl(fid);
    elseif (strcmp(tline, '*Node, nset= NALL'))
        node_area = node_area + 1;
        tline = fgetl(fid);
    elseif (node_area == 1 && ~strncmp(tline, '*', 1))
        N_NODES = N_NODES + 1;
        tline = fgetl(fid);
    elseif (node_area == 1 && strncmp(tline, '*', 1))
        node_area = 2;
        break;
    end
end
frewind(fid);

% Read Input file and write NODES incl. coordinates
NODES = zeros(N_NODES,4);
node_area = 0;
C_line = 1;
tline = fgetl(fid);
while(node_area<2)
    if (tline == -1) % end of file
        break;
    elseif (node_area == 0 && ~strcmp(tline, '*Node, nset= NALL'))
        tline = fgetl(fid);
    elseif (strcmp(tline, '*Node, nset= NALL'))
        node_area = node_area + 1;
        tline = fgetl(fid);
    elseif (node_area == 1 && ~strncmp(tline, '*', 1))
        C = sscanf(tline, '%d,%f,%f,%f');
        for i=1:numel(C)
          NODES(C_line,i) = C(i);
        end 
        C_line = C_line + 1;
        tline = fgetl(fid);
    elseif (node_area == 1 && strncmp(tline, '*', 1))
        node_area = 2;
        break;
    end
end
fclose(fid); 