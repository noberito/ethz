function [ELEMENTS] = getelements(inp);
% file-import - Abaqus *.inp file
fid = fopen(inp, 'r');
if (fid<0)
    warndlg('File broken or not available!',...
        'File ERROR');
end

% Read Input file and count node-number (N_NODES)
node_area = 0;
N_LINES=0;
tline = fgetl(fid);
while(node_area<2)
    if (tline == -1) % end of file
        break;
    elseif (node_area == 0 && ~strcmp(tline, '*Element, type=C3D8R, elset= ELALL'))
        tline = fgetl(fid);
    elseif (strcmp(tline, '*Element, type=C3D8R, elset= ELALL'))
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

% Read Input file and write NODES incl. coordinates
ELEMENTS = zeros(N_LINES,9);
node_area = 0;
C_line = 1;
tline = fgetl(fid);
while(node_area<2)
    if (tline == -1) % end of file
        break;
    elseif (node_area == 0 && ~strcmp(tline, '*Element, type=C3D8R, elset= ELALL'))
        tline = fgetl(fid);
    elseif (strcmp(tline, '*Element, type=C3D8R, elset= ELALL'))
        node_area = node_area + 1;
        tline = fgetl(fid);
    elseif (node_area == 1 && ~strncmp(tline, '*', 1))
        C = sscanf(tline, '%d,%d,%d,%d,%d,%d,%d,%d,%d');
        for i=1:numel(C)
          ELEMENTS(C_line,i) = C(i);
        end
        C_line = C_line + 1;
        tline = fgetl(fid);
    elseif (node_area == 1 && strncmp(tline, '*', 1))
        node_area = 2;
        break;
    end
end
fclose(fid); 