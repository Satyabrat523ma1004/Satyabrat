function pg32_lineset_search_5_16
clc;

%% =========================
%  POINTS OF PG(3,2)
%% =========================

F = [0 1];
P = [];
for a = F
    for b = F
        for c = F
            for d = F
                v = [a b c d];
                if any(v)
                    P = [P; v];
                end
            end
        end
    end
end
nPts = size(P,1);   % 15 points

%% =========================
%  LINES OF PG(3,2)
%% =========================

Lines = {};
for i = 1:nPts
    for j = i+1:nPts
        if rank([P(i,:); P(j,:)]) == 2
            L = unique(mod([P(i,:);
                            P(j,:);
                            P(i,:)+P(j,:)],2),'rows');
            if size(L,1) == 3
                found = false;
                for k = 1:length(Lines)
                    if isequal(sortrows(L), sortrows(Lines{k}))
                        found = true; break;
                    end
                end
                if ~found
                    Lines{end+1} = L;
                end
            end
        end
    end
end
nL = length(Lines);   % 35 lines

%% =========================
%  PLANES OF PG(3,2)
%% =========================

Planes = {};
for i = 1:nPts
    for j = i+1:nPts
        for k = j+1:nPts
            if rank([P(i,:); P(j,:); P(k,:)]) == 3
                B = [P(i,:); P(j,:); P(k,:)];
                Pl = [];
                for a = F
                    for b = F
                        for c = F
                            v = mod(a*B(1,:) + b*B(2,:) + c*B(3,:),2);
                            if any(v)
                                Pl = [Pl; v];
                            end
                        end
                    end
                end
                Pl = unique(Pl,'rows');
                if size(Pl,1) == 7
                    found = false;
                    for t = 1:length(Planes)
                        if isequal(sortrows(Pl), sortrows(Planes{t}))
                            found = true; break;
                        end
                    end
                    if ~found
                        Planes{end+1} = Pl;
                    end
                end
            end
        end
    end
end
nPl = length(Planes);   % 15 planes

%% =========================
%  INCIDENCE MATRICES
%% =========================

LP = zeros(nL, nPts);
for i = 1:nL
    for j = 1:nPts
        LP(i,j) = ismember(P(j,:), Lines{i}, 'rows');
    end
end

LPl = zeros(nL, nPl);
for i = 1:nL
    for j = 1:nPl
        LPl(i,j) = all(ismember(Lines{i}, Planes{j}, 'rows'));
    end
end

%% =========================
%  SEARCH PARAMETERS
%% =========================

target_size   = 16;    % <<< CORRECT TARGET
max_solutions = 5;

S = zeros(1,nL);
solutions = {};
SolutionsPointCounts = {};
SolutionsPlaneCounts = {};

stopSearch = false;

fprintf('Searching for at most %d solutions of size %d...\n', ...
        max_solutions, target_size);

backtrack(1,0);

%% =========================
%  WRITE RESULTS TO WORD FILE
%% =========================

fid = fopen('PG32_LineSets_16.doc','w');

fprintf(fid,'Line sets S in PG(3,2)\n');
fprintf(fid,'|S| = 16\n');
fprintf(fid,'Point condition: 3 or 4 lines through each point\n');
fprintf(fid,'Plane condition: 3 or 4 lines in each plane\n\n');

for s = 1:length(solutions)
    sol = solutions{s};

    fprintf(fid,'====================================\n');
    fprintf(fid,'Solution %d\n', s);
    fprintf(fid,'Line indices:\n%s\n\n', mat2str(sol));

    fprintf(fid,'Explicit lines:\n');
    for i = sol
        fprintf(fid,'%s\n', mat2str(Lines{i}));
    end

    fprintf(fid,'\nLines through each point:\n%s\n', ...
            mat2str(SolutionsPointCounts{s}));

    fprintf(fid,'Lines in each plane:\n%s\n\n', ...
            mat2str(SolutionsPlaneCounts{s}));
end

fclose(fid);

fprintf('Finished. %d solutions written to PG32_LineSets_16.doc\n', ...
        length(solutions));

%% =========================
%  NESTED BACKTRACK FUNCTION
%% =========================

function backtrack(idx, chosen)

    if stopSearch
        return
    end

    if chosen > target_size
        return
    end

    if idx > nL
        if chosen == target_size

            ptcount = S * LP;
            if any(ptcount < 3 | ptcount > 4)
                return
            end

            plcount = S * LPl;
            if any(plcount < 3 | plcount > 4)
                return
            end

            solutions{end+1} = find(S);
            SolutionsPointCounts{end+1} = ptcount;
            SolutionsPlaneCounts{end+1} = plcount;

            fprintf('Solution %d found\n', length(solutions));

            if length(solutions) >= max_solutions
                stopSearch = true;
            end
        end
        return
    end

    if chosen + (nL - idx + 1) < target_size
        return
    end

    S(idx) = 1;
    backtrack(idx+1, chosen+1);

    S(idx) = 0;
    backtrack(idx+1, chosen);
end

end
