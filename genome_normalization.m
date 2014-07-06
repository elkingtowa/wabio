%Read the data

M = importdata('12.15 caenorhabditis_elegans-scores.csv');
B = importdata('12.15 caenorhabditis_elegans-ids.csv');

% data referred to as A.data or B.data

% Normalize blast scores to C. elegans scores

for i = 1:length(M.data)
	M.data(i,:) = M.data(i,:)/M.data(i,7)
end;

% Add up rows

N = sum(M.data,2);

% Ranked top 200 rows using excel