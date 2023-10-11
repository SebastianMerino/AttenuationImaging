%%
function [diagSum] = diagSumCalc(squareMatrix, LLUR0_ULLR1)
%
% Input: squareMatrix: A square matrix.
%        LLUR0_ULLR1:  LowerLeft to UpperRight addition = 0
%                      UpperLeft to LowerRight addition = 1
%
% Output: diagSum: A vector of the sum of the diagnols of the matrix.
%
% Example:
%
% >> squareMatrix = [1 2 3;
%                    4 5 6;
%                    7 8 9];
%
% >> diagSum = diagSumCalc(squareMatrix, 0);
%
% diagSum =
%
%       1 6 15 14 9
%
% >> diagSum = diagSumCalc(squareMatrix, 1);
%
% diagSum =
%
%       7 12 15 8 3
%
% Written by M. Phillips
% Oct. 16th, 2013
% MIT Open Source Copywrite
% Contact mphillips@hmc.edu fmi.
%

if (nargin < 2)
    disp('Error on input. Needs two inputs.');
    return;
end

if (LLUR0_ULLR1 ~= 0 && LLUR0_ULLR1~= 1)
    disp('Error on input. Only accepts 0 or 1 as input for second condition.');
    return;
end

[M, N] = size(squareMatrix);

if (M ~= N)
    disp('Error on input. Only accepts a square matrix as input.');
    return;
end

diagSum = zeros(1, M+N-1);

if LLUR0_ULLR1 == 1
    squareMatrix = rot90(squareMatrix, -1);
end

for i = 1:length(diagSum)
    if i <= M
        countUp = 1;
        countDown = i;
        while countDown ~= 0
            diagSum(i) = squareMatrix(countUp, countDown) + diagSum(i);
            countUp = countUp+1;
            countDown = countDown-1;
        end
    end
    if i > M
        countUp = i-M+1;
        countDown = M;
        while countUp ~= M+1
            diagSum(i) = squareMatrix(countUp, countDown) + diagSum(i);
            countUp = countUp+1;
            countDown = countDown-1;
        end
    end
end

end
