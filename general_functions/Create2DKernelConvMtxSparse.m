function [ mK ] = Create2DKernelConvMtxSparse( mH, numRows, numCols, convShape )
CONVOLUTION_SHAPE_FULL  = 1;
CONVOLUTION_SHAPE_SAME  = 2;
CONVOLUTION_SHAPE_VALID = 3;
numColsKernel   = size(mH, 2);
numBlockMtx     = numColsKernel;
cBlockMtx = cell(numBlockMtx, 1);
for ii = 1:numBlockMtx
    cBlockMtx{ii} = CreateConvMtxSparse(mH(:, ii), numRows, convShape);
end
switch(convShape)
    case(CONVOLUTION_SHAPE_FULL)
        diagIdx     = 0;
        numRowsKron = numCols + numColsKernel - 1;
    case(CONVOLUTION_SHAPE_SAME)
        diagIdx     = floor(numColsKernel / 2);
        numRowsKron = numCols;
    case(CONVOLUTION_SHAPE_VALID)
        diagIdx     = numColsKernel - 1;
        numRowsKron = numCols - numColsKernel + 1;
end
vI = ones(min(numRowsKron, numCols), 1);
mK = kron(spdiags(vI, diagIdx, numRowsKron, numCols), cBlockMtx{1});
for ii = 2:numBlockMtx
    diagIdx = diagIdx - 1;
    mK = mK + kron(spdiags(vI, diagIdx, numRowsKron, numCols), cBlockMtx{ii});
end
end