%https://fr.mathworks.com/matlabcentral/answers/439928-creating-convolution-matrix-of-2d-kernel-for-different-shapes-of-convolution
function [ mK ] = CreateImageConvMtx( mH, nRows, nCols, convShape )
%convShape is 'full', 'same', or 'valid'
   impulse=zeros(nRows,nCols);
   
   for i=numel(impulse):-1:1
  
         impulse(i)=1;  %Create impulse image corresponding to i-th output matrix column
         
         tmp=sparse( conv2(impulse,mH,convShape) );  %impulse response
   
         Column{i}=tmp(:);
         
         impulse(i)=0;
   end
   mK=cell2mat(Column);   
end