function mat = np2mat(npary)

% convert python (Numpy) ndarray to matlab
sh = double(py.array.array('d',npary.shape));
if numel(sh) == 1
    sh = [sh 1];
end
npary2 = double(py.array.array('d',py.numpy.nditer(npary)));
mat = reshape(npary2,fliplr(sh))';  % matlab 2d array 

end