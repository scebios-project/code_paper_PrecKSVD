function ndarray = toNdarray(M)

[rows, cols] = size(M);
ndarray = py.numpy.array(M(:));
ndarray = py.numpy.reshape(ndarray, [uint32(rows), uint32(cols)], 'F');