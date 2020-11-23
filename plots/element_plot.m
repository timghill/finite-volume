function element_plot(dmesh,Cdata)
% element_plot is a thin wrapper to the matlab patch function. This wrapper
% makes it easier to work with this function, since its signature is quite
% complicated
patch('Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes,...
        'FaceVertexCData', Cdata, 'FaceColor', 'flat');