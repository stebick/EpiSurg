function [verticeslist]=SearchProjectionOnPial(mesh_total,mesh_outer,perim, step)

% step is typically set between 5 and 10. Default is 7, increasing it will
% limit redundancies in the resulting path file.

verticeslist=[];
<<<<<<< HEAD
for t=1:step:size(perim,2)
=======
si=max(size(perim));
for t=1:step:si
>>>>>>> epiSurg/master
    [nearestIndexMT,nearestValuesMT]=mesh_vertex_nearest(mesh_total.vertices,mesh_outer.vertices(perim(t),:));
    verticeslist= [verticeslist nearestIndexMT];
end
verticeslist=unique(verticeslist);
