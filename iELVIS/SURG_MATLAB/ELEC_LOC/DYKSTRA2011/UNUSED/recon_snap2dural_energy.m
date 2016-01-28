function coord_snapped=recon_snap2dural_energy(coord,surf,pairs)

% get starting coordinates
coord0=coord+repmat([5 0 0],size(coord,1),1);

ind=zeros(size(coord0,1),1);
for c=1:size(coord0,1)
    [d,ind(c)]=min(sum((repmat(coord0(c,:),size(surf.vert,1),1)-surf.vert).^2,2));
end
% coord0=surf.vert(ind,:);
% coord0=coord+repmat([10 0 -5],size(coord,1),1);
coord0=coord;
% coord0(30,:)=[0 0 0];

% compute pairs of neighbors
if(nargin<3 || isempty(pairs))
    pairs=knn_pairs(coord,4);
end

% anonymous function handles
efun=@(coord_snapped)energy_electrodesnap(coord_snapped,coord,pairs);
cfun=@(coord_snapped)surface_constraint(coord_snapped,surf);

% options
options=optimset('Algorithm','active-set',...
                 'Display','iter',...
                 'MaxIter',50,...
                 'MaxFunEvals',Inf,...
                 'UseParallel','always',...
                 'GradObj','off',...
                 'TypicalX',coord,...
                 'DiffMaxChange',2,...
                 'DiffMinChange',0.3,...
                 'TolFun',0.3,...
                 'TolCon',0.01*size(coord0,1),...
                 'TolX',0.5,...
                 'Diagnostics','off',...
                 'RelLineSrchBnd',1,...
                 'PlotFcns',{@optimplotfval,@optimplotstepsize,@optimplotfirstorderopt,@plotCoordFun});
             
% history
history.coords=[];
last_handle=[];
             
% run minimization
coord_snapped = fmincon(efun,coord0,[],[],[],[],[],[],cfun,options);

% patternsearch minimization
% psoptions=psoptimset(options,'PlotFcn',{@plotCoordFunPS,@psplotbestf,@psplotmeshsize});
% 
% optimproblem.objective=efun;
% optimproblem.x0=coord0;
% optimproblem.Aineq=[];
% optimproblem.bineq=[];
% optimproblem.Aeq=[];
% optimproblem.beq=[];
% optimproblem.lb=[];
% optimproblem.up=[];
% optimproblem.nonlcon=cfun;
% optimproblem.solver='patternsearch';
% optimproblem.options=psoptions;
% 
% coord_snapped=patternsearch(optimproblem);

% other functions
    function stop = plotCoordFun(coord,optimValues,state)

        switch state
            case 'init'
                hold on; axis vis3d; view([3 1 -0.5]);
                last_handle=scatter3(coord(:,1),coord(:,2),coord(:,3),'marker','o','markerfacecolor','r','markeredgecolor','r','sizedata',15);
                history.coords=coord;
            case 'iter'
                set(last_handle,'marker','.','sizedata',10,'markerfacecolor','b','markeredgecolor','b');
                last_handle=scatter3(coord(:,1),coord(:,2),coord(:,3),10,'marker','o','markerfacecolor','r','markeredgecolor','r','sizedata',15);
                plot3([coord(:,1) history.coords(:,1,end)]',[coord(:,2) history.coords(:,2,end)]',[coord(:,3) history.coords(:,3,end)]','k');
                history.coords=cat(3,history.coords,coord);
            case 'done'
                hold off;
            otherwise
        end
        stop = false;
    end

    function stop = plotCoordFunPS(optimValues,state)
        coord=optimValues.x;
        stop=plotCoordFun(coord,optimValues,state);
    end
end