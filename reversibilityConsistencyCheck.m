%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reversibilityConsistencyCheck(model)
%
% Function that checks the consistency between the reversibility and lower 
% bounds fields of a GEM matlab structure. 
%
% Ivan Domenzain.   Last edited: 2017-12-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inconsistencies = reversibilityConsistencyCheck(model)
    inconsistencies = [];
    for i=1:length(model.rev)
        if model.rev(i) == 1 && model.lb(i) == 0
          disp(['WARNING:Rxn ' num2str(i) ' reversibility insconsistency'])
          inconsistencies = [inconsistencies, i];
        else
          disp(['Rxn ' num2str(i) ' reversibility consistency --> OK'])
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%