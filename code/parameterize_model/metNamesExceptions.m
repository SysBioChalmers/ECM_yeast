function newMet = metNamesExceptions(met,EC)
    newMet = met;
    if strcmpi(EC,'EC2.7.1.6')
        if strcmpi(met,'galactose')
            newMet = {'d-galactose'};
        end

    elseif strcmpi(EC,'EC5.1.3.3')
        if strcmpi(met,'UDP-D-galactose')
            newMet = {'d-galactose'};
        elseif strcmpi(met,'UDP-glucose')
            newMet = {'d-glucose'};
        end

    elseif strcmpi(EC,'EC2.7.7.12')
        if strcmpi(met,'D-galactose 1-phosphate')
            newMet = {'alpha-d-galactose 1-phosphate'};
        elseif strcmpi(met,'UDP-D-galactose')
            newMet = {'udp-alpha-d-galactose'};
        elseif strcmpi(met,'D-glucose 1-phosphate')
            newMet = {'alpha-d-glucose 1-phosphate'};
        end

    elseif strcmpi(EC,'EC5.4.2.2')
        if strcmpi(met,'alpha-D-glucose 6-phosphate')
            newMet = {'d-glucose 6-phosphate'};
        end

    elseif strcmpi(EC,'EC2.7.1.1')
        if strcmpi(met,'alpha-D-glucose')
            newMet = {'d-glucose'};
        end

    elseif strcmpi(EC,'EC5.3.1.9')
        if strcmpi(met,'alpha-D-glucose 6-phosphate')
            newMet = {'d-glucose 6-phosphate'};
        elseif strcmpi(met,'beta-D-fructofuranose 6-phosphate')
            newMet = {'d-fructose 6-phosphate'};
        end

    elseif strcmpi(EC,'EC2.7.1.11')
        if strcmpi(met,'beta-D-fructofuranose 6-phosphate')
            newMet = {'d-fructose 6-phosphate'};
        elseif strcmpi(met,'beta-D-fructofuranose 1,6-bisphosphate')
            newMet = {'d-fructose 1,6-bisphosphate'};
        end

    elseif strcmpi(EC,'EC4.1.2.13')
        if strcmpi(met,'beta-D-fructofuranose 1,6-bisphosphate')
            newMet = {'d-fructose 1,6-bisphosphate'};
        end
    
    elseif strcmpi(EC,'EC5.4.2.11')
        if strcmpi(met,'2-phospho-D-glycerate')
            newMet = {'2-phosphoglycerate'};
        end

    elseif strcmpi(EC,'EC1.1.1.49')
        if strcmpi(met,'alpha-D-glucose 6-phosphate')
            newMet = {'d-glucose 6-phosphate'};
        end

    elseif strcmpi(EC,'EC1.1.1.44')
        if strcmpi(met,'D-ribulose 5-phosphate')
            newMet = {'ribulose-5-phosphate'};
        end

    elseif strcmpi(EC,'EC1.2.1.12')
        if strcmpi(met,'1,3-bisphospho-D-glycerate')
            newMet = {'1,3-diphosphoglyceric acid'};
        end

    elseif strcmpi(EC,'EC2.2.1.2')
        if strcmpi(met,'D-glyceraldehyde 3-phosphate')
            newMet = {'glyceraldehyde 3-phosphate'};
        elseif strcmpi(met,'beta-D-fructofuranose 6-phosphate')
            newMet = {'d-fructose 6-phosphate'};
        end

    elseif strcmpi(EC,'EC3.1.3.21')
        if strcmpi(met,'glycerol-3-phosphate')
            newMet = {'glycerol 3-phosphate'};
        end

    elseif strcmpi(EC,'EC4.1.3.30')
        if strcmpi(met,'methylisocitrate')
            newMet = {'isocitrate'};
        end

    elseif strcmpi(EC,'EC6.2.1.1')
        if strcmpi(met,'coenzyme A')
            newMet = {'coa'};
        end

    elseif strcmpi(EC,'EC6.2.1.5')
        if strcmpi(met,'coenzyme A')
            newMet = {'coa'};
        end

    elseif strcmpi(EC,'EC1.3.1.6')
        if strcmpi(met,'FADH2')
            newMet = {'fadh2'};
        end

    end
end
