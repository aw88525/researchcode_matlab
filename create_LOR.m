function [LOR_blur] = create_LOR(LOR_ID_default, crystal_pos_1_ID, crystal_pos_2_ID, N_crystal, event_num)
    if(length(crystal_pos_1_ID) == 0 || length(crystal_pos_2_ID) == 0)
        LOR_blur = {LOR_ID_default, 0, -1, -1};
        return;
    end
    if(length(crystal_pos_1_ID) == 1)
        if(crystal_pos_1_ID == 0)
           LOR_blur = {LOR_ID_default, 0, -1, -1};
           return;
        end
        
    end
    if(length(crystal_pos_2_ID) == 1)
        if(crystal_pos_2_ID == 0)
           LOR_blur = {LOR_ID_default, 0, -1, -1};
           return;
        end
        
    end
    crystal_pos_get_rid_of_zero = crystal_pos_1_ID == 0 | crystal_pos_2_ID == 0;
    crystal_pos_1_ID_new = crystal_pos_1_ID(~crystal_pos_get_rid_of_zero);
    crystal_pos_2_ID_new = crystal_pos_2_ID(~crystal_pos_get_rid_of_zero);
    crystal_pos_1_new = [crystal_pos_1_ID_new(crystal_pos_1_ID_new < crystal_pos_2_ID_new)-1 crystal_pos_2_ID_new(crystal_pos_2_ID_new < crystal_pos_1_ID_new)-1];
    crystal_pos_2_new = [crystal_pos_2_ID_new(crystal_pos_1_ID_new < crystal_pos_2_ID_new)-1 crystal_pos_1_ID_new(crystal_pos_2_ID_new < crystal_pos_1_ID_new)-1];
%    if(sum(ismember(crystal_MEM,crystal_ID_1)  == 1))
%      figure;
%      crystal_ID_s = [crystal_pos_1_new' crystal_pos_2_new'];
%      hist3(crystal_ID_s)
%    end
    %     crystal_f = (N_crystal-1).*ones(length(crystal_pos_1_new),1);
%     incremental_f = {crystal_pos_1_new(1,:):1:crystal_f(1,:)'};
%     
%     LOR_ID = sum(crystal_pos_1_new:1:crystal_f', 1);
    
    %index_crystal_pos_1_odd = mod(crystal_pos_1_new,2) == 1;
    %index_crystal_pos_1_even = mod(crystal_pos_1_new,2) == 0 & crystal_pos_1_new > 0;
    
%     lookup_table_new = [lookup_table(:,2) lookup_table(:,3)];
%     crystal_pos_new = [crystal_pos_1_new', crystal_pos_2_new'];
   
   



    LOR_ID(crystal_pos_1_new>0) = N_crystal.*crystal_pos_1_new(crystal_pos_1_new>0)-(1+crystal_pos_1_new(crystal_pos_1_new>0)).*crystal_pos_1_new(crystal_pos_1_new>0)./2;
    LOR_ID(crystal_pos_1_new == 0) = 0;
    LOR_ID = LOR_ID + crystal_pos_2_new-crystal_pos_1_new-1; % will have to change %this is due to an error in MC code
%     [loa,LOR_ID] = ismember(crystal_pos_new,lookup_table_new,'rows');
%     LOR_ID = lookup_table(LOR_ID_number,1);

    %LOR_ID(index_crystal_pos_1_odd) = (N_crystal-crystal_pos_1_new(index_crystal_pos_1_odd)-1+N_crystal-1)./2.*(crystal_pos_1_new(index_crystal_pos_1_odd)+1);
    %LOR_ID(index_crystal_pos_1_even) = (N_crystal-crystal_pos_1_new(index_crystal_pos_1_even)-1+N_crystal-1).*(crystal_pos_1_new(index_crystal_pos_1_even)+1)./2;
    
    LOR_ID_uni = unique(LOR_ID);
    counts = histc(LOR_ID, LOR_ID_uni);
%     [counts,LOR_ID_uni] = hist(LOR_ID,unique(LOR_ID));
%      counts = LOR_ID(ismember(LOR_ID, unique(LOR_ID))); 
     
     LOR_ID_uni = int32(LOR_ID_uni);
%      LOR_prob = counts./length(crystal_pos_1_ID_new);
     LOR_prob = counts./event_num;
     
     if(size(LOR_ID_uni,2) == 0)
         LOR_blur = {LOR_ID_default size(LOR_ID_uni,2) -1 -1};
     else
%          LOR_prob_save = 1/sum(LOR_prob(LOR_prob>0.01*max(LOR_prob))).*LOR_prob(LOR_prob>0.01*max(LOR_prob));
%          LOR_prob_save = sum(LOR_prob)./sum(LOR_prob(LOR_prob>0.01*max(LOR_prob))).*LOR_prob(LOR_prob>0.01*max(LOR_prob));
%          LOR_ID_uni_save = LOR_ID_uni(LOR_prob>0.01*max(LOR_prob));
         LOR_prob_save = LOR_prob;
         LOR_ID_uni_save = LOR_ID_uni;
         LOR_blur = {LOR_ID_default size(LOR_ID_uni_save,2) LOR_ID_uni_save LOR_prob_save};
     end
end