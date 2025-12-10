function [m_out] = first_time(m_in)
    m_out = NaN([size(m_in,1) size(m_in,2)]);

    % for i=1:size(m_in,1)
    %     for j= 1:size(m_in,2)
    %         time = find(m_in(i,j,:) > 0, 1,'first');
    %         if ~isempty(time)
    %             m_out(i,j) = time;
    %         end
    %     end
    % end

    [vals m_out] = max(m_in>0, [], 3);
end