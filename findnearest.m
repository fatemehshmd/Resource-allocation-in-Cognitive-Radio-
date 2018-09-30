function [su_num] = findnearest_fati(p_su_x, p_su_y,p_su_x_new,p_su_y_new)


dis = sqrt((p_su_x - p_su_x_new).^2+ (p_su_y - p_su_y_new).^2);

su_num = find(dis == min(dis));


end