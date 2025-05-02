function [RGB] = map_to_rgb(x)
% x 是一个实数，它的范围是 [0, 1]。
% R、G、B 是从 x 映射到 RGB 颜色空间的红、绿、蓝分量。
color1 = [0 0 255];
color2 = [255,0,0];
% 将 x 映射到 RGB 空间的红、绿、蓝分量
R = (color1(1) + x * (color2(1)-color1(1)))/255;
G = -4*(x-0.5)^2 + 1;
B = (color1(3) + x * (color2(3)-color1(3)))/255;
RGB = [R G B];
end
