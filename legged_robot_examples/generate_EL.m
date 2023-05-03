clear;
clc;

syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 real
syms up1 up2 up3 up4 up5 up6 up7 up8 up9 up10 up11 up12 up13 up14 real
syms t k real


out = get_EL(t,u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, up1, up2, up3, up4, up5, up6, up7, up8, up9, up10, up11, up12, up13, up14, k, 1);
pLx = out(:,1);
pLxd = out(:,2);