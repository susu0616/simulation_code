function D = BregDiv(para1, para2, type)
sigma1 = 10;
sigma2 = 1;
if type ==1
    D = sigma1*(para1-para2)'*(para1-para2);
elseif type ==2
    D = sigma2*(para1-para2)'*(para1-para2);
end
end