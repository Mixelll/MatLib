function str = fit2str(fit)
str = [];
for i = fieldnames(fit)'
    str = [str i{:} '=' num2str(eval(['fit.' i{:}])) ' '];
end

end

