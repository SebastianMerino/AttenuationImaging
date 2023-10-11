
function save_all_figures_to_directory(dir_name)

figlist=findobj('type','figure');
number=get(figlist,'Number');
for i=1:numel(figlist)
    %saveas(figlist(i),fullfile(dir_name,['figure' num2str(figlist(i)) '.fig']));
    figure(figlist(i))
    set(gcf,'PaperPositionMode','auto')
    %pause(2)
    %saveas(figlist(i),fullfile(dir_name,['figure' num2str(figlist(i)) '.png']));
    %pause(2)
    saveas(figlist(i),fullfile(dir_name,['figure' char(string(number(i))) '.png']));
    saveas(figlist(i),fullfile(dir_name,['figure' char(string(number(i))) '.fig']));
    
end

end