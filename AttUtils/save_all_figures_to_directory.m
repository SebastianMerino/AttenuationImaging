
function save_all_figures_to_directory(dir_name,title, varargin)

figlist=findobj('type','figure');
number=get(figlist,'Number');
for i=1:numel(figlist)
    %saveas(figlist(i),fullfile(dir_name,['figure' num2str(figlist(i)) '.fig']));
    %figure(figlist(i))
    %set(gcf,'PaperPositionMode','auto')
    %pause(2)
    %saveas(figlist(i),fullfile(dir_name,['figure' num2str(figlist(i)) '.png']));
    %pause(2)
    saveas(figlist(i),fullfile(dir_name,[title char(string(number(i))) '.png']));
    if nargin == 2
        saveas(figlist(i),fullfile(dir_name,[title char(string(number(i))) '.fig']));
    elseif varargin{1}
        saveas(figlist(i),fullfile(dir_name,[title char(string(number(i))) '.fig']));
    end
end

end