
function save_all_figures_to_directory(dir_name,title, varargin)
% optional argument is other format to save the image

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
    if nargin == 3
        saveas(figlist(i),fullfile(dir_name, ...
            [title,char(string(number(i))),'.',varargin{1}]));
    end
end

end